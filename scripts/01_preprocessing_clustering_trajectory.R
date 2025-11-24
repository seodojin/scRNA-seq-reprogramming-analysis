# ============================================================================
# Complete pipeline from raw data to trajectory inference
#              (Preprocessing → ScType annotation → Slingshot trajectory)
# 
# NOTE: This script requires raw count matrices from NCBI SRA BioProject PRJNA1256192
# Raw data files (not included):
#   - 1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv
#   - 1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(HGNChelper)
})

SEED <- 123
set.seed(SEED)

# ============================================================================
# PART 1: Data Loading and Preprocessing
# ============================================================================

message("=== PART 1: Data Loading and Preprocessing ===\n")

# Load control sample
message("Loading shCtrl sample...")
counts <- read.table("1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv",
                     skip = 5, sep = ",", header = TRUE, row.names = 1)
shCtrl <- CreateSeuratObject(counts = t(counts), project = "shCtrl")
message("  shCtrl: ", ncol(shCtrl), " cells")

# Load PTBP1 knockdown sample
message("Loading shPTBP1 sample...")
counts_sh <- read.table("1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv",
                        skip = 5, sep = ",", header = TRUE, row.names = 1)
shPTBP1 <- CreateSeuratObject(counts = t(counts_sh), project = "shPTBP1")
message("  shPTBP1: ", ncol(shPTBP1), " cells")

# Merge datasets
plus <- merge(shCtrl, shPTBP1, add.cell.ids = c("shCtrl", "shPTBP1"), project = "both")
message("  Merged: ", ncol(plus), " cells")

# ============================================================================
# Quality Control
# ============================================================================

message("\n--- Quality Control ---")

# Calculate mitochondrial percentage
plus[["percent.mt"]] <- PercentageFeatureSet(plus, pattern = "^MT.")

# Create QC violin plots
p_qc <- VlnPlot(plus,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3) & 
  labs(x = NULL) &
  theme(legend.position = "none")

ggsave("figures/QC_ViolinPlots.svg", p_qc, width = 12, height = 4)
message("  ✓ Saved: figures/QC_ViolinPlots.svg")

# Filter cells based on QC metrics
plus <- subset(plus, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 20)
message("  After QC filtering: ", ncol(plus), " cells retained")

# ============================================================================
# Normalization and Dimensionality Reduction
# ============================================================================

message("\n--- Normalization and Dimensionality Reduction ---")

plus <- NormalizeData(plus, normalization.method = "LogNormalize", 
                      scale.factor = 10000, verbose = FALSE)
message("  ✓ Normalization completed")

plus <- FindVariableFeatures(plus, selection.method = "vst", 
                             nfeatures = 2000, verbose = FALSE)
message("  ✓ Found ", length(VariableFeatures(plus)), " variable features")

plus <- ScaleData(plus, features = VariableFeatures(plus), verbose = FALSE)
message("  ✓ Data scaled")

plus <- RunPCA(plus, features = VariableFeatures(plus), verbose = FALSE)
message("  ✓ PCA completed")

# ============================================================================
# Clustering
# ============================================================================

message("\n--- Clustering ---")

plus <- FindNeighbors(plus, dims = 1:10, verbose = FALSE)
plus <- FindClusters(plus, resolution = 0.2, verbose = FALSE)
message("  ✓ Found ", length(unique(plus$seurat_clusters)), " clusters")

plus <- RunUMAP(plus, dims = 1:10, n.neighbors = 15, min.dist = 0.3, verbose = FALSE)
message("  ✓ UMAP completed")

# Visualize
p_umap <- DimPlot(plus, reduction = "umap", split.by = "orig.ident") +
  ggtitle("UMAP by Condition")
ggsave("figures/UMAP_by_condition.svg", p_umap, width = 10, height = 5)
message("  ✓ Saved: figures/UMAP_by_condition.svg")

# ============================================================================
# PART 2: ScType Cell Type Annotation
# ============================================================================

message("\n=== PART 2: ScType Cell Type Annotation ===\n")

# Load ScType functions
source("sctype_function_1.R")
source("sctype_function_2.R")

# ScType database
db_file <- "ScTypeDB_full.xlsx"
tissue <- "Brain"

message("Loading gene sets from: ", db_file)
gs_list <- gene_sets_prepare(db_file, tissue)

# Get scaled data
scale_data <- GetAssayData(plus, layer = "scale.data")

# Calculate ScType scores
message("Calculating ScType scores...")
es.max <- sctype_score(
  scRNAseqData = scale_data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)

# Calculate scores per cluster
cL_results <- do.call("rbind", lapply(unique(plus@meta.data$seurat_clusters), function(cl) {
  es.max.cl <- sort(
    rowSums(es.max[, rownames(plus@meta.data[plus@meta.data$seurat_clusters == cl, ])]),
    decreasing = TRUE
  )
  head(data.frame(
    cluster = cl,
    type = names(es.max.cl),
    scores = es.max.cl,
    ncells = sum(plus@meta.data$seurat_clusters == cl)
  ), 10)
}))

# Select top cell type per cluster
sctype_scores <- cL_results %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)

# Set low-confidence clusters to "Unknown"
sctype_scores$type[as.numeric(sctype_scores$scores) < sctype_scores$ncells / 4] <- "Unknown"

message("\nScType results:")
print(sctype_scores[, 1:3])

# Assign cell types to cells
plus@meta.data$customclassif <- ""
for(j in unique(sctype_scores$cluster)){
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  plus@meta.data$customclassif[plus@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# Replace empty with "Unknown"
plus@meta.data$customclassif[plus@meta.data$customclassif == ""] <- "Unknown"

# ============================================================================
# Filter and merge cell types for visualization
# ============================================================================

message("\n--- Filtering and merging cell types for Figure 2 ---")

# Remove Unknown cells
message("Before filtering - Total cells: ", ncol(plus))
message("Unknown cells: ", sum(plus@meta.data$customclassif == "Unknown"))

plus <- subset(plus, subset = customclassif != "Unknown")
message("After filtering Unknown - Remaining cells: ", ncol(plus))

# Merge GABAergic and Glutamatergic neurons into unified "Neurons" category
plus@meta.data$customclassif <- ifelse(
  plus@meta.data$customclassif %in% c("GABAergic neurons", "Glutamatergic neurons"),
  "Neurons",
  plus@meta.data$customclassif
)

message("\nCell type composition after merging:")
print(table(plus@meta.data$customclassif, plus@meta.data$orig.ident))

# ============================================================================
# FIGURE 2 GENERATION
# ============================================================================

message("\n=== GENERATING FIGURE 2 PANELS (a), (b), (c) ===\n")

# Define consistent colors
cell_colors <- c(
  "Fibroblasts" = "#19c3a3",
  "Immature neurons" = "#00bef3",
  "Neurons" = "#d4a600",
  "Myofibroblasts" = "#ff8c8c"
)

# Panel (a): Split UMAP by condition
message("Generating Figure 2a: Split UMAP...")
p_fig2a <- DimPlot(plus, 
                   reduction = "umap",
                   group.by = "customclassif",
                   split.by = "orig.ident",
                   cols = cell_colors,
                   pt.size = 0.5) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))

ggsave("figures/UMAP_split.png", p_fig2a, width = 12, height = 5, dpi = 600)

# Panel (c): Cell composition
message("Generating Figure 2c: Cell Composition...")
composition_data <- plus@meta.data %>%
  group_by(orig.ident, customclassif) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count) * 100)

composition_data$customclassif <- factor(composition_data$customclassif,
                                         levels = c("Neurons", "Immature neurons", 
                                                    "Myofibroblasts", "Fibroblasts"))

p_fig2c <- ggplot(composition_data, aes(x = orig.ident, y = proportion, 
                                        fill = customclassif)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))

print(p_fig2c)
ggsave("figures/cell_type_composition.png", p_fig2c, width = 6, height = 5, dpi = 600)

# Panel (b): PTBP1 Expression
message("Generating Figure 2b: PTBP1 Expression...")

# Filter for shPTBP1 cells only
shPTBP1_cells <- rownames(plus@meta.data[plus@meta.data$orig.ident == "shPTBP1", ])

# Verify cell type distribution for shPTBP1 cells
message("Cell type distribution in shPTBP1:")
print(table(plus@meta.data[shPTBP1_cells, "customclassif"]))

# Set active identity (using already-merged customclassif)
Idents(plus) <- plus@meta.data$customclassif


# Generate violin plot
vln_plot2 <- VlnPlot(
  object = subset(plus, cells = shPTBP1_cells),
  features = "PTBP1",
  group.by = "customclassif",
  pt.size = 0,
  cols = c("#19c3a3", "#ff8c8c", "#00bef3", "#d4a600")
) +
  ggtitle("PTBP1 Expression in shPTBP1") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(vln_plot2)
ggsave("figures/ptbp1_vlnplot_shPTBP1.png", plot = vln_plot2, width = 7, height = 6, dpi = 600)

# ============================================================================
# PART 3: Slingshot Trajectory Inference (shPTBP1 only)
# ============================================================================

message("\n=== PART 3: Slingshot Trajectory Inference ===\n")
message("Rationale: Neuronal differentiation occurs only in shPTBP1 group\n")

# Filter for shPTBP1 cells only
plus_shPTBP1 <- subset(plus, subset = orig.ident == "shPTBP1")
message("Total shPTBP1 cells: ", ncol(plus_shPTBP1))

# Remove "Unknown" cells for cleaner trajectory
plus_filtered <- subset(plus_shPTBP1, subset = customclassif != "Unknown")
message("Cells after filtering Unknown: ", ncol(plus_filtered))

message("\nCell type composition:")
print(table(plus_filtered$customclassif))

# Save preprocessed object for downstream analysis
saveRDS(plus_filtered, "data/seurat_preprocessed.rds")

# ============================================================================
# Convert to SingleCellExperiment
# ============================================================================

message("\n--- Preparing SingleCellExperiment object ---")

# Extract and merge count matrices from both assay slots
counts_shCtrl <- plus[["RNA"]]$counts.shCtrl
counts_shPTBP1 <- plus[["RNA"]]$counts.shPTBP1
counts <- cbind(counts_shCtrl, counts_shPTBP1)

# Extract and merge normalized data
data_shCtrl <- plus[["RNA"]]$data.shCtrl
data_shPTBP1 <- plus[["RNA"]]$data.shPTBP1
data <- cbind(data_shCtrl, data_shPTBP1)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    logcounts = data
  ),
  colData = plus@meta.data,
  reducedDims = list(
    PCA = Embeddings(plus, "pca"),
    UMAP = Embeddings(plus, "umap")
  )
)

# Add active identity
sce$cell_type <- plus@active.ident

# Perform normalization and PCA
sce <- logNormCounts(sce)
sce <- scater::runPCA(sce, exprs_values = "logcounts")

# Update cluster labels: merge neuron subtypes
update_cluster_labels <- function(label) {
  if (label %in% c("GABAergic neurons", "Glutamatergic neurons")) {
    return("Neurons")
  } else {
    return(label)
  }
}

sce$updated_customclassif <- sapply(sce$customclassif, update_cluster_labels)

# Verify unique clusters
unique_clusters <- unique(sce$updated_customclassif)
print(unique_clusters)

# ============================================================================
# Prepare filtered SCE for trajectory analysis
# ============================================================================

# Convert filtered Seurat object to SCE
sce_full <- as.SingleCellExperiment(plus_filtered)

# Add reduced dimensions
reducedDim(sce_full, "UMAP") <- Embeddings(plus_filtered, "umap")
reducedDim(sce_full, "PCA") <- Embeddings(plus_filtered, "pca")

# Merge neuron subtypes
message("Merging GABAergic and Glutamatergic neurons into 'Neurons'...")
sce_full$celltype_merged <- ifelse(
  sce_full$customclassif %in% c("GABAergic neurons", "Glutamatergic neurons"),
  "Neurons",
  sce_full$customclassif
)

message("\nCell type composition after merging:")
print(table(sce_full$celltype_merged))

# ============================================================================
# Gene and cell filtering
# ============================================================================

message("\n--- Filtering for trajectory inference ---")

# Gene filtering (memory optimization)
counts_sparse <- as(assay(sce_full, "counts"), "sparseMatrix")
gene_filter <- rowSums(counts_sparse >= 10) >= 3
sce_full <- sce_full[gene_filter, ]
message("  Genes after filtering: ", nrow(sce_full))

# Cell filtering
cell_filter <- colSums(assay(sce_full, "counts") > 0) >= 500
sce_full <- sce_full[, cell_filter]
message("  Cells after filtering: ", ncol(sce_full))

# ============================================================================
# Run Slingshot
# ============================================================================

message("\n--- Running Slingshot trajectory inference ---")

sds <- slingshot(
  sce_full,
  clusterLabels = 'celltype_merged',
  reducedDim = 'UMAP',
  start.clus = "Fibroblasts",
  end.clus = c("Neurons", "Immature neurons", "Myofibroblasts"),
  allow.breaks = FALSE,
  extend = 'n',
  approx_points = 300
)

message("✓ Slingshot completed!")

# ============================================================================
# Panel (d): UMAP with trajectory lines
# ============================================================================

# Extract UMAP coordinates
umap_coords <- reducedDims(sce)$UMAP

# Create plot data frame
plot_data <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2],
  CellType = sce$updated_customclassif
)

# Define trajectory colors
trajectory_colors <- c("#21908c", "#fde725", "#440154")

# Initialize plot
p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cell_colors) +
  theme_minimal() +
  labs(title = "Cell Types with Trajectories")

# Add trajectory curves
for (i in seq_along(slingCurves(sds))) {
  curve_data <- slingCurves(sds)[[i]]$s[slingCurves(sds)[[i]]$ord, ]
  
  # Find terminal cluster for this lineage
  end_cluster <- slingLineages(sds)[[i]][length(slingLineages(sds)[[i]])]
  
  # Find centroid of terminal cluster
  end_cluster_cells <- which(sce$updated_customclassif == end_cluster)
  end_point <- colMeans(umap_coords[end_cluster_cells,])
  
  # Calculate distances to terminal point
  distances <- sqrt(rowSums((curve_data - matrix(end_point, nrow = nrow(curve_data), ncol = 2, byrow = TRUE))^2))
  
  # Find closest point to terminal cluster
  closest_point <- which.min(distances)
  
  # Plot trajectory up to terminal cluster
  p <- p + geom_path(data = data.frame(UMAP1 = curve_data[1:closest_point,1],
                                       UMAP2 = curve_data[1:closest_point,2]),
                     aes(x = UMAP1, y = UMAP2),
                     color = trajectory_colors[i],
                     linewidth = 1,
                     alpha = 0.7)
}

print(p)
ggsave("figures/trajectory_umap_plot.png", plot = p, width = 5, height = 4, dpi = 1000)

# ============================================================================
# Extract pseudotime and weights
# ============================================================================

pseudotime_all <- slingPseudotime(sds)
weights_all <- slingCurveWeights(sds)

n_lineages <- ncol(pseudotime_all)
message("\nNumber of lineages detected: ", n_lineages)

# Check lineage coverage
for (i in 1:n_lineages) {
  n_valid <- sum(!is.na(pseudotime_all[, i]))
  pct_valid <- round(n_valid / nrow(pseudotime_all) * 100, 1)
  message(sprintf("  Lineage %d: %d cells (%.1f%%)", i, n_valid, pct_valid))
}

# ============================================================================
# Primary lineage selection
# ============================================================================

message("\n--- Applying primary lineage selection ---")
message("Strategy: Assign each cell to its highest-weight lineage\n")

# Find primary lineage for each cell
primary_lineage <- apply(weights_all, 1, function(w) {
  if (all(is.na(w))) return(NA)
  which.max(w)
})

message("Cells assigned to each lineage:")
for (i in 1:n_lineages) {
  n_cells <- sum(primary_lineage == i, na.rm = TRUE)
  pct <- round(n_cells / length(primary_lineage) * 100, 1)
  message(sprintf("  Primary lineage %d: %d cells (%.1f%%)", i, n_cells, pct))
}

# Extract pseudotime from primary lineage
pseudotime_primary <- sapply(1:nrow(pseudotime_all), function(i) {
  lin <- primary_lineage[i]
  if (is.na(lin)) return(NA)
  pseudotime_all[i, lin]
})

# Check for NAs
n_na <- sum(is.na(pseudotime_primary))
n_valid <- sum(!is.na(pseudotime_primary))

message("\nAfter primary lineage selection:")
message(sprintf("  Valid cells: %d (%.1f%%)", n_valid, 
                round(n_valid/length(pseudotime_primary)*100, 1)))
message("  NA cells: ", n_na)

# Keep only valid cells
valid_cells <- !is.na(pseudotime_primary)
pseudotime_final <- pseudotime_primary[valid_cells]
primary_lineage_final <- primary_lineage[valid_cells]
sce_final <- sce_full[, valid_cells]

message("\n=== Final dataset for downstream analysis ===")
message("Cells: ", ncol(sce_final))
message("Genes: ", nrow(sce_final))

# ============================================================================
# Panel (e): Normalized pseudotime distribution
# ============================================================================

message("Generating Figure 2e: Pseudotime Distribution...")

# Normalize pseudotime
normalize_pseudotime <- function(pseudotime) {
  apply(pseudotime, 2, function(x) {
    valid_idx <- !is.na(x)
    if(sum(valid_idx) > 0) {
      x[valid_idx] <- (x[valid_idx] - min(x[valid_idx])) / 
        (max(x[valid_idx]) - min(x[valid_idx]))
    }
    return(x)
  })
}

normalized_pseudotime <- normalize_pseudotime(pseudotime_all)

# Assign cells to lineages
cell_assignments <- apply(weights_all, 1, which.max)

# Prepare data
plot_data <- data.frame(
  pseudotime = normalized_pseudotime[cbind(1:nrow(normalized_pseudotime), 
                                           cell_assignments)],
  lineage = cell_assignments
)

plot_data <- plot_data[!is.na(plot_data$pseudotime), ]

# Map lineage numbers to cell types
lineage_mapping <- c(
  "1" = "Neurons",
  "2" = "Myofibroblasts", 
  "3" = "Immature neurons"
)

plot_data$lineage_label <- factor(
  lineage_mapping[as.character(plot_data$lineage)],
  levels = c("Immature neurons", "Myofibroblasts", "Neurons")
)

# Colors for lineages
lineage_colors <- c(
  "Neurons" = "#440154",
  "Myofibroblasts" = "#21908c",
  "Immature neurons" = "#fde725"
)

# Create violin plot
p_fig2e <- ggplot(plot_data, aes(x = pseudotime, y = lineage_label, 
                                 fill = lineage_label)) +
  geom_violin(scale = "width", trim = FALSE) +
  scale_fill_manual(values = lineage_colors) +
  labs(
    x = "Normalized Pseudotime",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    panel.grid.major.y = element_blank()
  )

print(p_fig2e)
ggsave("figures/Figure2e_pseudotime_distribution.png", p_fig2e, 
       width = 8, height = 6, dpi = 600)

# ============================================================================
# Prepare for tradeSeq
# ============================================================================

message("\n--- Preparing for tradeSeq ---")

# Normalize pseudotime to 0-1 range
pseudotime_norm <- (pseudotime_final - min(pseudotime_final)) /
  (max(pseudotime_final) - min(pseudotime_final))

# Add to colData
colData(sce_final)$curve1 <- pseudotime_norm
colData(sce_final)$lineageWeight1 <- rep(1, ncol(sce_final))
colData(sce_final)$primary_lineage <- primary_lineage_final

# Verification
n_na_pst <- sum(is.na(colData(sce_final)$curve1))
n_na_weight <- sum(is.na(colData(sce_final)$lineageWeight1))

message("NAs in curve1: ", n_na_pst)
message("NAs in lineageWeight1: ", n_na_weight)

if (n_na_pst == 0 && n_na_weight == 0) {
  message("✓ No NAs detected - Ready for tradeSeq!")
} else {
  warning("⚠ NAs detected in pseudotime/weights!")
}

# ============================================================================
# Lineage composition analysis
# ============================================================================

message("\n--- Lineage composition analysis ---\n")

for (i in 1:n_lineages) {
  message(sprintf("Lineage %d:", i))
  lin_cells <- colData(sce_final)$primary_lineage == i
  
  if (sum(lin_cells) > 0) {
    ct_dist <- table(sce_final$celltype_merged[lin_cells])
    print(ct_dist)
    
    # Find dominant cell type
    dominant <- names(ct_dist)[which.max(ct_dist)]
    pct_dominant <- round(max(ct_dist) / sum(ct_dist) * 100, 1)
    message(sprintf("  → Dominant: %s (%.1f%%)\n", dominant, pct_dominant))
  }
}

# Final cell type distribution
message("Final cell type distribution:")
ct_table <- table(sce_final$celltype_merged)
for (ct in names(ct_table)) {
  pct <- round(ct_table[ct] / sum(ct_table) * 100, 1)
  message(sprintf("  %s: %d (%.1f%%)", ct, ct_table[ct], pct))
}


# ============================================================================
# Save final object
# ============================================================================

message("\n=== Saving final results ===")
saveRDS(sce_final, "data/sce_trajectory.rds")
message("✓ Saved: data/sce_trajectory.rds")
message("  - Contains: ", ncol(sce_final), " cells × ", nrow(sce_final), " genes")
message("  - Pseudotime: curve1 (normalized 0-1)")
message("  - Primary lineage assignments")
message("  - Cell types: ", paste(names(table(sce_final$celltype_merged)), collapse=", "))
