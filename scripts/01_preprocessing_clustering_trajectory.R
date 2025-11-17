# ============================================================================
# 01_preprocessing_clustering_trajectory.R
# Description: Complete pipeline from raw data to trajectory inference
#              (Preprocessing → ScType annotation → Slingshot trajectory)
# 
# NOTE: This script requires raw count matrices from NCBI SRA BioProject PRJNA1256192
# Raw data files (not included):
#   - 1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv
#   - 1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv
# 
# For reviewers: If starting from preprocessed data, load sce_trajectory.rds directly
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

SEED <- 123
set.seed(SEED)

# ============================================================================
# PART 1: Data Loading and Preprocessing
# ============================================================================

message("=== PART 1: Data Loading and Preprocessing ===\n")

# Load control sample
message("Loading shCtrl sample...")
counts <- read.table("data/raw/1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv",
                     skip = 5, sep = ",", header = TRUE, row.names = 1)
shCtrl <- CreateSeuratObject(counts = t(counts), project = "shCtrl")
message("  shCtrl: ", ncol(shCtrl), " cells")

# Load PTBP1 knockdown sample
message("Loading shPTBP1 sample...")
counts_sh <- read.table("data/raw/1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv",
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
source("scripts/utils/sctype_function_1.R")
source("scripts/utils/sctype_function_2.R")

# ScType database
db_file <- "data/ScTypeDB_full.xlsx"
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

# Visualize
p_sctype <- DimPlot(plus, reduction = "umap", label = TRUE, repel = TRUE, 
                    group.by = 'customclassif') +
  ggtitle("ScType Cell Type Annotation") +
  theme(legend.position = "right")
ggsave("figures/ScType_annotation.svg", p_sctype, width = 10, height = 8)
message("✓ Saved: figures/ScType_annotation.svg")

message("\nCell type distribution:")
print(table(plus$customclassif, plus$orig.ident))

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
message("\n=== CHECKPOINT 1: Saved seurat_preprocessed.rds ===")
message("Purpose: Gene set analysis (SynGO, pathway analysis)")
message("Content: ", nrow(plus_filtered), " genes × ", ncol(plus_filtered), " cells")
message("Cell types: ", paste(names(table(plus_filtered$customclassif)), collapse=", "))

# ============================================================================
# Convert to SingleCellExperiment
# ============================================================================

message("\n--- Preparing SingleCellExperiment object ---")

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

# Extract pseudotime and weights
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

message("\n--- Saving results ---")

saveRDS(sce_final, "data/sce_trajectory.rds")

message("\n=== CHECKPOINT 2: Saved sce_trajectory.rds ===")
message("Purpose: Trajectory inference (Slingshot, tradeSeq)")
message("Content: ", nrow(sce_final), " genes × ", ncol(sce_final), " cells")
message("Note: Subset to highly variable genes for trajectory analysis")
message("\nThis object contains:")
message("  - Gene expression counts")
message("  - Cell type annotations")
message("  - Trajectory pseudotime (curve1)")
message("  - Primary lineage assignments")
message("  - Ready for Countsplit DEG analysis")

message("\n=== Pipeline completed successfully! ===")
message("Next step: Run 02_DEG_countsplit.R")