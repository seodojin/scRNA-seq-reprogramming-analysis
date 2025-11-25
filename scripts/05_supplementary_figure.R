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

# Visualize feature plots
FeaturePlot(plus, features = c("GRIN2A", "GABBR2", "ACHE", "TH"))

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
# Filter Unknown cells
# ============================================================================

message("\n--- Filtering cell types ---")

# Remove Unknown cells
message("Before filtering - Total cells: ", ncol(plus))
message("Unknown cells: ", sum(plus@meta.data$customclassif == "Unknown"))

plus <- subset(plus, subset = customclassif != "Unknown")
message("After filtering Unknown - Remaining cells: ", ncol(plus))

# ============================================================================
# PART 3: Differential Expression Analysis (FindAllMarkers)
# ============================================================================

message("\n=== PART 3: Differential Expression Analysis ===\n")

# Set identity to cell type annotation
Idents(plus) <- "customclassif"

# Define cluster order
cluster_order <- c("Fibroblasts", "GABAergic neurons", "Glutamatergic neurons",
                   "Immature neurons", "Myofibroblasts")

# Join layers for differential expression analysis
plus <- JoinLayers(plus)

# Find all markers
message("Running FindAllMarkers...")
markers_all <- FindAllMarkers(
  plus,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  p.adjust.method = "BH"
) %>% filter(p_val_adj < 0.05)

message("  ✓ Found ", nrow(markers_all), " significant markers")

# ============================================================================
# Select Top Markers
# ============================================================================

# Select top 5 markers per cluster for heatmap
heatmap_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

# Extract gene names in cluster order
features <- heatmap_markers %>%
  arrange(match(cluster, cluster_order), desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

# ============================================================================
# Balanced Sampling for Heatmap (max 300 cells per cluster)
# ============================================================================

cells_bal <- unlist(lapply(cluster_order, function(cl){
  cs <- WhichCells(plus, idents = cl)
  if (length(cs) > 300) sample(cs, 300) else cs
}))

message("  ✓ Sampled ", length(cells_bal), " cells for heatmap")

# ============================================================================
# Generate Heatmap
# ============================================================================

message("\n--- Generating heatmap ---")

p_heatmap <- DoHeatmap(plus[, cells_bal],
                       features = features,
                       group.by = "customclassif",
                       raster = TRUE) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(size = 8))

# Save heatmap
ggsave("figures/TopMarkers_Heatmap.png",
       plot = p_heatmap,
       width = 12, height = 10,
       dpi = 1000)

message("  ✓ Saved: TopMarkers_Heatmap.png")

# ============================================================================
# Save Top 10 Markers per Cluster as CSV
# ============================================================================

message("\n--- Saving top 10 markers per cluster ---")

# Select top 10 markers per cluster
top10_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

# Sort by cluster order
top10_markers_sorted <- top10_markers %>%
  arrange(match(cluster, cluster_order), desc(avg_log2FC))

# Save as CSV
write.csv(top10_markers_sorted, 
          file = "Top10_Markers_per_Cluster.csv",
          row.names = FALSE)

message("  ✓ Saved ", nrow(top10_markers_sorted), " markers to Top10_Markers_per_Cluster.csv")

# ============================================================================
# PART 4: Canonical Marker Analysis
# ============================================================================

message("\n=== PART 4: Canonical Marker Analysis ===\n")

# Define canonical markers (6 genes)
canonical_markers <- c(
  "TUBA1C",    # Immature neurons
  "ACTA2",     # Myofibroblasts  
  "COL1A1",    # Fibroblasts
  "STMN2",     # Pan-neuronal (mature)
  "GRIN2A",    # Glutamatergic-specific
  "GABBR2"     # GABAergic-specific
)

# Check availability of markers
available_markers <- canonical_markers[canonical_markers %in% rownames(plus)]

message("Canonical markers:")
marker_info <- data.frame(
  Marker = canonical_markers,
  Target = c("Immature neurons", "Myofibroblasts", "Fibroblasts",
             "Pan-neuronal", "Glutamatergic", "GABAergic"),
  Available = canonical_markers %in% rownames(plus)
)
print(marker_info)

# ============================================================================
# Statistical Analysis of Canonical Markers
# ============================================================================

message("\n--- Statistical analysis of canonical markers ---")

canonical_stats <- data.frame()

for(marker in available_markers) {
  expr_data <- FetchData(plus, vars = c(marker, "customclassif"))
  colnames(expr_data) <- c("expression", "cluster")
  
  # Calculate statistics per cluster
  cluster_stats <- expr_data %>%
    group_by(cluster) %>%
    summarise(
      n_cells = n(),
      mean_expr = round(mean(expression), 4),
      median_expr = round(median(expression), 4),
      pct_expressing = round(sum(expression > 0) / n() * 100, 2),
      .groups = 'drop'
    )
  
  # Kruskal-Wallis test
  kruskal_result <- kruskal.test(expression ~ cluster, data = expr_data)
  
  # Store results
  for(i in 1:nrow(cluster_stats)) {
    canonical_stats <- rbind(canonical_stats, data.frame(
      Gene = marker,
      Cell_Type = cluster_stats$cluster[i],
      N_Cells = cluster_stats$n_cells[i],
      Mean_Expression = cluster_stats$mean_expr[i],
      Median_Expression = cluster_stats$median_expr[i],
      Pct_Expressing = cluster_stats$pct_expressing[i],
      Kruskal_Pvalue = ifelse(i == 1,
                              format(kruskal_result$p.value, scientific = TRUE, digits = 3),
                              ""),
      stringsAsFactors = FALSE
    ))
  }
}

# Save statistics
write.csv(canonical_stats,
          "Supp_Table_Canonical_Markers_Stats.csv",
          row.names = FALSE)

message("  ✓ Saved: Supp_Table_Canonical_Markers_Stats.csv")

# ============================================================================
# Generate Canonical Marker Violin Plots (Supp Fig 2A)
# ============================================================================

message("\n--- Generating canonical marker violin plots ---")

violin_plots <- list()

for(marker in available_markers) {
  p <- VlnPlot(plus,
               features = marker,
               group.by = "customclassif",
               pt.size = 0) +
    NoLegend() +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 10, color = "black"),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold", color = "black"),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(15, 15, 15, 15),
      legend.position = "none"
    ) +
    ggtitle(marker) +
    ylab("Expression Level")
  
  violin_plots[[marker]] <- p
}

# Combine plots (3x2 layout)
combined_violins <- wrap_plots(violin_plots, ncol = 3)

# Save violin plots
ggsave("figures/Supp_Fig2A_Canonical_Markers_Violins.png",
       combined_violins,
       width = 15, height = 10, dpi = 600)

message("  ✓ Saved: Supp_Fig2A_Canonical_Markers_Violins.png")
