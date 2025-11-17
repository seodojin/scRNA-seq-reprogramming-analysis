# ============================================================================
# 04_figures.R
# Description: Generate all publication-quality figures
#              - Heatmap: Predicted gene expression patterns
#              - Enrichment dotplots: GO BP and KEGG
#              - TF smoother plots: Individual and combined
# ============================================================================

suppressPackageStartupMessages({
  library(tradeSeq)
  library(slingshot)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(forcats)
  library(cowplot)
  library(patchwork)
})

message("=== Publication Figure Generation Pipeline ===\n")

# ============================================================================
# Configuration
# ============================================================================

USE_FULL_COUNTS <- TRUE  # Use full counts for cleaner visualization
N_GENES_HEATMAP <- 50
N_TERMS_DOTPLOT <- 20
N_TOP_TFS <- 6           # Number of top TFs to plot
SEED <- 123

message("Configuration:")
message("  Use full counts: ", USE_FULL_COUNTS)
message("  Genes in heatmap: ", N_GENES_HEATMAP)
message("  Terms in dotplot: ", N_TERMS_DOTPLOT)
message("  Top TFs to plot: ", N_TOP_TFS)

# ============================================================================
# PART 1: Load Data
# ============================================================================

message("\n=== PART 1: Loading Data ===\n")

# Load validated genes and TFs
sig_pat <- read.csv("results/patternTest_splitB_testing_validated.csv")
validated_tfs <- read.csv("results/validated_TFs_annotated.csv")

message("Validated genes loaded: ", nrow(sig_pat))
message("Validated TFs loaded: ", nrow(validated_tfs))

# Load trajectory objects
sce_final <- readRDS("data/sce_trajectory.rds")
message("SingleCellExperiment loaded: ", ncol(sce_final), " cells")

# Load GAM model or prepare to refit
if (USE_FULL_COUNTS) {
  message("Will refit GAM on full counts for visualization")
} else {
  sce_gam_both_B <- readRDS("data/sce_gam_lineages.rds")
  message("Loaded Split B GAM model")
}

# ============================================================================
# PART 2: Prepare Lineage Data
# ============================================================================

message("\n=== PART 2: Preparing Lineage Data ===\n")

# Extract lineage information
primary_lineage <- colData(sce_final)$primary_lineage
pseudotime_col <- "curve1"

# Identify mature and immature lineages
lineage_composition <- table(
  sce_final$celltype_merged,
  primary_lineage
)

lineage_props <- prop.table(lineage_composition, margin = 2)
mature_lin <- which.max(lineage_props["Neurons", ])
immature_lin <- which.max(lineage_props["Immature neurons", ])

message("Mature lineage: ", mature_lin)
message("Immature lineage: ", immature_lin)

# Select cells in both lineages
lin1_cells <- primary_lineage == mature_lin
lin3_cells <- primary_lineage == immature_lin
both_cells <- lin1_cells | lin3_cells

message("Cells: Mature=", sum(lin1_cells), ", Immature=", sum(lin3_cells))

# Prepare pseudotime and weights matrices
n_both <- sum(both_cells)
pseudotime_both <- matrix(0, nrow = n_both, ncol = 2)
weights_both <- matrix(0, nrow = n_both, ncol = 2)

colnames(pseudotime_both) <- c("lineage_mature", "lineage_immature")
colnames(weights_both) <- c("lineage_mature", "lineage_immature")

all_pseudotime <- colData(sce_final)[[pseudotime_col]]

lin1_in_both <- which(lin1_cells[both_cells])
pseudotime_both[lin1_in_both, 1] <- all_pseudotime[both_cells][lin1_in_both]
weights_both[lin1_in_both, 1] <- 1

lin3_in_both <- which(lin3_cells[both_cells])
pseudotime_both[lin3_in_both, 2] <- all_pseudotime[both_cells][lin3_in_both]
weights_both[lin3_in_both, 2] <- 1

# ============================================================================
# PART 3: Fit GAM for Visualization
# ============================================================================

message("\n=== PART 3: Fitting GAM Model ===\n")

if (USE_FULL_COUNTS) {
  message("Fitting GAM on FULL counts...")
  
  counts_full <- assay(sce_final, "counts")
  counts_both <- counts_full[, both_cells]
  
  min_cells <- ceiling(ncol(counts_both) * 0.05)
  keep_genes <- rowSums(counts_both > 0) >= min_cells
  counts_filtered <- counts_both[keep_genes, ]
  
  message("Genes after filtering: ", nrow(counts_filtered))
  message("Fitting GAM (1-2 minutes)...")
  
  set.seed(SEED)
  sce_gam_vis <- fitGAM(
    counts = counts_filtered,
    pseudotime = pseudotime_both,
    cellWeights = weights_both,
    nknots = 4,
    verbose = FALSE
  )
  
  message("✓ GAM fitting completed")
  
} else {
  sce_gam_vis <- sce_gam_both_B
}

# ============================================================================
# PART 4: Heatmap
# ============================================================================

message("\n=== PART 4: Generating Heatmap ===\n")

# Select genes with diverse representation
sig_pat_diverse <- sig_pat %>%
  mutate(first_letter = substr(gene, 1, 1)) %>%
  group_by(first_letter) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  arrange(padj) %>%
  slice_head(n = N_GENES_HEATMAP)

top_genes <- sig_pat_diverse$gene
message("Selected ", length(top_genes), " genes for heatmap")

# Generate predictions
pred_expr_tidy <- predictSmooth(
  models = sce_gam_vis,
  gene = top_genes,
  nPoints = 100,
  tidy = TRUE
)

# Convert to wide format
expr_wide <- pred_expr_tidy %>%
  mutate(
    lineage_label = ifelse(lineage == 1, "Mature", "Immature"),
    time_rounded = round(time, 2),
    time_point = paste0(lineage_label, "_", time_rounded)
  ) %>%
  dplyr::select(gene, time_point, yhat) %>%
  group_by(gene, time_point) %>%
  summarise(yhat = mean(yhat, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = time_point, values_from = yhat) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Z-score normalization
expr_scaled <- t(scale(t(expr_wide)))

# Column annotation
col_df <- data.frame(
  Lineage = ifelse(grepl("^Mature_", colnames(expr_scaled)),
                   "Mature_Neurons", "Immature_Neurons"),
  row.names = colnames(expr_scaled)
)

ann_colors <- list(
  Lineage = c(
    "Mature_Neurons" = "purple",
    "Immature_Neurons" = "yellow"
  )
)

# Save heatmap
svg("figures/Figure_heatmap_patternTest.svg", width = 12, height = 14)
pheatmap(
  expr_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = col_df,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "Differentially Expressed Genes\n(Mature vs Immature Neuronal Lineages)",
  border_color = NA
)
dev.off()

message("✓ Saved: figures/Figure_heatmap_patternTest.svg\n")

# ============================================================================
# PART 5: Enrichment Analysis
# ============================================================================

message("=== PART 5: GO Enrichment Analysis ===\n")

all_genes <- sig_pat$gene

# Convert to Entrez IDs
entrez_conversion <- clusterProfiler::bitr(
  all_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

message("Converted: ", nrow(entrez_conversion), "/", length(all_genes), " genes")

# GO BP enrichment
ego_bp <- clusterProfiler::enrichGO(
  gene = entrez_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
  message("✓ Found ", nrow(ego_bp), " GO BP terms")
  write.csv(as.data.frame(ego_bp),
            "results/enrichment_GO_BP.csv", row.names = FALSE)
}

# KEGG enrichment
ekegg <- clusterProfiler::enrichKEGG(
  gene = entrez_conversion$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

if (!is.null(ekegg) && nrow(ekegg) > 0) {
  ekegg <- clusterProfiler::setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  message("✓ Found ", nrow(ekegg), " KEGG pathways")
  write.csv(as.data.frame(ekegg),
            "results/enrichment_KEGG.csv", row.names = FALSE)
}

# ============================================================================
# PART 6: Enrichment Dotplots
# ============================================================================

message("\n=== PART 6: Creating Enrichment Dotplots ===\n")

make_dotplot <- function(enrich_result, title = "Enrichment", n_terms = 20) {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) return(NULL)
  
  plot_data <- as.data.frame(enrich_result) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_terms) %>%
    mutate(
      Description = fct_reorder(Description, -p.adjust),
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        nums <- as.numeric(strsplit(x, "/")[[1]])
        nums[1] / nums[2]
      })
    )
  
  ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradientn(
      colors = c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"),
      name = "Adjusted\nP-value",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
}

# GO dotplot
if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
  p_go <- make_dotplot(ego_bp, "GO Biological Process Enrichment", N_TERMS_DOTPLOT)
  ggsave("figures/Figure_dotplot_GO_BP.svg", p_go, width = 10, height = 8, device = "svg")
  message("✓ Saved: figures/Figure_dotplot_GO_BP.svg")
}

# KEGG dotplot
if (!is.null(ekegg) && nrow(ekegg) > 0) {
  p_kegg <- make_dotplot(ekegg, "KEGG Pathway Enrichment", N_TERMS_DOTPLOT)
  ggsave("figures/Figure_dotplot_KEGG.svg", p_kegg, width = 10, height = 8, device = "svg")
  message("✓ Saved: figures/Figure_dotplot_KEGG.svg")
}

# Combined
if (!is.null(ego_bp) && nrow(ego_bp) > 0 && !is.null(ekegg) && nrow(ekegg) > 0) {
  p_combined <- plot_grid(
    make_dotplot(ego_bp, "GO Biological Process", 15),
    make_dotplot(ekegg, "KEGG Pathway", 15),
    labels = c("B", "C"), label_size = 16, ncol = 1, align = "v"
  )
  ggsave("figures/Figure_combined_GO_KEGG.svg", p_combined,
         width = 10, height = 14, device = "svg")
  message("✓ Saved: figures/Figure_combined_GO_KEGG.svg")
}

# ============================================================================
# PART 7: Transcription Factor Smoother Plots
# ============================================================================

message("\n=== PART 7: TF Smoother Plots ===\n")

# Custom function for ggplot smoothers
plot_smoothers_ggplot <- function(models, gene, log_scale = TRUE, normalize_time = TRUE) {
  
  pred_data <- predictSmooth(models = models, gene = gene, nPoints = 100, tidy = TRUE)
  
  if (normalize_time) {
    pred_data <- pred_data %>%
      group_by(lineage) %>%
      mutate(time_normalized = (time - min(time)) / (max(time) - min(time))) %>%
      ungroup()
    time_var <- "time_normalized"
    x_label <- "Normalized Pseudotime"
  } else {
    time_var <- "time"
    x_label <- "Pseudotime"
  }
  
  if (log_scale) {
    pred_data <- pred_data %>% mutate(yhat_transformed = log1p(yhat))
    y_var <- "yhat_transformed"
    y_label <- "Log(Expression + 1)"
  } else {
    y_var <- "yhat"
    y_label <- "Expression"
  }
  
  pred_data <- pred_data %>%
    mutate(Lineage = factor(paste0("Lineage ", lineage)))
  
  ggplot(pred_data, aes(x = .data[[time_var]], y = .data[[y_var]], color = Lineage)) +
    geom_line(size = 1.5, alpha = 0.75) +
    scale_color_manual(values = c("Lineage 1" = "purple", "Lineage 2" = "yellow")) +
    scale_x_continuous(limits = if(normalize_time) c(0, 1) else NULL) +
    labs(x = x_label, y = y_label, title = gene) +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
}

# Select top TFs
if (nrow(validated_tfs) > 0) {
  top_tfs <- head(validated_tfs$gene, N_TOP_TFS)
  message("Plotting ", length(top_tfs), " top TFs:")
  print(top_tfs)
  
  # Individual plots
  message("\nGenerating individual TF plots...")
  for (tf in top_tfs) {
    tryCatch({
      p <- plot_smoothers_ggplot(sce_gam_vis, tf, log_scale = TRUE, normalize_time = TRUE) +
        ggtitle(tf) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14))
      
      filename <- paste0("figures/TF_", tf, ".svg")
      ggsave(filename, p, width = 4, height = 3, device = "svg")
      message("  ✓ ", filename)
    }, error = function(e) {
      message("  ✗ Error plotting ", tf)
    })
  }
  
  # Combined plot
  message("\nGenerating combined TF plot...")
  tf_plots_list <- lapply(top_tfs, function(tf) {
    plot_smoothers_ggplot(sce_gam_vis, tf, log_scale = TRUE, normalize_time = TRUE) +
      ggtitle(tf) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 0.6)
      )
  })
  
  combined <- wrap_plots(tf_plots_list, ncol = 3, nrow = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Transcription Factor Expression Patterns",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  
  ggsave("figures/TF_combined.svg", combined, width = 12, height = 8, device = "svg")
  message("✓ Saved: figures/TF_combined.svg")
}

# ============================================================================
# Summary
# ============================================================================

message("\n=== SUMMARY ===")
message("Generated figures:")
message("  1. Heatmap: figures/Figure_heatmap_patternTest.svg")
message("  2. GO dotplot: figures/Figure_dotplot_GO_BP.svg")
message("  3. KEGG dotplot: figures/Figure_dotplot_KEGG.svg")
message("  4. Combined: figures/Figure_combined_GO_KEGG.svg")
if (nrow(validated_tfs) > 0) {
  message("  5. Individual TF plots: figures/TF_*.svg (", length(top_tfs), " files)")
  message("  6. Combined TF plot: figures/TF_combined.svg")
}

message("\n✓ All publication figures generated successfully!")