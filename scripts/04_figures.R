# ============================================================================
# Publication Figure Generation Pipeline 
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(forcats)
})

message("=== Publication Figure Generation ===\n")

# Configuration
N_GENES_HEATMAP <- 50
N_TERMS_DOTPLOT <- 20
N_TOP_TFS <- 6
SEED <- 123
set.seed(SEED)

# ============================================================================
# PART 1: Load Data
# ============================================================================

message("\n=== Loading Data ===")

validated_degs <- read.csv("data/20251120_DEGs_Complete_Countsplit_Validated.csv")
message("✓ Validated genes: ", nrow(validated_degs))

all_results <- readRDS("data/20251120_Complete_Countsplit_All_Results.rds")
sce_gam_BA <- all_results$sce_gam_BA

sce_original <- readRDS("data/sce_trajectory.rds")
message("✓ Data loaded")

# ============================================================================
# PART 2: Gene Selection
# ============================================================================

message("\n=== Gene Selection ===")

pred_expr_all <- predictSmooth(
  models = sce_gam_BA,
  gene = validated_degs$gene,
  nPoints = 100,
  tidy = TRUE
)
message("✓ Generated predictions for all genes")

# Select top genes by effect size
top_genes <- validated_degs %>%
  arrange(desc(abs(mean_fcMedian))) %>%
  slice_head(n = N_GENES_HEATMAP) %>%
  pull(gene)

message("✓ Selected ", length(top_genes), " genes by effect size")

# ============================================================================
# PART 3: Heatmap
# ============================================================================

message("\n=== Generating Heatmap ===")

pred_expr_selected <- pred_expr_all %>%
  filter(gene %in% top_genes)

expr_wide <- pred_expr_selected %>%
  mutate(
    lineage_label = ifelse(lineage == 1, "Mature", "Immature"),
    time_point = paste0(lineage_label, "_", sprintf("%05.2f", time))
  ) %>%
  dplyr::select(gene, time_point, yhat) %>%
  pivot_wider(names_from = time_point, values_from = yhat) %>%
  column_to_rownames("gene") %>%
  as.matrix()

message("Expression matrix: ", nrow(expr_wide), " genes × ", ncol(expr_wide), " timepoints")

# Z-score normalization
expr_scaled <- t(scale(t(expr_wide)))

# Column annotation
col_df <- data.frame(
  Lineage = ifelse(grepl("^Mature_", colnames(expr_scaled)),
                   "Mature_Neurons", "Immature_Neurons"),
  row.names = colnames(expr_scaled)
)

ann_colors <- list(
  Lineage = c("Mature_Neurons" = "purple", "Immature_Neurons" = "yellow")
)

# Save heatmap
heatmap_file <- "Figure_Heatmap_Validated_DEGs.png"
message("Saving: ", heatmap_file)

png(heatmap_file, width = 13, height = 14, units = "in", res = 600)
pheatmap(
  expr_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = col_df,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "Validated DEGs (Mature vs Immature Neuronal Lineages)",
  border_color = NA
)
dev.off()

message("✓ Heatmap saved: ", heatmap_file)

# ============================================================================
# PART 4: Enrichment Analysis
# ============================================================================

message("\n=== Enrichment Analysis ===")

all_genes <- validated_degs %>%
  filter(consistent_direction == TRUE) %>%
  pull(gene)

message("Using ", length(all_genes), " genes with consistent direction")

entrez_conv <- tryCatch({
  bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) {
  message("Gene conversion error: ", e$message)
  return(NULL)
})

ego_bp <- NULL
ekegg <- NULL

if (!is.null(entrez_conv) && nrow(entrez_conv) > 0) {
  message("✓ Converted ", nrow(entrez_conv), "/", length(all_genes), " genes")
  
  # GO Biological Process
  ego_bp <- tryCatch({
    enrichGO(
      gene = entrez_conv$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }, error = function(e) {
    message("GO enrichment error: ", e$message)
    return(NULL)
  })
  
  if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
    message("✓ Found ", nrow(ego_bp), " GO BP terms")
    write.csv(as.data.frame(ego_bp), "Enrichment_GO_BP.csv", row.names = FALSE)
    message("  Saved: Enrichment_GO_BP.csv")
  } else {
    message("⚠ No significant GO terms found")
  }
  
  # KEGG Pathway
  ekegg <- tryCatch({
    enrichKEGG(
      gene = entrez_conv$ENTREZID,
      organism = "hsa",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    message("KEGG enrichment error: ", e$message)
    return(NULL)
  })
  
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    message("✓ Found ", nrow(ekegg), " KEGG pathways")
    write.csv(as.data.frame(ekegg), "Enrichment_KEGG.csv", row.names = FALSE)
    message("  Saved: Enrichment_KEGG.csv")
  } else {
    message("⚠ No significant KEGG pathways found")
  }
}

# ============================================================================
# PART 5: Dotplots
# ============================================================================

message("\n=== Creating Individual Dotplots ===")

make_dotplot <- function(enrich_result, title = "Enrichment", n_terms = 20) {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    message("  No results for ", title)
    return(NULL)
  }
  
  plot_data <- as.data.frame(enrich_result) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_terms) %>%
    mutate(
      Description = fct_reorder(Description, -p.adjust),
      GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"), 
                                 function(x) as.numeric(x[1]) / as.numeric(x[2]))
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

# Individual GO plot
if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
  p_go <- make_dotplot(ego_bp, "GO Biological Process Enrichment", N_TERMS_DOTPLOT)
  if (!is.null(p_go)) {
    ggsave("Figure_Dotplot_GO_BP.png", p_go, width = 10, height = 8, dpi = 600)
    message("✓ Saved: Figure_Dotplot_GO_BP.png")
  }
}

# Individual KEGG plot
if (!is.null(ekegg) && nrow(ekegg) > 0) {
  p_kegg <- make_dotplot(ekegg, "KEGG Pathway Enrichment", N_TERMS_DOTPLOT)
  if (!is.null(p_kegg)) {
    ggsave("Figure_Dotplot_KEGG.png", p_kegg, width = 10, height = 8, dpi = 600)
    message("✓ Saved: Figure_Dotplot_KEGG.png")
  }
}

# ============================================================================
# PART 6: TF Expression Plots
# ============================================================================

message("\n=== Individual TF Expression Plots ===")

tf_file <- "20251120_Validated_TFs_Complete_Annotation.csv"

if (file.exists(tf_file)) {
  validated_tfs <- read.csv(tf_file)
  message("✓ Loaded ", nrow(validated_tfs), " validated TFs")
  
  top_tfs <- validated_tfs %>%
    filter(consistent_direction == TRUE) %>%
    arrange(max_padj) %>%
    slice_head(n = N_TOP_TFS) %>%
    pull(gene)
  
  if (length(top_tfs) > 0) {
    message("Top ", length(top_tfs), " TFs: ", paste(top_tfs, collapse = ", "))
    
    plot_tf_smoother <- function(models, gene) {
      pred_data <- predictSmooth(
        models = models,
        gene = gene,
        nPoints = 100,
        tidy = TRUE
      )
      
      pred_data <- pred_data %>%
        group_by(lineage) %>%
        mutate(time_norm = (time - min(time)) / (max(time) - min(time))) %>%
        ungroup() %>%
        mutate(
          expr = log1p(yhat),
          Lineage = factor(ifelse(lineage == 1, "Mature", "Immature"))
        )
      
      ggplot(pred_data, aes(x = time_norm, y = expr, color = Lineage)) +
        geom_line(linewidth = 1.5, alpha = 0.8) +
        scale_color_manual(values = c("Mature" = "#7570B3", "Immature" = "#E7B800")) +
        scale_x_continuous(limits = c(0, 1)) +
        labs(x = "Normalized Pseudotime", y = "Log(Expression + 1)", title = gene) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
        )
    }
    
    # Individual TF plots only
    message("\nGenerating individual TF plots...")
    n_success <- 0
    for (tf in top_tfs) {
      if (tf %in% rownames(sce_gam_BA)) {
        tryCatch({
          p <- plot_tf_smoother(sce_gam_BA, tf)
          tf_filename <- paste0("TF_", tf, ".png")
          ggsave(tf_filename, p, width = 5, height = 4, dpi = 600)
          message("  ✓ ", tf_filename)
          n_success <- n_success + 1
        }, error = function(e) {
          message("  ✗ Error plotting ", tf, ": ", e$message)
        })
      } else {
        message("  ✗ ", tf, " not in GAM model")
      }
    }
    
    if (n_success == 0) {
      message("⚠ No TFs could be plotted")
    }
  } else {
    message("⚠ No TFs with consistent direction found")
  }
} else {
  message("⚠ TF annotation file not found: ", tf_file)
  message("  Run TF analysis first")
}

message("\n=== Complete ===")