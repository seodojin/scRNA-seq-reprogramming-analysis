library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(Matrix)
library(dplyr)

SEED <- 123
set.seed(SEED)

message("=== Complete Countsplit Validation Pipeline ===\n")

# ============================================================
# Load Data
# ============================================================

sce_original <- readRDS("data/sce_trajectory.rds")

message("Loaded: ", nrow(sce_original), " genes × ", 
        ncol(sce_original), " cells")
message("Cell types:")
print(table(sce_original$celltype_merged))

# Extract components
counts_original <- assay(sce_original, "counts")
umap_coords <- reducedDim(sce_original, "UMAP")
pca_coords <- reducedDim(sce_original, "PCA")
metadata <- colData(sce_original)

# ============================================================
# Split Counts
# ============================================================

split_binomial_sparse <- function(mat, p = 0.5, seed = 123) {
  set.seed(seed)
  mt <- as(mat, "TsparseMatrix")
  if (length(mt@x) == 0L) return(list(A = mat, B = mat))
  xA <- rbinom(length(mt@x), size = mt@x, prob = p)
  xB <- mt@x - xA
  A <- sparseMatrix(i = mt@i + 1, j = mt@j + 1, x = xA,
                    dims = dim(mat), dimnames = dimnames(mat))
  B <- sparseMatrix(i = mt@i + 1, j = mt@j + 1, x = xB,
                    dims = dim(mat), dimnames = dimnames(mat))
  list(A = as(A, "dgCMatrix"), B = as(B, "dgCMatrix"))
}

message("\n=== Splitting Counts ===")
splits <- split_binomial_sparse(counts_original, p = 0.5, seed = SEED)

message("✓ Split complete:")
message("  Original: ", sum(counts_original), " UMIs")
message("  Split A:  ", sum(splits$A), " UMIs (", 
        round(sum(splits$A)/sum(counts_original)*100, 1), "%)")
message("  Split B:  ", sum(splits$B), " UMIs (", 
        round(sum(splits$B)/sum(counts_original)*100, 1), "%)")

# ============================================================
# Helper Function
# ============================================================

create_sce_from_counts <- function(counts, metadata, 
                                   umap_coords, pca_coords,
                                   min_pct = 0.05) {
  # Gene filtering (keep genes expressed in at least 5% of cells)
  min_cells <- ceiling(ncol(counts) * min_pct)
  keep_genes <- rowSums(counts > 0) >= min_cells
  counts_filtered <- counts[keep_genes, ]
  
  message("    Genes: ", nrow(counts_filtered), "/", nrow(counts), 
          " retained (", round(mean(keep_genes)*100, 1), "%)")
  
  # Create SCE
  sce <- SingleCellExperiment(
    assays = list(counts = counts_filtered),
    colData = metadata
  )
  
  # Add dimensionality reductions
  reducedDim(sce, "UMAP") <- umap_coords
  reducedDim(sce, "PCA") <- pca_coords
  
  return(sce)
}

# ============================================================
# PATH 1: Trajectory A → DE B
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PATH 1: Trajectory inference with Split A")
message("        Differential expression with Split B")
message(paste(rep("=", 60), collapse = ""))

## Step 1A: Create SCE with Split A
message("\n[1A] Creating SCE with Split A counts...")
sce_A <- create_sce_from_counts(
  splits$A, 
  metadata,
  umap_coords,
  pca_coords,
  min_pct = 0.05
)

## Step 1B: Run Slingshot on Split A
message("\n[1B] Running Slingshot on Split A...")
sds_A <- slingshot(
  sce_A,
  clusterLabels = 'celltype_merged',
  reducedDim = 'UMAP',
  start.clus = "Fibroblasts",
  end.clus = c("Neurons", "Immature neurons", "Myofibroblasts"),
  allow.breaks = FALSE,
  extend = 'n',
  approx_points = 300
)

n_lineages <- ncol(slingPseudotime(sds_A))
message("      ✓ ", n_lineages, " lineages inferred")

## Step 1C: Extract trajectories
pseudotime_A <- slingPseudotime(sds_A)
weights_A <- slingCurveWeights(sds_A)

# Identify lineages (1=Neurons, 3=Immature)
mature_lin <- 1
immature_lin <- 3

lin1_cells <- !is.na(pseudotime_A[, mature_lin])
lin3_cells <- !is.na(pseudotime_A[, immature_lin])
both_cells <- lin1_cells | lin3_cells

message("      Lineage 1 (Neurons): ", sum(lin1_cells), " cells")
message("      Lineage 3 (Immature): ", sum(lin3_cells), " cells")
message("      Total in both: ", sum(both_cells), " cells")

## Step 1D: Prepare pseudotime/weights matrices
n_both <- sum(both_cells)
pseudotime_both_A <- matrix(0, nrow = n_both, ncol = 2)
weights_both_A <- matrix(0, nrow = n_both, ncol = 2)
colnames(pseudotime_both_A) <- c("lineage1", "lineage3")
colnames(weights_both_A) <- c("lineage1", "lineage3")

both_cells_idx <- which(both_cells)
lin1_in_both <- which(lin1_cells[both_cells])
lin3_in_both <- which(lin3_cells[both_cells])

pseudotime_both_A[lin1_in_both, 1] <- pseudotime_A[both_cells_idx[lin1_in_both], mature_lin]
weights_both_A[lin1_in_both, 1] <- weights_A[both_cells_idx[lin1_in_both], mature_lin]
pseudotime_both_A[lin3_in_both, 2] <- pseudotime_A[both_cells_idx[lin3_in_both], immature_lin]
weights_both_A[lin3_in_both, 2] <- weights_A[both_cells_idx[lin3_in_both], immature_lin]

## Step 1E: DE analysis with Split B counts
message("\n[1C] Fitting GAM with Split B counts...")
counts_B_both <- splits$B[, both_cells]
keep_genes_B <- rowSums(counts_B_both > 0) >= ceiling(n_both * 0.05)
counts_B_filtered <- counts_B_both[keep_genes_B, ]

message("      Cells: ", ncol(counts_B_filtered))
message("      Genes: ", nrow(counts_B_filtered))

set.seed(SEED)
sce_gam_BA <- fitGAM(
  counts = counts_B_filtered,
  pseudotime = pseudotime_both_A,  # From Split A
  cellWeights = weights_both_A,    # From Split A
  nknots = 4,
  verbose = FALSE
)

message("      ✓ GAM fitted")

## Step 1F: Pattern test
message("\n[1D] Running patternTest...")
pattern_BA <- patternTest(sce_gam_BA, l2fc = log2(1.5)) %>%
  as.data.frame() %>%
  mutate(
    gene = rownames(.),
    padj_BA = p.adjust(pvalue, "BH")
  ) %>%
  arrange(pvalue)

n_sig_BA <- sum(pattern_BA$padj_BA < 0.05)
message("      ✓ Pattern test complete")
message("      Significant genes (FDR < 0.05): ", n_sig_BA)

# ============================================================
# PATH 2: Trajectory B → DE A
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PATH 2: Trajectory inference with Split B")
message("        Differential expression with Split A")
message(paste(rep("=", 60), collapse = ""))

## Step 2A: Create SCE with Split B
message("\n[2A] Creating SCE with Split B counts...")
sce_B <- create_sce_from_counts(
  splits$B,
  metadata,
  umap_coords,
  pca_coords,
  min_pct = 0.05
)

## Step 2B: Run Slingshot on Split B
message("\n[2B] Running Slingshot on Split B...")
sds_B <- slingshot(
  sce_B,
  clusterLabels = 'celltype_merged',
  reducedDim = 'UMAP',
  start.clus = "Fibroblasts",
  end.clus = c("Neurons", "Immature neurons", "Myofibroblasts"),
  allow.breaks = FALSE,
  extend = 'n',
  approx_points = 300
)

n_lineages_B <- ncol(slingPseudotime(sds_B))
message("      ✓ ", n_lineages_B, " lineages inferred")

## Step 2C: Extract trajectories
pseudotime_B <- slingPseudotime(sds_B)
weights_B <- slingCurveWeights(sds_B)

## Step 2D: Prepare matrices (using same cell selection)
pseudotime_both_B <- matrix(0, nrow = n_both, ncol = 2)
weights_both_B <- matrix(0, nrow = n_both, ncol = 2)
colnames(pseudotime_both_B) <- c("lineage1", "lineage3")
colnames(weights_both_B) <- c("lineage1", "lineage3")

pseudotime_both_B[lin1_in_both, 1] <- pseudotime_B[both_cells_idx[lin1_in_both], mature_lin]
weights_both_B[lin1_in_both, 1] <- weights_B[both_cells_idx[lin1_in_both], mature_lin]
pseudotime_both_B[lin3_in_both, 2] <- pseudotime_B[both_cells_idx[lin3_in_both], immature_lin]
weights_both_B[lin3_in_both, 2] <- weights_B[both_cells_idx[lin3_in_both], immature_lin]

## Step 2E: DE analysis with Split A counts
message("\n[2C] Fitting GAM with Split A counts...")
counts_A_both <- splits$A[, both_cells]
keep_genes_A <- rowSums(counts_A_both > 0) >= ceiling(n_both * 0.05)
counts_A_filtered <- counts_A_both[keep_genes_A, ]

message("      Cells: ", ncol(counts_A_filtered))
message("      Genes: ", nrow(counts_A_filtered))

set.seed(SEED)
sce_gam_AB <- fitGAM(
  counts = counts_A_filtered,
  pseudotime = pseudotime_both_B,  # From Split B
  cellWeights = weights_both_B,    # From Split B
  nknots = 4,
  verbose = FALSE
)

message("      ✓ GAM fitted")

## Step 2F: Pattern test
message("\n[2D] Running patternTest...")
pattern_AB <- patternTest(sce_gam_AB, l2fc = log2(1.5)) %>%
  as.data.frame() %>%
  mutate(
    gene = rownames(.),
    padj_AB = p.adjust(pvalue, "BH")
  ) %>%
  arrange(pvalue)

n_sig_AB <- sum(pattern_AB$padj_AB < 0.05)
message("      ✓ Pattern test complete")
message("      Significant genes (FDR < 0.05): ", n_sig_AB)

# ============================================================
# VALIDATION: Intersection
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("VALIDATION: Computing Intersection")
message(paste(rep("=", 60), collapse = ""))

genes_BA_sig <- pattern_BA %>% filter(padj_BA < 0.05) %>% pull(gene)
genes_AB_sig <- pattern_AB %>% filter(padj_AB < 0.05) %>% pull(gene)

validated_genes <- intersect(genes_BA_sig, genes_AB_sig)
union_genes <- union(genes_BA_sig, genes_AB_sig)

message("\nResults:")
message("  Path A→B (Traj A, DE B): ", length(genes_BA_sig), " genes")
message("  Path B→A (Traj B, DE A): ", length(genes_AB_sig), " genes")
message("  Union: ", length(union_genes), " genes")
message("  INTERSECTION: ", length(validated_genes), " genes")
message("  Validation rate: ", 
        round(length(validated_genes) / length(union_genes) * 100, 1), "%")

# ============================================================
# Combine Results
# ============================================================

message("\n=== Creating Final Validated Gene List ===")

validated_results <- inner_join(
  pattern_BA %>% select(gene, 
                        waldStat_BA = waldStat, 
                        pvalue_BA = pvalue, 
                        padj_BA,
                        fcMedian_BA = fcMedian),
  pattern_AB %>% select(gene, 
                        waldStat_AB = waldStat, 
                        pvalue_AB = pvalue, 
                        padj_AB,
                        fcMedian_AB = fcMedian),
  by = "gene"
) %>%
  filter(gene %in% validated_genes) %>%
  mutate(
    max_padj = pmax(padj_BA, padj_AB),
    mean_waldStat = (waldStat_BA + waldStat_AB) / 2,
    mean_fcMedian = (fcMedian_BA + fcMedian_AB) / 2,
    both_positive_fc = (fcMedian_BA > 0 & fcMedian_AB > 0),
    both_negative_fc = (fcMedian_BA < 0 & fcMedian_AB < 0),
    consistent_direction = both_positive_fc | both_negative_fc
  ) %>%
  arrange(max_padj)

message("Final validated genes: ", nrow(validated_results))
message("  Consistent direction: ", sum(validated_results$consistent_direction))

# ============================================================
# Save Results
# ============================================================

message("\n=== Saving Results ===")

# Main results
write.csv(validated_results,
          "20251120_DEGs_Complete_Countsplit_Validated.csv",
          row.names = FALSE)

# Individual paths (for supplementary)
write.csv(pattern_BA,
          "20251120_DEGs_Path_A_to_B_Full.csv",
          row.names = FALSE)

write.csv(pattern_AB,
          "20251120_DEGs_Path_B_to_A_Full.csv",
          row.names = FALSE)

# Complete R object
saveRDS(list(
  parameters = list(
    seed = SEED,
    n_knots = 4,
    l2fc_threshold = log2(1.5),
    fdr_threshold = 0.05
  ),
  sds_A = sds_A,
  sds_B = sds_B,
  sce_gam_BA = sce_gam_BA,
  sce_gam_AB = sce_gam_AB,
  pattern_BA = pattern_BA,
  pattern_AB = pattern_AB,
  validated = validated_results
), "20251120_Complete_Countsplit_All_Results.rds")

message("\n✓ Files saved:")
message("  1. 20251120_DEGs_Complete_Countsplit_Validated.csv")
message("  2. 20251120_DEGs_Path_A_to_B_Full.csv")
message("  3. 20251120_DEGs_Path_B_to_A_Full.csv")
message("  4. 20251120_Complete_Countsplit_All_Results.rds")








# ============================================================
# Transcription Factor Analysis from Validated DEGs
# ============================================================

library(decoupleR)
library(OmnipathR)
library(dplyr)
library(tidyr)

message("\n", paste(rep("=", 70), collapse = ""))
message("TRANSCRIPTION FACTOR ANALYSIS")
message(paste(rep("=", 70), collapse = ""))

# ============================================================
# Step 1: Load Validated DEGs
# ============================================================

message("\n=== Step 1: Load Validated DEGs ===")

validated_degs <- read.csv("20251120_DEGs_Complete_Countsplit_Validated.csv")

message("Total validated genes: ", nrow(validated_degs))
message("  With consistent direction: ", sum(validated_degs$consistent_direction))

# ============================================================
# Step 2: Load DoRothEA Database
# ============================================================

message("\n=== Step 2: Load DoRothEA TF Database ===")

net <- get_dorothea(organism = 'human', levels = c('A', 'B', 'C'))

message("DoRothEA database loaded:")
message("  Total TF-target interactions: ", nrow(net))
message("  Unique TFs: ", length(unique(net$source)))

# Confidence level distribution
conf_dist <- table(net$confidence)
message("\nConfidence level distribution:")
for (level in names(conf_dist)) {
  message("  Level ", level, ": ", conf_dist[level], " interactions")
}

# ============================================================
# Step 3: Identify Validated TFs
# ============================================================

message("\n=== Step 3: Identify Validated TFs ===")

all_tfs <- unique(net$source)
message("Total TFs in DoRothEA: ", length(all_tfs))

# Find TFs in validated gene list
validated_tfs <- validated_degs %>%
  filter(gene %in% all_tfs) %>%
  arrange(max_padj)

message("\nValidated transcription factors:")
message("  Validated genes (FDR < 0.05): ", nrow(validated_degs))
message("  TFs among validated genes: ", nrow(validated_tfs))
message("  Proportion: ", round(nrow(validated_tfs)/nrow(validated_degs)*100, 1), "%")

if (nrow(validated_tfs) == 0) {
  message("\n⚠ No TFs found in validated gene list!")
  message("This may indicate:")
  message("  1. TFs are not differentially expressed")
  message("  2. Gene name mismatch")
  message("  3. Very stringent filtering")
  
  # Check for near-significant TFs
  message("\nChecking for TFs with relaxed threshold (max_padj < 0.1):")
  relaxed_tfs <- validated_degs %>%
    filter(gene %in% all_tfs, max_padj < 0.1) %>%
    arrange(max_padj)
  
  if (nrow(relaxed_tfs) > 0) {
    message("  Found ", nrow(relaxed_tfs), " TFs at FDR < 0.1")
    print(relaxed_tfs[, c("gene", "max_padj", "mean_fcMedian", "consistent_direction")])
  }
  
} else {
  
  # ============================================================
  # Step 4: Annotate TF Confidence and Targets
  # ============================================================
  
  message("\n=== Step 4: Annotate TF Confidence and Targets ===")
  
  ## 4A: Confidence levels
  tf_confidence <- net %>%
    filter(source %in% validated_tfs$gene) %>%
    group_by(source, confidence) %>%
    summarise(n = n(), .groups = 'drop') %>%
    group_by(source) %>%
    summarise(
      confidence_levels = paste(unique(confidence), collapse = ","),
      n_targets_total = sum(n),
      highest_confidence = min(confidence),
      .groups = 'drop'
    ) %>%
    rename(gene = source)
  
  message("✓ Confidence levels annotated")
  
  ## 4B: Mode of regulation
  tf_mor_raw <- net %>%
    filter(source %in% validated_tfs$gene) %>%
    group_by(source, mor) %>%
    summarise(n = n(), .groups = 'drop')
  
  message("  Mode of regulation values found:")
  print(unique(tf_mor_raw$mor))
  
  tf_mor_summary <- tf_mor_raw %>%
    mutate(
      mor_type = case_when(
        mor == 1 ~ "activator",
        mor == -1 ~ "repressor",
        TRUE ~ "unknown"
      )
    ) %>%
    group_by(source, mor_type) %>%
    summarise(n = sum(n), .groups = 'drop') %>%
    pivot_wider(
      names_from = mor_type,
      values_from = n,
      values_fill = 0
    ) %>%
    rename(gene = source)
  
  message("✓ Mode of regulation annotated")
  message("  Columns: ", paste(colnames(tf_mor_summary), collapse = ", "))
  
  ## 4C: Merge and annotate
  validated_tfs_annotated <- validated_tfs %>%
    left_join(tf_confidence, by = "gene") %>%
    left_join(tf_mor_summary, by = "gene")
  
  validated_tfs_annotated <- validated_tfs_annotated %>%
    mutate(
      validation_strength = case_when(
        max_padj < 0.001 ~ "Very strong",
        max_padj < 0.01 ~ "Strong",
        max_padj < 0.05 ~ "Significant",
        TRUE ~ "Marginal"
      ),
      expression_change = case_when(
        mean_fcMedian > 0.5 ~ "Upregulated",
        mean_fcMedian < -0.5 ~ "Downregulated",
        TRUE ~ "Modest change"
      )
    )
  
  # Add regulation_type safely
  if ("activator" %in% names(validated_tfs_annotated) && 
      "repressor" %in% names(validated_tfs_annotated)) {
    validated_tfs_annotated <- validated_tfs_annotated %>%
      mutate(
        regulation_type = case_when(
          activator > 0 & repressor == 0 ~ "Activator",
          repressor > 0 & activator == 0 ~ "Repressor",
          activator > 0 & repressor > 0 ~ "Dual",
          TRUE ~ "Unknown"
        )
      )
  } else if ("activator" %in% names(validated_tfs_annotated)) {
    validated_tfs_annotated$regulation_type <- "Activator"
  } else if ("repressor" %in% names(validated_tfs_annotated)) {
    validated_tfs_annotated$regulation_type <- "Repressor"
  } else {
    validated_tfs_annotated$regulation_type <- "Unknown"
  }
  
  validated_tfs_annotated <- validated_tfs_annotated %>%
    arrange(max_padj)
  
  message("✓ Complete annotation finished")
  
  # ============================================================
  # Step 5: Display Results
  # ============================================================
  
  message("\n=== Step 5: Validated TF Summary ===\n")
  
  summary_display <- validated_tfs_annotated %>%
    select(gene, max_padj, mean_fcMedian, 
           consistent_direction, validation_strength, 
           expression_change, confidence_levels, 
           n_targets_total, regulation_type)
  
  print(summary_display)
  
  message("\n=== TF Statistics ===")
  message("Total validated TFs: ", nrow(validated_tfs_annotated))
  
  message("\nBy validation strength:")
  print(table(validated_tfs_annotated$validation_strength))
  
  message("\nBy expression change:")
  print(table(validated_tfs_annotated$expression_change))
  
  message("\nBy regulation type:")
  print(table(validated_tfs_annotated$regulation_type))
  
  message("\nBy highest confidence:")
  print(table(validated_tfs_annotated$highest_confidence))
  
  message("\nBy direction consistency:")
  print(table(validated_tfs_annotated$consistent_direction))
  
  # ============================================================
  # Step 6: Save Results
  # ============================================================
  
  message("\n=== Step 6: Saving Results ===")
  
  write.csv(validated_tfs_annotated, 
            "20251120_Validated_TFs_Complete_Annotation.csv", 
            row.names = FALSE)
  
  top_tfs <- validated_tfs_annotated %>%
    filter(consistent_direction == TRUE) %>%
    arrange(max_padj)
  
  write.csv(top_tfs,
            "20251120_Top_Validated_TFs_Consistent.csv",
            row.names = FALSE)
  
  message("\n✓ Saved:")
  message("  - 20251120_Validated_TFs_Complete_Annotation.csv")
  message("  - 20251120_Top_Validated_TFs_Consistent.csv")
}

# End of script
message("\n=== Analysis Complete ===")