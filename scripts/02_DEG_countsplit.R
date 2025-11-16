# ============================================================================
# 02_DEG_countsplit.R
# Description: Trajectory-based differential expression with Countsplit validation
#              - patternTest: Compare gene expression patterns between lineages
#              - Countsplit: Split A (selection) → Split B (testing)
#              - TF analysis: Identify validated transcription factors
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(Matrix)
  library(dplyr)
  library(decoupleR)
  library(OmnipathR)
  library(ggplot2)
})

SEED <- 123
set.seed(SEED)

# ============================================================================
# PART 1: Load trajectory data
# ============================================================================

message("=== PART 1: Load Trajectory Data ===\n")

sce_final <- readRDS("data/sce_trajectory.rds")

message("Loaded SingleCellExperiment object:")
message("  Cells: ", ncol(sce_final))
message("  Genes: ", nrow(sce_final))

# Extract pseudotime and weights from colData
# Note: Adjust column names based on your sce_trajectory.rds structure
pseudotime_col <- "curve1"  # or the actual pseudotime column name
primary_lineage_col <- "primary_lineage"

if (!pseudotime_col %in% colnames(colData(sce_final))) {
  stop("Pseudotime column '", pseudotime_col, "' not found in colData")
}

# ============================================================================
# PART 2: Identify lineages
# ============================================================================

message("\n=== PART 2: Identify Lineages ===\n")

# Analyze cell type composition per lineage
lineage_composition <- table(
  sce_final$celltype_merged,
  colData(sce_final)[[primary_lineage_col]]
)

message("Cell type composition per lineage:")
print(lineage_composition)

# Identify which lineage is which based on dominant cell type
# Mature lineage: dominated by "Neurons"
# Immature lineage: dominated by "Immature neurons"

lineage_props <- prop.table(lineage_composition, margin = 2)
mature_lin <- which.max(lineage_props["Neurons", ])
immature_lin <- which.max(lineage_props["Immature neurons", ])

message("\nLineage assignments:")
message("  Mature neuron lineage: ", mature_lin)
message("  Immature neuron lineage: ", immature_lin)

# ============================================================================
# PART 3: Count Splitting
# ============================================================================

message("\n=== PART 3: Count Splitting ===\n")

# Binomial split function
split_binomial_sparse <- function(mat, p = 0.5, seed = 123){
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

# Perform split
splits <- split_binomial_sparse(assay(sce_final, "counts"), p = 0.5, seed = SEED)
counts_A <- splits$A
counts_B <- splits$B

message("Split A: ", sum(counts_A), " total counts")
message("Split B: ", sum(counts_B), " total counts")
message("Split ratio: ", round(sum(counts_A)/(sum(counts_A)+sum(counts_B)), 3))

# ============================================================================
# PART 4: Prepare combined lineage data
# ============================================================================

message("\n=== PART 4: Prepare Combined Lineage Data ===\n")

# Get lineage assignments
primary_lineages <- colData(sce_final)[[primary_lineage_col]]

# Cells belonging to either mature or immature lineage
lin1_cells <- primary_lineages == mature_lin
lin3_cells <- primary_lineages == immature_lin
both_cells <- lin1_cells | lin3_cells

message("Lineage ", mature_lin, " (Mature neurons): ", sum(lin1_cells), " cells")
message("Lineage ", immature_lin, " (Immature neurons): ", sum(lin3_cells), " cells")
message("Total cells in both lineages: ", sum(both_cells))

# Create pseudotime and weights matrices for both lineages
n_both <- sum(both_cells)
pseudotime_both <- matrix(0, nrow = n_both, ncol = 2)
weights_both <- matrix(0, nrow = n_both, ncol = 2)

colnames(pseudotime_both) <- c("lineage_mature", "lineage_immature")
colnames(weights_both) <- c("lineage_mature", "lineage_immature")

# Get pseudotime for all cells
all_pseudotime <- colData(sce_final)[[pseudotime_col]]

# Fill in pseudotime and weights
# For cells in mature lineage
lin1_in_both <- which(lin1_cells[both_cells])
pseudotime_both[lin1_in_both, 1] <- all_pseudotime[both_cells][lin1_in_both]
weights_both[lin1_in_both, 1] <- 1

# For cells in immature lineage
lin3_in_both <- which(lin3_cells[both_cells])
pseudotime_both[lin3_in_both, 2] <- all_pseudotime[both_cells][lin3_in_both]
weights_both[lin3_in_both, 2] <- 1

# Verify
message("\nVerification:")
message("  NAs in pseudotime: ", sum(is.na(pseudotime_both)))
message("  NAs in weights: ", sum(is.na(weights_both)))

n_lin1_only <- sum(weights_both[, 1] > 0 & weights_both[, 2] == 0)
n_lin3_only <- sum(weights_both[, 1] == 0 & weights_both[, 2] > 0)

message("\nCell assignment:")
message("  Mature lineage only: ", n_lin1_only)
message("  Immature lineage only: ", n_lin3_only)

# ============================================================================
# PART 5A: SELECTION using Split A
# ============================================================================

message("\n=== PART 5A: SELECTION (Split A) ===\n")

# Extract counts for both lineages from Split A
counts_A_both <- counts_A[, both_cells]

# Gene filtering (keep genes expressed in at least 5% of cells)
min_cells_both <- ceiling(ncol(counts_A_both) * 0.05)
keep_genes_A <- rowSums(counts_A_both > 0) >= min_cells_both
counts_A_both_filtered <- counts_A_both[keep_genes_A, ]

message("Genes retained in Split A: ", nrow(counts_A_both_filtered), "/", nrow(counts_A))
message("Cells: ", ncol(counts_A_both_filtered))

# Fit GAM on Split A
message("\nFitting GAM on Split A...")
set.seed(SEED)
sce_gam_both_A <- fitGAM(
  counts = counts_A_both_filtered,
  pseudotime = pseudotime_both,
  cellWeights = weights_both,
  nknots = 4,
  verbose = FALSE
)
message("✓ fitGAM completed on Split A")

# patternTest on Split A (for candidate selection)
message("\nRunning patternTest on Split A for candidate selection...")
patA <- patternTest(sce_gam_both_A, l2fc = log2(1.5)) %>%
  as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  arrange(pvalue)

message("✓ patternTest completed on Split A")

# Selection threshold (liberal for discovery phase)
PATTERN_SELECT_P <- 0.1
cand_pat <- patA %>% 
  filter(pvalue < PATTERN_SELECT_P) %>% 
  pull(gene)

message("\nCandidate genes selected (p < ", PATTERN_SELECT_P, "): ", length(cand_pat), " genes")

# Save selection results
write.csv(patA, "results/patternTest_splitA_selection.csv", row.names = FALSE)

# ============================================================================
# PART 5B: TESTING using Split B
# ============================================================================

message("\n=== PART 5B: TESTING (Split B) ===\n")

# Extract counts for both lineages from Split B
counts_B_both <- counts_B[, both_cells]

# Use same genes as Split A for consistency
counts_B_both_filtered <- counts_B_both[keep_genes_A, ]

message("Genes retained in Split B: ", nrow(counts_B_both_filtered), "/", nrow(counts_B))
message("Cells: ", ncol(counts_B_both_filtered))

# Fit GAM on Split B
message("\nFitting GAM on Split B...")
set.seed(SEED)
sce_gam_both_B <- fitGAM(
  counts = counts_B_both_filtered,
  pseudotime = pseudotime_both,
  cellWeights = weights_both,
  nknots = 4,
  verbose = FALSE
)
message("✓ fitGAM completed on Split B")

# Save GAM model for downstream analyses
saveRDS(sce_gam_both_B, "data/sce_gam_lineages.rds")
message("✓ Saved GAM model: data/sce_gam_lineages.rds")

# patternTest on Split B (independent testing)
message("\nRunning patternTest on Split B for independent testing...")
patB <- patternTest(sce_gam_both_B, l2fc = log2(1.5)) %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))

message("✓ patternTest completed on Split B")

# ============================================================================
# PART 6: Final Validated Results
# ============================================================================

message("\n=== PART 6: Final Validated Results ===\n")

# Final: genes selected in A AND significant in B
PATTERN_TEST_FDR <- 0.05

pat_final <- patB %>%
  filter(gene %in% cand_pat) %>%
  mutate(padj = p.adjust(pvalue, "BH")) %>%
  arrange(padj)

sig_pat <- pat_final %>% 
  filter(padj < PATTERN_TEST_FDR)

message("Validated pattern genes (FDR < ", PATTERN_TEST_FDR, "):")
message("  Candidates from Split A: ", length(cand_pat))
message("  Validated in Split B: ", nrow(sig_pat), " genes")
message("  Validation rate: ", 
        round(nrow(sig_pat) / length(cand_pat) * 100, 1), "%")

# Save results
write.csv(pat_final, "results/patternTest_splitB_testing_all.csv", row.names = FALSE)
write.csv(sig_pat, "results/patternTest_splitB_testing_validated.csv", row.names = FALSE)

# ============================================================================
# PART 7: Transcription Factor Analysis
# ============================================================================

message("\n=== PART 7: Transcription Factor Analysis ===\n")

# Load DoRothEA database
message("Loading DoRothEA TF database...")
net <- get_dorothea(organism = 'human', levels = c('A', 'B', 'C'))

message("DoRothEA database loaded:")
message("  Total TF-target interactions: ", nrow(net))
message("  Unique TFs: ", length(unique(net$source)))

# Identify validated TFs
all_tfs <- unique(net$source)

validated_tfs <- sig_pat %>%
  filter(gene %in% all_tfs) %>%
  arrange(padj)

message("\nValidated transcription factors:")
message("  Total validated genes (FDR < 0.05): ", nrow(sig_pat))
message("  TFs among validated genes: ", nrow(validated_tfs))

if (nrow(validated_tfs) > 0) {
  message("\nTop validated TFs:")
  print(head(validated_tfs[, c("gene", "pvalue", "padj")], 10))
}

# Annotate with DoRothEA confidence information
if (nrow(validated_tfs) > 0) {
  message("\nAnnotating TF confidence levels...")
  
  # Get confidence level distribution
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
  
  # Get mode of regulation distribution
  tf_mor <- net %>%
    filter(source %in% validated_tfs$gene) %>%
    group_by(source, mor) %>%
    summarise(n = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(
      names_from = mor,
      values_from = n,
      values_fill = 0,
      names_prefix = "mor_"
    ) %>%
    rename(gene = source)
  
  # Merge annotations
  validated_tfs_annotated <- validated_tfs %>%
    left_join(tf_confidence, by = "gene") %>%
    left_join(tf_mor, by = "gene")
  
  message("✓ TF annotation completed")
  
  # Save annotated TF results
  write.csv(validated_tfs_annotated, 
            "results/validated_TFs_annotated.csv", 
            row.names = FALSE)
  
  message("\nAnnotated TF summary:")
  print(validated_tfs_annotated[, c("gene", "padj", "confidence_levels", "n_targets_total")])
}


# ============================================================================
# Summary
# ============================================================================

message("\n=== SUMMARY ===")
message("Selection phase (Split A, p < 0.1): ", length(cand_pat), " candidates")
message("Testing phase (Split B, FDR < 0.05): ", nrow(sig_pat), " validated")
message("Validated transcription factors: ", nrow(validated_tfs))
message("\nKey validated TFs:")
if (nrow(validated_tfs) > 0) {
  top_tfs <- head(validated_tfs$gene, 5)
  for (tf in top_tfs) {
    message("  - ", tf)
  }
}

message("\n✓ Differential expression analysis with Countsplit completed!")
message("\nFiles saved:")
message("  - results/patternTest_splitA_selection.csv")
message("  - results/patternTest_splitB_testing_all.csv")
message("  - results/patternTest_splitB_testing_validated.csv")
message("  - results/validated_TFs_annotated.csv")
message("  - data/sce_gam_lineages.rds")
