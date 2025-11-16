# ============================================================================
# 03_SynGO_analysis.R
# Description: SynGO gene ontology enrichment with Countsplit validation
#              - AUCell scoring on synaptic gene sets
#              - Split A (selection) → Split B (testing)
#              - Statistical testing with permutation and bootstrap
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(AUCell)
  library(Matrix)
  library(ggplot2)
  library(rstatix)
  library(effsize)
  library(ggrepel)
})

SEED <- 123
set.seed(SEED)

message("=== SynGO AUCell Analysis with Countsplit ===\n")

# ============================================================================
# PART 1: Load Data and Prepare
# ============================================================================

message("=== PART 1: Data Preparation ===\n")

# Load trajectory data
sce_final <- readRDS("data/sce_trajectory.rds")
message("Loaded SingleCellExperiment: ", ncol(sce_final), " cells")

# Convert to Seurat if needed for easier manipulation
seurat_obj <- as.Seurat(sce_final, data = NULL)

# Define cell groups (Neurons, Immature neurons, Myofibroblasts)
message("\nDefining cell groups...")

cell_types <- seurat_obj$celltype_merged
is_neuron <- cell_types %in% c("Neurons")
is_imm <- cell_types %in% c("Immature neurons")
is_myo <- cell_types %in% c("Myofibroblasts")

seurat_obj$group3 <- ifelse(is_neuron, "Neurons",
                            ifelse(is_imm, "Immature neurons",
                                   ifelse(is_myo, "Myofibroblasts", NA_character_)))

# Keep only cells in the three groups
keep_cells <- seurat_obj$group3 %in% c("Neurons", "Immature neurons", "Myofibroblasts")
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_cells])

seurat_obj$group3 <- factor(seurat_obj$group3,
                            levels = c("Neurons", "Immature neurons", "Myofibroblasts"))

message("Cells retained: ", ncol(seurat_obj))
message("Group composition:")
print(table(seurat_obj$group3))

# ============================================================================
# PART 2: Count Splitting
# ============================================================================

message("\n=== PART 2: Count Splitting ===\n")

# Get count matrix
counts_original <- GetAssayData(seurat_obj, layer = "counts")
if (is.null(counts_original)) {
  counts_original <- GetAssayData(seurat_obj, slot = "counts")
}

message("Original counts: ", nrow(counts_original), " genes × ", ncol(counts_original), " cells")

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
set.seed(SEED)
splits <- split_binomial_sparse(counts_original, p = 0.5, seed = SEED)
counts_A <- splits$A
counts_B <- splits$B

message("Split A: ", sum(counts_A), " total counts")
message("Split B: ", sum(counts_B), " total counts")
message("Split ratio: ", round(sum(counts_A)/(sum(counts_A)+sum(counts_B)), 3))

# ============================================================================
# PART 3: Helper Functions
# ============================================================================

message("\n=== PART 3: Loading Helper Functions ===\n")

# Function to find column in dataframe
.syngo_pick_col <- function(df, candidates){
  here <- tolower(names(df))
  hit <- match(tolower(candidates), here)
  nm <- names(df)[hit[!is.na(hit)]]
  if (length(nm) == 0) stop("Required column missing")
  nm[1]
}

# Map gene symbols to object genes
.map_to_obj <- function(sym_vec, obj_genes){
  obj_genes[match(toupper(sym_vec), toupper(obj_genes), nomatch = 0)]
}

# Jaccard similarity
.jaccard <- function(a, b){
  u <- union(a, b)
  if (length(u) == 0) return(0)
  length(intersect(a, b)) / length(u)
}

# Sequence helper
.seq2 <- function(a, b) if (a <= b) seq.int(a, b) else integer(0)

# Build pruned SynGO gene sets
build_pruned_syngo_sets <- function(syngo_xlsx, obj_genes,
                                    min_size = 10, max_size = 1000,
                                    jaccard_cut = 0.7){
  
  message("Reading SynGO annotations from: ", syngo_xlsx)
  syngo_ann <- readxl::read_xlsx(syngo_xlsx)
  
  term_col <- .syngo_pick_col(syngo_ann, c("go_name", "name", "term", "term_name"))
  gene_col <- .syngo_pick_col(syngo_ann, c("hgnc_symbol", "gene_symbol", "symbol"))
  
  syngo_ann2 <- syngo_ann %>%
    dplyr::mutate(term = .data[[term_col]], gene = .data[[gene_col]]) %>%
    dplyr::filter(!is.na(term), !is.na(gene),
                  nzchar(as.character(term)), nzchar(as.character(gene)))
  
  # Build gene sets
  sets_all <- syngo_ann2 %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(genes = list(unique(.map_to_obj(gene, obj_genes))),
                     .groups = "drop") %>%
    dplyr::mutate(n = lengths(genes)) %>%
    dplyr::filter(n >= min_size, n <= max_size, n > 0)
  
  sets <- setNames(sets_all$genes, sets_all$term)
  
  # Remove redundant sets based on Jaccard similarity
  nm <- names(sets)
  L <- vapply(sets, length, integer(1))
  keep <- rep(TRUE, length(nm))
  
  for (i in seq_along(nm)){
    if (!keep[i]) next
    for (j in .seq2(i + 1, length(nm))){
      if (!keep[j]) next
      if (.jaccard(sets[[i]], sets[[j]]) > jaccard_cut){
        if (L[i] >= L[j]) keep[j] <- FALSE else { keep[i] <- FALSE; break }
      }
    }
  }
  
  message("Original terms: ", length(sets))
  message("After pruning: ", sum(keep))
  
  sets[keep]
}

# AUCell scoring function
aucell_scores <- function(count_matrix, sets_list, rankings = NULL,
                          auc_top_frac = 0.05, auc_cap_frac = 0.2){
  
  # Normalize counts to log-scale
  expr <- log1p(t(t(count_matrix) / colSums(count_matrix)) * 10000)
  
  if (is.null(rankings)) rankings <- AUCell_buildRankings(expr, plotStats = FALSE)
  
  max_size <- max(vapply(sets_list, length, 1L))
  aucMaxRank <- min(max(ceiling(1.1 * max_size), ceiling(auc_top_frac * nrow(rankings))),
                    floor(auc_cap_frac * nrow(rankings)))
  
  auc <- AUCell_calcAUC(sets_list, rankings, aucMaxRank = aucMaxRank)
  
  list(auc_mat = t(assay(auc)), rankings = rankings, aucMaxRank = aucMaxRank)
}

# Statistical testing function
term_stats <- function(scores, groups, nperm = 2000, boot_n = 200, boot_B = 300){
  
  g <- factor(groups, levels = c("Neurons", "Immature neurons", "Myofibroblasts"))
  
  # Medians
  med <- tapply(scores, g, median, na.rm = TRUE)
  
  # Cliff's Delta effect sizes
  dNI <- effsize::cliff.delta(scores[g == "Neurons"],
                              scores[g == "Immature neurons"])$estimate
  dIM <- effsize::cliff.delta(scores[g == "Immature neurons"],
                              scores[g == "Myofibroblasts"])$estimate
  dNM <- effsize::cliff.delta(scores[g == "Neurons"],
                              scores[g == "Myofibroblasts"])$estimate
  
  # Kruskal-Wallis test
  kwp <- tryCatch(rstatix::kruskal_test(scores ~ g)$p, error = function(e) NA_real_)
  
  # Pairwise Wilcoxon tests
  pw <- suppressWarnings(rstatix::pairwise_wilcox_test(
    data.frame(scores, g), scores ~ g, p.adjust.method = "BH"))
  
  getp <- function(a, b){
    i <- which((pw$group1 == a & pw$group2 == b) | (pw$group1 == b & pw$group2 == a))
    if (length(i) > 0) pw$p.adj[i[1]] else NA_real_
  }
  
  pNI <- getp("Neurons", "Immature neurons")
  pIM <- getp("Immature neurons", "Myofibroblasts")
  pNM <- getp("Neurons", "Myofibroblasts")
  
  # Permutation test (Neurons vs Myofibroblasts)
  obs <- med["Neurons"] - med["Myofibroblasts"]
  more <- 0L
  for (k in seq_len(nperm)) {
    lab <- sample(g)
    stat <- tapply(scores, lab, median, na.rm = TRUE)
    if (abs(stat["Neurons"] - stat["Myofibroblasts"]) >= abs(obs)) more <- more + 1L
  }
  pperm <- (more + 1) / (nperm + 1)
  
  # Bootstrap CI
  idx <- split(seq_along(scores), g)
  bn <- min(boot_n, length(idx$Neurons))
  bm <- min(boot_n, length(idx$Myofibroblasts))
  boot <- numeric(boot_B)
  
  for (b in seq_len(boot_B)){
    mN <- median(scores[sample(idx$Neurons, bn, replace = TRUE)], na.rm = TRUE)
    mM <- median(scores[sample(idx$Myofibroblasts, bm, replace = TRUE)], na.rm = TRUE)
    boot[b] <- mN - mM
  }
  ci <- stats::quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  
  list(med = med, dNI = dNI, dIM = dIM, dNM = dNM,
       kwp = kwp, pNI = pNI, pIM = pIM, pNM = pNM,
       pperm = pperm, ciNM = unname(ci))
}

# ============================================================================
# PART 4: Build SynGO Gene Sets
# ============================================================================

message("\n=== PART 4: Building SynGO Gene Sets ===\n")

obj_genes <- rownames(seurat_obj)
sets_pruned <- build_pruned_syngo_sets(
  "data/syngo_annotations.xlsx",
  obj_genes,
  min_size = 10,
  max_size = 1000,
  jaccard_cut = 0.7
)

message("Pruned SynGO terms: ", length(sets_pruned))
stopifnot(length(sets_pruned) > 0)

# ============================================================================
# PART 5A: SELECTION using Split A
# ============================================================================

message("\n=== PART 5A: SELECTION (Split A) ===\n")

# Calculate AUCell scores on Split A
message("Calculating AUCell scores on Split A...")
auc_res_A <- aucell_scores(counts_A, sets_pruned)
auc_mat_A <- auc_res_A$auc_mat

message("AUCell matrix: ", nrow(auc_mat_A), " cells × ", ncol(auc_mat_A), " terms")

# Create dataframe with group labels
df_syn_A <- cbind.data.frame(
  group3 = seurat_obj$group3,
  as.data.frame(auc_mat_A, check.names = FALSE)
)

# Run statistical tests on all terms (Split A)
message("Running statistical tests on Split A...")

one_term_A <- function(term){
  st <- term_stats(df_syn_A[[term]], df_syn_A$group3,
                   nperm = 2000, boot_n = 200, boot_B = 300)
  
  data.frame(
    term = term,
    n_genes = length(sets_pruned[[term]]),
    med_N = st$med["Neurons"],
    med_I = st$med["Immature neurons"],
    med_M = st$med["Myofibroblasts"],
    Cliff_NI = st$dNI,
    Cliff_IM = st$dIM,
    Cliff_NM = st$dNM,
    KW_p = st$kwp,
    p_NI = st$pNI,
    p_IM = st$pIM,
    p_NM = st$pNM,
    perm_p_NM = st$pperm,
    boot_CI_NM_low = st$ciNM[1],
    boot_CI_NM_high = st$ciNM[2],
    stringsAsFactors = FALSE
  )
}

tab_A <- do.call(rbind, lapply(colnames(auc_mat_A), one_term_A))

# Multiple testing correction
tab_A$KW_q <- p.adjust(tab_A$KW_p, method = "BH")
tab_A$q_NI <- p.adjust(tab_A$p_NI, method = "BH")
tab_A$q_IM <- p.adjust(tab_A$p_IM, method = "BH")
tab_A$q_NM <- p.adjust(tab_A$p_NM, method = "BH")

# Determine order (N > I > M, etc.)
tab_A$order <- apply(tab_A[, c("med_N", "med_I", "med_M")], 1, function(v){
  nm <- c("N", "I", "M")[order(v, decreasing = TRUE)]
  paste0(nm, collapse = "")
})

# Selection threshold (liberal for discovery)
SELECTION_P <- 0.1
sig_A <- (tab_A$q_NM < SELECTION_P) | (tab_A$perm_p_NM < SELECTION_P)

message("Terms selected (p < ", SELECTION_P, "): ", sum(sig_A))

# Save selection results
write.csv(tab_A, "results/SynGO_splitA_selection.csv", row.names = FALSE)
message("✓ Saved: results/SynGO_splitA_selection.csv")

# ============================================================================
# PART 5B: TESTING using Split B
# ============================================================================

message("\n=== PART 5B: TESTING (Split B) ===\n")

# Calculate AUCell scores on Split B
message("Calculating AUCell scores on Split B...")
auc_res_B <- aucell_scores(counts_B, sets_pruned)
auc_mat_B <- auc_res_B$auc_mat

# Create dataframe
df_syn_B <- cbind.data.frame(
  group3 = seurat_obj$group3,
  as.data.frame(auc_mat_B, check.names = FALSE)
)

# Test only selected terms (Split B)
selected_terms <- tab_A$term[sig_A]
message("Testing ", length(selected_terms), " selected terms on Split B...")

one_term_B <- function(term){
  st <- term_stats(df_syn_B[[term]], df_syn_B$group3,
                   nperm = 2000, boot_n = 200, boot_B = 300)
  
  data.frame(
    term = term,
    n_genes = length(sets_pruned[[term]]),
    med_N = st$med["Neurons"],
    med_I = st$med["Immature neurons"],
    med_M = st$med["Myofibroblasts"],
    Cliff_NI = st$dNI,
    Cliff_IM = st$dIM,
    Cliff_NM = st$dNM,
    KW_p = st$kwp,
    p_NI = st$pNI,
    p_IM = st$pIM,
    p_NM = st$pNM,
    perm_p_NM = st$pperm,
    boot_CI_NM_low = st$ciNM[1],
    boot_CI_NM_high = st$ciNM[2],
    stringsAsFactors = FALSE
  )
}

tab_B <- do.call(rbind, lapply(selected_terms, one_term_B))

# Multiple testing correction
tab_B$KW_q <- p.adjust(tab_B$KW_p, method = "BH")
tab_B$q_NI <- p.adjust(tab_B$p_NI, method = "BH")
tab_B$q_IM <- p.adjust(tab_B$p_IM, method = "BH")
tab_B$q_NM <- p.adjust(tab_B$p_NM, method = "BH")

tab_B$order <- apply(tab_B[, c("med_N", "med_I", "med_M")], 1, function(v){
  nm <- c("N", "I", "M")[order(v, decreasing = TRUE)]
  paste0(nm, collapse = "")
})

# ============================================================================
# PART 6: Final Validated Results
# ============================================================================

message("\n=== PART 6: Final Validated Results ===\n")

# Final: significant in BOTH splits
FINAL_FDR <- 0.05
sig_B <- (tab_B$q_NM < FINAL_FDR) | (tab_B$perm_p_NM < FINAL_FDR)

tab_final <- tab_B[sig_B, ]
tab_final <- tab_final[order(tab_final$q_NM), ]

message("Validated terms (FDR < ", FINAL_FDR, "): ", nrow(tab_final))

# Save results
write.csv(tab_B, "results/SynGO_splitB_testing_all.csv", row.names = FALSE)
write.csv(tab_final, "results/SynGO_splitB_testing_validated.csv", row.names = FALSE)

message("✓ Saved: results/SynGO_splitB_testing_all.csv")
message("✓ Saved: results/SynGO_splitB_testing_validated.csv")

# ============================================================================
# PART 7: Visualization (Volcano Plot)
# ============================================================================

message("\n=== PART 7: Creating Volcano Plot ===\n")

dfp <- tab_B
dfp$neglog10q <- -log10(pmax(dfp$q_NM, 1e-300))
dfp$class <- ifelse(dfp$order %in% c("NIM", "NMI") & sig_B, "Neuron-enriched",
                    ifelse(dfp$order %in% c("MIN", "MNI") & sig_B, "Myofibroblast-enriched",
                           ifelse(dfp$order %in% c("INM", "IMN") & sig_B, "Immature-enriched",
                                  "Not significant")))

# Select top terms for labeling
lab_candidates <- dfp %>%
  dplyr::filter(class != "Not significant") %>%
  dplyr::mutate(score = abs(Cliff_NM) * neglog10q) %>%
  dplyr::group_by(class) %>%
  dplyr::slice_max(score, n = 8, with_ties = FALSE) %>%
  dplyr::ungroup()

# Create volcano plot
p_volcano <- ggplot(dfp, aes(x = Cliff_NM, y = neglog10q)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = class), alpha = 0.7, size = 2.5) +
  ggrepel::geom_text_repel(
    data = lab_candidates,
    aes(label = term, color = class),
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    size = 3,
    seed = SEED
  ) +
  scale_color_manual(
    values = c(
      "Neuron-enriched" = "#1f77b4",
      "Myofibroblast-enriched" = "#d62728",
      "Immature-enriched" = "#2ca02c",
      "Not significant" = "grey60"
    ),
    name = "Classification"
  ) +
  labs(
    x = "Cliff's Δ (Neurons - Myofibroblasts)",
    y = "-log10(FDR)",
    title = "SynGO Analysis with Countsplit",
    subtitle = "Selection: Split A (p < 0.1) | Testing: Split B (FDR < 0.05)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

ggsave("figures/SynGO_countsplit_volcano.svg", p_volcano,
       width = 12, height = 9, device = "svg")

message("✓ Saved: figures/SynGO_countsplit_volcano.svg")

# ============================================================================
# Summary
# ============================================================================

message("\n=== SUMMARY ===")

n_neuron <- sum(tab_B$order %in% c("NIM", "NMI") & sig_B)
n_myo <- sum(tab_B$order %in% c("MIN", "MNI") & sig_B)
n_imm <- sum(tab_B$order %in% c("INM", "IMN") & sig_B)

message("Selection phase (Split A, p < 0.1): ", sum(sig_A), " candidates")
message("Testing phase (Split B, FDR < 0.05): ", nrow(tab_final), " validated")
message("\nValidated terms by category:")
message("  Neuron-enriched: ", n_neuron)
message("  Myofibroblast-enriched: ", n_myo)
message("  Immature-enriched: ", n_imm)

if (nrow(tab_final) > 0) {
  message("\nTop 5 validated neuron-enriched terms:")
  neuronal_top <- tab_final[tab_final$order %in% c("NIM", "NMI"), ]
  if (nrow(neuronal_top) > 0) {
    neuronal_top <- neuronal_top[order(-neuronal_top$Cliff_NM), ]
    print(head(neuronal_top[, c("term", "n_genes", "Cliff_NM", "q_NM", "perm_p_NM")], 5))
  }
}

message("\n✓ SynGO analysis with countsplit completed!")
message("\nFiles saved:")
message("  - results/SynGO_splitA_selection.csv")
message("  - results/SynGO_splitB_testing_all.csv")
message("  - results/SynGO_splitB_testing_validated.csv")
message("  - figures/SynGO_countsplit_volcano.svg")