# ============================================================================
# 03_SynGO_analysis.R
# Description: SynGO synaptic gene enrichment with Countsplit cross-validation
#              - Implements rigorous two-fold cross-validation
#              - Split A → Select genes | Split B → Validate enrichment
#              - Split B → Select genes | Split A → Validate enrichment
#              - Prevents "double dipping" circularity
# Reference: Neufeld et al., Biostatistics (2024)
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
})

SEED <- 123
set.seed(SEED)

message("=== SynGO Cross-Validation Analysis with Countsplit ===\n")

# ============================================================================
# PART 1: Load Data and Prepare Cell Groups
# ============================================================================

message("=== PART 1: Data Preparation ===\n")

# Load preprocessed Seurat object
seurat_obj <- readRDS("data/seurat_preprocessed.rds")
message("Loaded: ", nrow(seurat_obj), " genes × ", ncol(seurat_obj), " cells")

# Define three-group comparison (Neurons, Immature neurons, Myofibroblasts)
message("\nDefining cell groups...")

cell_types <- seurat_obj$customclassif
is_neuron <- cell_types %in% c("GABAergic neurons", "Glutamatergic neurons")
is_imm <- cell_types %in% c("Immature neurons")
is_myo <- cell_types %in% c("Myofibroblasts")

seurat_obj$group3 <- ifelse(is_neuron, "Neurons",
                            ifelse(is_imm, "Immature neurons",
                                   ifelse(is_myo, "Myofibroblasts", NA_character_)))

# Retain only cells in the three groups
keep_cells <- seurat_obj$group3 %in% c("Neurons", "Immature neurons", "Myofibroblasts")
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_cells])

seurat_obj$group3 <- factor(seurat_obj$group3,
                            levels = c("Neurons", "Immature neurons", "Myofibroblasts"))

message("Cells retained: ", ncol(seurat_obj))
message("Group composition:")
print(table(seurat_obj$group3))

# ============================================================================
# PART 2: Count Splitting for Independent Validation
# ============================================================================

message("\n=== PART 2: Count Splitting ===\n")

# Extract original count matrix
counts_original <- GetAssayData(seurat_obj, layer = "counts")
if (is.null(counts_original)) {
  counts_original <- GetAssayData(seurat_obj, slot = "counts")
}

message("Original counts: ", nrow(counts_original), " genes × ", ncol(counts_original), " cells")

# Binomial splitting function
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
# PART 3: Create Separate Seurat Objects for Each Split
# ============================================================================

message("\n=== PART 3: Creating Split-Specific Seurat Objects ===\n")

# Copy metadata from original object
meta_data <- seurat_obj@meta.data

# Split A object
seurat_A <- CreateSeuratObject(counts = counts_A, meta.data = meta_data)
seurat_A <- NormalizeData(seurat_A, verbose = FALSE)
Idents(seurat_A) <- "group3"

# Split B object
seurat_B <- CreateSeuratObject(counts = counts_B, meta.data = meta_data)
seurat_B <- NormalizeData(seurat_B, verbose = FALSE)
Idents(seurat_B) <- "group3"

message("Split A object created: ", ncol(seurat_A), " cells")
message("Split B object created: ", ncol(seurat_B), " cells")

# ============================================================================
# PART 4: Load SynGO Gene Ontology Database
# ============================================================================

message("\n=== PART 4: Loading SynGO Database ===\n")

syngo_df <- read_xlsx("syngo_annotations.xlsx")
syngo_genes_all <- unique(na.omit(syngo_df$hgnc_symbol))

message("Total SynGO genes available: ", length(syngo_genes_all))

# ============================================================================
# PART 5: Cross-Validation Round 1 (Select on A, Validate on B)
# ============================================================================

message("\n=== PART 5: Cross-Validation Round 1 (A→B) ===\n")
message("Selection on Split A | Validation on Split B\n")

# --- 5.1: Gene Selection on Split A ---
message("Step 1: Identifying neuron-specific genes on Split A...")

markers_A <- FindMarkers(seurat_A,
                         ident.1 = "Neurons",
                         ident.2 = "Myofibroblasts",
                         group.by = "group3",
                         logfc.threshold = 0.2,
                         min.pct = 0.1,
                         verbose = FALSE)

# Select upregulated genes in neurons
neuron_genes_A <- rownames(markers_A[markers_A$avg_log2FC > 0 &
                                       markers_A$p_val_adj < 0.05, ])

# Intersect with SynGO to get neuron-specific synaptic genes
refined_genes_A <- intersect(neuron_genes_A, syngo_genes_all)

message("Selected genes from Split A: ", length(refined_genes_A))

# --- 5.2: AUCell Scoring on Split B (Validation) ---
message("Step 2: Computing AUCell scores on Split B (independent validation)...")

gene_sets_A <- list(SynGO_Neuron_Specific = refined_genes_A)

exprMatrix_B <- GetAssayData(seurat_B, layer = "counts")
rankings_B <- AUCell_buildRankings(exprMatrix_B, plotStats=FALSE, verbose=FALSE)
auc_B <- AUCell_calcAUC(gene_sets_A, rankings_B,
                        aucMaxRank = ceiling(0.05 * nrow(rankings_B)))

scores_B <- as.numeric(getAUC(auc_B)["SynGO_Neuron_Specific", ])
seurat_B <- AddMetaData(seurat_B, metadata = scores_B,
                        col.name = "SynGO_Score_ValidatedB")

# --- 5.3: Statistical Testing on Split B ---
message("Step 3: Statistical testing on Split B...")

stat_test_B <- pairwise.wilcox.test(seurat_B$SynGO_Score_ValidatedB,
                                    seurat_B$group3,
                                    p.adjust.method = "bonferroni")

message("\nRound 1 (A→B) Statistical Results:")
print(stat_test_B)

# ============================================================================
# PART 6: Cross-Validation Round 2 (Select on B, Validate on A)
# ============================================================================

message("\n=== PART 6: Cross-Validation Round 2 (B→A) ===\n")
message("Selection on Split B | Validation on Split A\n")

# --- 6.1: Gene Selection on Split B ---
message("Step 1: Identifying neuron-specific genes on Split B...")

markers_B <- FindMarkers(seurat_B,
                         ident.1 = "Neurons",
                         ident.2 = "Myofibroblasts",
                         group.by = "group3",
                         logfc.threshold = 0.2,
                         min.pct = 0.1,
                         verbose = FALSE)

# Select upregulated genes in neurons
neuron_genes_B <- rownames(markers_B[markers_B$avg_log2FC > 0 &
                                       markers_B$p_val_adj < 0.05, ])

# Intersect with SynGO
refined_genes_B <- intersect(neuron_genes_B, syngo_genes_all)

message("Selected genes from Split B: ", length(refined_genes_B))

# --- 6.2: AUCell Scoring on Split A (Validation) ---
message("Step 2: Computing AUCell scores on Split A (independent validation)...")

gene_sets_B <- list(SynGO_Neuron_Specific = refined_genes_B)

exprMatrix_A <- GetAssayData(seurat_A, layer = "counts")
rankings_A <- AUCell_buildRankings(exprMatrix_A, plotStats=FALSE, verbose=FALSE)
auc_A <- AUCell_calcAUC(gene_sets_B, rankings_A,
                        aucMaxRank = ceiling(0.05 * nrow(rankings_A)))

scores_A <- as.numeric(getAUC(auc_A)["SynGO_Neuron_Specific", ])
seurat_A <- AddMetaData(seurat_A, metadata = scores_A,
                        col.name = "SynGO_Score_ValidatedA")

# --- 6.3: Statistical Testing on Split A ---
message("Step 3: Statistical testing on Split A...")

stat_test_A <- pairwise.wilcox.test(seurat_A$SynGO_Score_ValidatedA,
                                    seurat_A$group3,
                                    p.adjust.method = "bonferroni")

message("\nRound 2 (B→A) Statistical Results:")
print(stat_test_A)

# ============================================================================
# PART 7: Visualization of Cross-Validated Results
# ============================================================================

message("\n=== PART 7: Generating Visualizations ===\n")

# Define consistent colors
group_colors <- c("Neurons" = "#E41A1C",
                  "Immature neurons" = "#377EB8",
                  "Myofibroblasts" = "#4DAF4A")

# --- Plot 1: Validation on Split B (genes from A) ---
p_validB <- VlnPlot(seurat_B,
                    features = "SynGO_Score_ValidatedB",
                    group.by = "group3",
                    pt.size = 0,
                    cols = group_colors) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  ggtitle("") +
  labs(y = "AUCell Score\n(Validated on Split B)", x = "") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 13))

print(p_validB)
ggsave("syngo_validated_splitB.png", plot = p_validB,
       dpi = 600, height = 6, width = 6)

# --- Plot 2: Validation on Split A (genes from B) ---
p_validA <- VlnPlot(seurat_A,
                    features = "SynGO_Score_ValidatedA",
                    group.by = "group3",
                    pt.size = 0,
                    cols = group_colors) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  ggtitle("") +
  labs(y = "AUCell Score\n(Validated on Split A)", x = "") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 13))

print(p_validA)
ggsave("syngo_validated_splitA.png", plot = p_validA,
       dpi = 600, height = 6, width = 6)