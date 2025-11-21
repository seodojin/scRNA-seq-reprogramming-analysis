# scRNA-seq Direct Neuronal Reprogramming Analysis

Analysis pipeline for PTBP1 knockdown-mediated direct neuronal reprogramming in human dermal fibroblasts.

**Data Availability:** NCBI BioProject [PRJNA1256192](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1256192)

## Analysis Summary

### 1. QC & Clustering

-   BD Rhapsody WTA pipeline → count matrix\
-   Filtered cells:
    -   **1,000–4,000 features**
    -   **\<20% MT genes**
-   Final dataset: **9,971 cells**
-   PCA → SNN clustering (resolution 0.2) → UMAP

### 2. Cell Type Annotation

-   Automated annotation using **ScType**\
-   Major populations:
    -   Fibroblasts\
    -   Myofibroblasts\
    -   Immature neurons\
    -   Neurons (GABAergic + Glutamatergic)

### 3. Trajectory Inference

-   Converted to **SingleCellExperiment**
-   Slingshot on UMAP embeddings\
-   Root: Fibroblasts\
-   Terminal: Neurons / Immature neurons / Myofibroblasts

### 4. Countsplit Differential Expression (tradeSeq)

-   Raw counts split into **Split A / Split B** (Binomial p=0.5)\
-   **Path 1:** Trajectory on A → DE on B\
-   **Path 2:** Trajectory on B → DE on A\
-   Final DE genes = intersection\
-   Used `patternTest` (FDR \< 0.05)

### 5. SynGO AUCell Analysis

-   72 non-redundant SynGO terms\
-   AUCell scoring performed **only on the validation split**\
-   Bidirectional Countsplit (A→B, B→A)

### 6. Transcription Factor Analysis

-   DoRothEA levels A–C\
-   Validated TFs include:
    -   **CEBPB, EGR1, PBX3, HMGA1, JUND, FOSL1**
