# scRNA-seq Direct Neuronal Reprogramming Analysis

Analysis pipeline for PTBP1 knockdown-mediated direct neuronal reprogramming in human dermal fibroblasts.

**Data Availability:** NCBI BioProject [PRJNA1256192](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1256192)

## Repository Structure

```
scRNA-seq-reprogramming-analysis/
├── scripts/
│   ├── utils/                                      # ScType helper functions
│   ├── 01_preprocessing_clustering_trajectory.R    # QC, clustering, trajectory
│   ├── 02_DEG_countsplit.R                         # DEG with Countsplit validation
│   ├── 03_SynGO_analysis.R                         # Synaptic gene enrichment
│   └── 04_figures.R                                # Heatmap, dotplots, TF plots
└── data/
    ├── ScTypeDB_full.xlsx                          # Cell type database
    ├── syngo_annotations.xlsx                      # SynGO gene sets
    └── sce_trajectory.rds                          # Preprocessed data (8.4 MB)
```

## Analysis Pipeline

### 1. Preprocessing & Trajectory (Script 01)
- Quality control, normalization
- ScType cell type annotation
- Slingshot trajectory inference
- **Output:** `sce_trajectory.rds`

**Note:** Requires raw data from NCBI (PRJNA1256192). Reviewers can skip and use provided `sce_trajectory.rds`.

### 2. Differential Expression (Script 02)
- tradeSeq patternTest (mature vs immature lineages)
- Countsplit validation (Split A selection → Split B testing)
- Transcription factor identification
- **Output:** Validated gene lists and TF annotations

### 3. SynGO Enrichment (Script 03)
- AUCell scoring on synaptic gene sets
- Countsplit validation across cell types
- **Output:** Validated SynGO terms, volcano plot

### 4. Publication Figures (Script 04)
- Gene expression heatmap
- GO/KEGG enrichment dotplots
- Transcription factor smoothers
- **Output:** All figures in SVG format

## Results Summary

**Cell composition (3,702 cells):**
- Fibroblasts: 188 (5.1%)
- Immature neurons: 1,657 (44.8%)
- Myofibroblasts: 1,213 (32.8%)
- Mature neurons: 644 (17.4%)

**Validated findings (FDR < 0.05):**
- Differential expression genes validated through Countsplit
- Key transcription factors: CEBPB, EGR1, PBX3
- Synaptic gene ontology enrichment patterns

## Contact

Dojin Seo
