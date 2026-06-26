# Code for: Single-cell multi-omic analysis of the mouse heart (scRNA-seq + scATAC-seq)

This repository contains all analysis and plotting scripts used in the study.
All scripts are written in **R** (notebooks run on the **IRkernel**).

## Directory structure

```
code/
├── analysis/                      # data processing & downstream analysis
│   ├── 1.Load_QC_scRNA.R                # load 10x Multiome, QC, doublet/decontX filtering, RNA+ATAC integration
│   ├── 2.scRNA+scATAC_annalysis.ipynb   # clustering, cell-type annotation, gene-set / motif activity (AUCell)
│   ├── 3.vCManalysis.ipynb              # ventricular cardiomyocyte (vCM) sub-clustering & scDist analysis
│   ├── 4.figR.ipynb                     # FigR regulatory network inference (DORC + GRN)
│   ├── 4.figR copy.ipynb                # (backup of 4.figR)
│   ├── 5.cell_subtype_analysis.ipynb    # subtype abundance / prevalence filtering
│   └── 6.sup_table.ipynb                # supplementary tables
└── plot/                          # figure & supplementary-figure plotting
    ├── 1.fig1plot.ipynb … 5.fig5plot.ipynb   # main figures 1–5
    └── 6.S1plot.ipynb … 10.S5plot.ipynb      # supplementary figures S1–S5
```

> Inputs are read from a relative `data/` folder (e.g. `data/mHeart.rds`,
> `data/vCM.Rds`) and figures/CSVs are written under `plot/` and `results/`.
> Adjust these paths to your local layout if needed.

---

## R environment (tested with R 4.3.x, IRkernel)

Recommended: **R ≥ 4.3.0** on Linux (the integration steps are multi-core and
memory intensive). Install [IRkernel](https://github.com/IRkernel/IRkernel) to
run the notebooks.

### CRAN packages

| Package | Tested version | Package | Tested version |
|---|---|---|---|
| Seurat | 5.1.0 | dplyr | 1.1.4 |
| Signac | 1.13.0 | tidyr | 1.3.1 |
| harmony | 0.1.7 | tidyverse | 2.0.0 |
| ggplot2 | 3.5.1 | reshape2 | 1.4.4 |
| patchwork | 1.3.0 | RColorBrewer | 1.1.3 |
| ggrepel | 0.9.6 | scales | 1.3.0 |
| ggraph | 2.2.1 | igraph | 2.0.3 |
| ggpubr | 0.6.0 | tidygraph | 1.3.1 |
| Hmisc | 5.2.1 | FNN | 1.1.4 |
| openxlsx | 4.2.7.1 | writexl | 1.5.0 |
| readxl | 1.4.3 | outliers | 0.15 |
| future | 1.34.0 | future.apply | 1.11.3 |
| doParallel | 1.0.17 | parallel | (base) |

### Bioconductor packages (Bioc 3.18, matching R 4.3)

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install(c(
  "GenomeInfoDb", "org.Mm.eg.db",
  "EnsDb.Mmusculus.v79", "BSgenome.Mmusculus.UCSC.mm10", "BSgenome",
  "SingleCellExperiment", "scater", "clusterProfiler",
  "chromVAR",          # installs motifmatchr, TFBSTools automatically
  "decontX",           # ambient-RNA decontamination
  "scDblFinder",       # doublet detection
  "BiocParallel", "ComplexHeatmap", "circlize"
))
```

| Package | Tested version | Note |
|---|---|---|
| GenomeInfoDb | 1.40.1 | |
| org.Mm.eg.db | 3.18.0 | mouse annotation |
| EnsDb.Mmusculus.v79 | 2.99.0 | gene model annotation |
| BSgenome.Mmusculus.UCSC.mm10 | 1.4.3 | mm10 genome |
| SingleCellExperiment | 1.24.0 | |
| scater | 1.32.1 | |
| clusterProfiler | 4.12.0 | |
| chromVAR | 1.24.0 | TF deviation (called via Signac); brings motifmatchr/TFBSTools |
| decontX | 1.16.0 | ambient-RNA decontamination |
| scDblFinder | 1.18.0 | doublet detection |
| ComplexHeatmap | 2.20.0 | |
| circlize | 0.4.16 | |
| BiocParallel | 1.38.0 | |

### GitHub-only packages

These are **not on CRAN/Bioc** and must be installed from source with
`devtools`:

```r
install.packages("devtools")

# scDist — distribution of cell-type-level shifts
devtools::install_github("phillipnicol/scDist")

# FigR — joint scRNA/scATAC regulatory network inference
devtools::install_github("buenrostrolab/FigR")

# irGSEA — integrated gene-set enrichment on single cells
devtools::install_github("chuiqin/irGSEA")

# scop — single-cell omics plotting helpers
devtools::install_github("mengxu98/scop")

# ROGUE — assessing cluster purity
devtools::install_github("PaulingLiu/ROGUE")

# miloR — differential abundance on KNN graphs
devtools::install_github("MarioniLab/miloR")

# ggunchained / BuenColors — colour palettes
devtools::install_github("doswald/ggunchained")
devtools::install_github("caleblareau/BuenColors")
```

| Package | Source | Purpose |
|---|---|---|
| scDist | `phillipnicol/scDist` | cell-type shift distance |
| FigR | `buenrostrolab/FigR` | GRN / DORC inference |
| irGSEA | `chuiqin/irGSEA` | gene-set enrichment (AUCell etc.) |
| scop | `mengxu98/scop` | plotting utilities (`GroupHeatmap`, `FeatureDimPlot`) |
| ROGUE | `PaulingLiu/ROGUE` | cluster purity |
| miloR | `MarioniLab/miloR` | differential abundance |
| ggunchained | `doswald/ggunchained` | colour palettes |
| BuenColors | `caleblareau/BuenColors` | colour palettes |

### Quick version export

After setting up the environment, dump the exact versions used to run the code:

```r
pkgs <- c("Seurat","Signac","harmony","chromVAR","decontX","scDblFinder",
          "ComplexHeatmap","clusterProfiler","FigR","irGSEA","scop",
          "miloR","ROGUE","scDist","ggplot2","dplyr")
write.csv(data.frame(Package = pkgs,
                     Version = sapply(pkgs, function(p) as.character(packageVersion(p)))),
          "session_versions.csv", row.names = FALSE)
```

---

## Python environment

> **Note:** the code in this repository is written entirely in R; no `.py`
> scripts or Python-kernel notebooks are used. A Python environment is **only**
> required for the upstream 10x Genomics pipeline that produced the inputs
> (i.e. **Cell Ranger ARC**), which is run on the raw FASTQ before these scripts.

| Tool | Tested version | Purpose |
|---|---|---|
| Python | 3.10+ | runtime for Cell Ranger ARC |
| Cell Ranger ARC | 2.0.x | 10x Multiome (GEX + ATAC) alignment & count |

If you start from already-processed `filtered_feature_bc_matrix.h5` /
`atac_fragments.tsv.gz` / `atac_peaks.bed` outputs (as this code does), no
Python environment is needed at all.

---

## Reproducibility

1. Place the processed Multiome outputs and the saved Seurat objects under
   `data/`.
2. Open the notebooks in `analysis/` first (in numerical order), then `plot/`.
3. Run on a machine with sufficient RAM (the merged-object integration benefits
   from ≥ 64 GB RAM) and multiple cores (`plan("multicore", ...)`).

For any questions regarding the analysis, please refer to the corresponding
author of the manuscript.
