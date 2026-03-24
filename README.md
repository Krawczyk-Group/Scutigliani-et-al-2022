# The Effects of Heat Stress on the Transcriptome of Human Cancer Cells

[![Paper: Cancers](https://img.shields.io/badge/Journal-Cancers-blue.svg)](https://www.mdpi.com/2072-6694/15/1/113)
[![DOI: 10.3390/cancers15010113](https://img.shields.io/badge/DOI-10.3390%2Fcancers15010113-brightgreen.svg)](https://doi.org/10.3390/cancers15010113)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

This repository contains the curated datasets, analysis scripts, and reproducible workflows for the meta-analysis published in **Cancers**: 

> **"The Effects of Heat Stress on the Transcriptome of Human Cancer Cells: A Meta-Analysis"** > *Enzo M. Scutigliani, Fernando Lobo-Cerna, Sergio Mingo Barba, Stephan Scheidegger, and Przemek M. Krawczyk.*

---

## Table of Contents
- [Project Overview](#project-overview)
- [Key Findings](#key-findings)
- [Repository Structure](#repository-structure)
- [Data Sources](#data-sources)
- [Getting Started](#getting-started)
- [Citation](#citation)

---

## Project Overview

Hyperthermia is a clinically applied cancer treatment used in conjunction with radio- and/or chemotherapy, wherein the tumor volume is exposed to supraphysiological temperatures (39–44 °C). This project provides a comprehensive meta-analysis of **18 transcriptomic datasets** across **9 human cancer cell lines**, aiming to delineate the core molecular shifts and the inherent variability of cellular responses to heat stress.

We analyzed global gene expression profiles to understand how hyperthermia dysregulates critical biological processes, including protein folding, cell cycle control, and programmed cell death.

---

## Key Findings

* **Pathway Reshaping:** Hyperthermia alters the expression of genes involved in **KRAS**, **TNF-alpha**, and **EMT** (Epithelial-to-Mesenchymal Transition) signaling.
* **Transcriptional Complexity:** Our results demonstrate a considerable inter-study variability and an apparent absence of a "universal" heat stress signature; response patterns are highly dependent on experimental context, including:
    * **Temperature** (°C) and **Duration** (minutes)
    * **Time of recovery** after treatment
    * **Cell line specificity**
* **Biological Impact:** Consistent alterations were observed in protein folding, cell cycle division, and mitosis. The small heat-shock protein CRYAB was found to be the most commonly overexpressed gene across datasets.

---

## Repository Structure

```text
├── .github/                # GitHub Actions & workflows
├── Scutigliani 2022a/      # Supplementary processed data
├── analysis.R              # Main meta-analysis processing script
├── data_wrangling_...R     # Data cleaning and merging script
├── *.csv                   # Processed datasets (18 total)
├── paper.pdf               # Full-text PDF of the publication
└── README.md               # Documentation (you are here)
---

## Data Sources

The analysis integrates data from multiple landmark studies in the field:

| Study | Datasets | Genes/Conditions |
| :--- | :--- | :--- |
| **Scutigliani 2022** | a, b | Internal study data |
| **Amaya 2014** | a, b, c | HT in breast cancer cell lines |
| **Andocs 2015** | a, b | In vitro HT effects |
| **Court 2017** | 1 | HT treatment profiles |
| **Tabuchi 2011** | a - h | Large-scale HT screening |
| **Yunoki 2016** | 1 | Time-resolved HT response |

---

## Getting Started

### Prerequisites
To reproduce the analysis, you will need **R (>= 4.0.0)** and the following libraries:

```R
# Visualization
install.packages(c("ggplot2", "tidyverse", "cowplot", "ggrepel", "RColorBrewer"))

# Bioinformatics & Pathway Analysis
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "msigdbr", "org.Hs.eg.db", "VennDiagram", "UpSetR", "ggupset", "pheatmap"))

# Statistics
install.packages(c("matrixTests", "ggfortify"))
```

### Reproducing the Analysis

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/Krawczyk-Group/Scutigliani-et-al-2022.git](https://github.com/Krawczyk-Group/Scutigliani-et-al-2022.git)
      cd Scutigliani-et-al-2022
    ```
2.  **Run Data Wrangling:**
    Execute `data_wrangling_before_analysis.R` to consolidate raw study files into the master expression matrix (`expression_all.csv`).
3.  **Run Main Analysis:**
    Execute `analysis.R` to generate statistical summaries, PCA plots, Heatmaps, and UpSet visualizations.

---

## 📑 Citation

If you use this repository or these datasets in your research, please cite the original publication:

> **Scutigliani, E. M., Lobo-Cerna, F., Mingo Barba, S., Scheidegger, S., & Krawczyk, P. M. (2022).** The Effects of Heat Stress on the Transcriptome of Human Cancer Cells: A Meta-Analysis. *Cancers (Basel)*, 14(24), 6022. [https://doi.org/10.3390/cancers14246022](https://doi.org/10.3390/cancers14246022)

---
*Maintained by the **P. Krawczyk Research Group**. For queries or collaboration, please open an issue.*
