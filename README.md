# CES1-Mediated Cell-Cell Communication in Colorectal Cancer

## Overview

This project investigates the role of CES1 (Carboxylesterase 1) in shaping intercellular communication within the colorectal cancer (CRC) tumour microenvironment (TME) using single-cell RNA sequencing (scRNA-seq) data.

By integrating transcriptomic analysis with cell–cell communication modelling, the study characterises ligand–receptor interactions and pathway-level signalling networks between tumour and immune cell populations.

---

## Research Objectives

- Analyse CES1 expression across tumour and immune cell subsets
- Infer ligand–receptor mediated communication networks
- Perform pathway-level signalling analysis using CellChat
- Identify upstream ligands and downstream targets using NicheNet
- Explore signalling heterogeneity at subcluster resolution

---

## Methodology

### 1. Data Processing
- Quality control filtering
- Normalisation
- Variable feature selection
- Dimensionality reduction (PCA, UMAP)
- Cell type annotation using Seurat

### 2. Cell–Cell Communication Inference
- CellChat for global communication network modelling
- Pathway-level signalling strength analysis
- Visualisation using circle plots, bubble plots and chord diagrams

### 3. Ligand–Target Modelling
- NicheNet to identify candidate ligands influencing gene expression
- Subcluster-level ligand activity comparison

---

## Repository Structure

scripts/
- 01_understanding_dataset.R
- 02_cellchat_analysis.R
- 03_cellchat_subcluster.R
- 04_nichenet_celltype.R
- 05_nichenet_each_celltype.R
- 06_nichenet_subcluster.R

data/
- Dataset description and access instructions

results/
- Visual outputs and communication summaries

---

## Key Findings

- CES1 expression varies across tumour-associated macrophages and immune subsets.
- Multiple ligand–receptor signalling pathways were enriched within tumour–immune interactions.
- Pathways including MIF, GALECTIN, CCL, CD40, CD99, LCK and ICAM demonstrated differential communication strength.
- Subcluster analysis revealed heterogeneity in signalling dynamics within immune populations.
- NicheNet identified potential upstream ligands influencing CES1-associated transcriptional activity.

---

## How to Run

1. Install required R packages.
2. Obtain the CRC scRNA-seq dataset.
3. Load the dataset into Seurat.
4. Run scripts sequentially from the scripts/ directory.

---

## Required R Packages

- Seurat
- CellChat
- NicheNet
- tidyverse
- ggplot2
- dplyr
- patchwork

---

## Skills Demonstrated

- Single-cell RNA-seq analysis
- Network inference modelling
- Pathway-level signalling interpretation
- Reproducible R workflows
- Biological data visualisation
- Systems-level tumour microenvironment analysis

---

## Author

Jacintha Beena Mathias
MSc Data Science Dissertation  
University of Bristol
