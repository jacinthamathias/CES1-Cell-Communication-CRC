# Large-Scale Data Analysis & Modelling (43k+ Profiles)

## Overview

This project analyses over 43,000 single-cell transcriptomic profiles to investigate interaction patterns between tumour and immune cell populations in colorectal cancer.

The focus of the project was structured data processing, feature engineering, network inference modelling, and reproducible analytical workflows using R. The objective was to quantify interaction strength, identify key signalling drivers, and generate interpretable visual outputs from high-dimensional data.

---

## Objectives

- Process and structure large, high-dimensional datasets
- Perform data cleaning, validation and normalisation
- Engineer features for interaction modelling
- Quantify communication strength between defined cell groups
- Produce clear visual summaries to support interpretation

---

## Methodology

### 1. Data Processing & Structuring
- Quality control filtering and validation
- Normalisation and variable feature selection
- Dimensionality reduction (PCA, UMAP)
- Structured cell-type annotation using Seurat

### 2. Interaction Network Modelling
- Modelled global communication networks using CellChat
- Quantified pathway-level signalling strength
- Compared interaction patterns across tumour and immune subsets
- Visualised results using circle plots, bubble plots and chord diagrams

### 3. Ligand–Target Modelling
- Applied NicheNet to identify candidate upstream ligands
- Ranked ligands influencing downstream transcriptional activity
- Performed subgroup-level comparison to explore signalling heterogeneity

---

## Key Findings

- CES1 expression varies across tumour-associated macrophages and immune subsets.
- Multiple ligand–receptor pathways showed differential interaction strength.
- Pathways including MIF, GALECTIN, CCL, CD40, CD99, LCK and ICAM demonstrated enriched tumour–immune signalling.
- Subgroup analysis revealed heterogeneity in signalling dynamics within immune populations.
- Candidate upstream ligands influencing CES1-associated transcriptional activity were identified.

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

## How to Run

1. Install required R packages.
2. Obtain the CRC scRNA-seq dataset.
3. Load the dataset into Seurat.
4. Run scripts sequentially from the `scripts/` directory.

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

- Large-scale data cleaning and transformation  
- Feature engineering  
- Network inference modelling  
- Dimensionality reduction  
- Reproducible R workflows  
- Data visualisation and structured result interpretation  

---

Author:  
Jacintha Beena Mathias  
MSc Data Science, University of Bristol
