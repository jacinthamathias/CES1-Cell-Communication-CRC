# Dataset Description

This project utilises colorectal cancer (CRC) single-cell RNA sequencing (scRNA-seq) data generated using 10x Genomics technology.

## Data Components

- Gene expression count matrix
- Cell-level metadata
- Cell type annotations
- Tumour microenvironment cell populations

Key gene of interest:
- CES1 (Carboxylesterase 1)

## Data Access

The raw dataset is not included in this repository due to:

- File size constraints
- Licensing restrictions
- Ethical data-sharing considerations

To reproduce the analysis:

1. Obtain the CRC scRNA-seq dataset from the original publication.
2. Load expression matrix and metadata into Seurat.
3. Execute scripts sequentially from the scripts/ directory.

## Preprocessing Summary

The following steps were performed:

- Quality control filtering
- Normalisation
- Variable feature selection
- PCA and UMAP dimensionality reduction
- Cell clustering and annotation
