#Script: Load and Explore CRC scRNA-seq Dataset (Seurat setup)

#1.Loading Libraries
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(RColorBrewer)

#2.Load Metadata and TPM Expression Matrix

#Set file paths
dataDir <- "E:/Desktop/Dissertation/Data/Colon_zhang"
outDir <- "E:/Desktop/Dissertation/Data/Colon_zhang/output"

#Load metadata and TPM
meta.10x <- read.csv(paste(dataDir,"/CRC.Leukocyte.10x.Metadata.txt",sep=""), header=T, as.is=T,sep="\t") #   
tpm_mat.10x <- read.table(paste(dataDir,"/CRC.Leukocyte.10x.TPM.txt",sep=""), header=T, as.is=T,sep=" ") #  33695 genes x 63689 cells

#3. Preprocess TPM Matrix
#this step stores non zero entries, space efficient format
mat_sparse <- Matrix(as.matrix(tpm_mat.10x), sparse = TRUE)
log_tpm_mat_sparse <- mat_sparse
#TPM matrix are highly skewed, with some very large values and many small ones therefore log2 is used.
#normally distributed, reduces extreme outliers, makes expression patterms easier to interpret visually.
log_tpm_mat_sparse@x <- log2(log_tpm_mat_sparse@x +1 ) #log2(TPM + 1)

#4.Create Seurat Object 
seurat <- CreateSeuratObject(
  counts = log_tpm_mat_sparse,
  min.cells = 0,
  min.features = 0,
  assay = "scRNA",
)

#Assign TPM values directly to normalized data slot
seurat@assays[["scRNA"]]$data <- seurat@assays[["scRNA"]]$counts #Assign log(TPM) to dataslot

#5.Add Metadata Columns
#here the metadata rows are attached to the seurat objext rows
meta_cols <- c("Sample", "Tissue", "raw.nUMI", "filter.nUMI", "filter.nGene",
               "ribo.per", "Global_Cluster", "Sub_Cluster", "Sub_ClusterID",
               "Global_tSNE_1", "Global_tSNE_2", "Sub_tSNE_1", "Sub_tSNE_2",
               "integrated.Global_tSNE_1", "integrated.Global_tSNE_2")

for (col in meta_cols){
  seurat <- AddMetaData(seurat, metadata = meta.10x[[col]], col.name = col)
}

#Rename tSNE fields to standard Seurat reduction names
names(seurat@meta.data) <- gsub("Sub_tSNE_", "sTSNE_", names(seurat@meta.data))
names(seurat@meta.data) <- gsub("Global_tSNE_", "gTSNE_", names(seurat@meta.data))

#6.Create Dimensionality Reductions
seurat@reductions[["sTSNE"]] <- CreateDimReducObject(
  embeddings = as.matrix(seurat@meta.data[, c("sTSNE_1", "sTSNE_2")]),
  key = "sTSNE_",
  assay = "scRNA"
)
seurat@reductions[["gTSNE"]] <- CreateDimReducObject(
  embeddings = as.matrix(seurat@meta.data[, c("gTSNE_1", "gTSNE_2")]),
  key = "gTSNE_",
  assay = "scRNA"
)

# 7. Quality Control Filtering
seurat <- subset(seurat, subset = filter.nGene > 300 & filter.nUMI > 500 & ribo.per < 50)

# 8. Identify Highly Variable Genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# 9. Scale Data (Z-score normalization)
seurat <- ScaleData(seurat, features = VariableFeatures(seurat))

# 10. Dimensionality Reduction (PCA & UMAP)
seurat <- RunPCA(seurat, features = VariableFeatures(seurat), verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:20)

#11. Summary Statistics

cat("Number of cells:", ncol(seurat), "\n")
cat("Number of genes:", nrow(seurat), "\n")
cat("Unique patients/samples:", length(unique(seurat$Sample)), "\n\n")

cat("Cell types (Global_Cluster):\n")
print(table(seurat$Global_Cluster))

cat("\nTissue types:\n")
print(table(seurat$Tissue))

#12.Save Processed Object
save(seurat, file = file.path(outDir, "seurat_10x_logTPM.RData"))

#Load your saved seurat object
load("E:/Desktop/Dissertation/Data/Colon_zhang/output/seurat_10x_logTPM.RData")

#Set Output Path
output_dir <- "E:/Desktop/Dissertation/Data/Colon_zhang/output/Understanding_Dataset"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


#VISUALIZATIONS
#1.Cell type composition barplot.
#how many cells per Global_Cluster?
library(ggplot2)
png(file.path(output_dir, "Cells per Global_Cluster.png"), width =1000, height=1000, res = 300)
ggplot(seurat@meta.data, aes(x = Global_Cluster)) +
  geom_bar(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Cell Type Composition", x = "Cell Type", y = "Number of Cells")
dev.off()

#2.Tissue Type Composition Barplot
#How many cells are from tumour (T) vs normal (N)?
png(file.path(output_dir, "Tissue Type Composition.png"), width =1000, height=1000, res = 300)
ggplot(seurat@meta.data, aes(x=Tissue))+
  geom_bar(fill = "red") +
  theme_minimal() +
  labs(title = "Tissue Type Composition", x = "Tissue", y = "Number of Cells")
dev.off()

#3.t-SNE Plot by Cell Type
#Are cell type well separated in expression space?
png(file.path(output_dir, "t-SNE - Cells by Global Cluster.png"), width =1600, height=1000, res = 300)
DimPlot(seurat, reduction = "gTSNE", group.by = "Global_Cluster", label = TRUE) + 
  ggtitle("t-SNE - Cells by Global Cluster")
#print(a2)
dev.off()

#4.t-SNE Plot by Tissue
#Do tumour and normal tissues form distinct clusters?
png(file.path(output_dir, "t-SNE - Cells by Tissue Type.png"), width =1600, height=1000, res = 300)
DimPlot(seurat, reduction = "gTSNE", group.by = "Tissue", cols = c("red", "blue", "green")) + 
  ggtitle("t-SNE - Cells by Tissue Type")
dev.off()

#5 PCA plot by cell types
# Run PCA on your Seurat object
# Step 1: Identify the most variable genes
seurat <- FindVariableFeatures(seurat)
# Step 2: Scale the data (recommended to only scale variable features)
seurat <- ScaleData(seurat, features = VariableFeatures(seurat))
# Step 3: Run PCA using only variable features
seurat <- RunPCA(seurat, features = VariableFeatures(seurat), verbose = FALSE)

#Visualize PCA with cell types(Global_Cluster)
png(file.path(output_dir, "PCA with Cell Types - Global Cluster.png"), width =1800, height=1400, res = 300)
DimPlot(seurat, reduction = "pca", group.by = "Global_Cluster", label = TRUE) +
  ggtitle("PCA with Cell Types - Global Cluster")
dev.off()
#Visualize PCA with tissue types
png(file.path(output_dir, "PCA - Cells Colored by Tissue Type.png"), width =1800, height=1400, res = 300)
DimPlot(seurat, reduction = "pca", group.by = "Tissue") + 
  ggtitle("PCA - Cells Colored by Tissue Type")
dev.off()
#PCA Elbow Plot (to choose number of Principle components)
png(file.path(output_dir, "PCA - Elbow Plot .png"), width =1600, height=1000, res = 300)
ElbowPlot(seurat)
dev.off()

#6 UMAP Plot
#Run UMAP on Seurat Object
seurat <- RunUMAP(seurat, dims = 1:20)
#Plot UMAP Colored by Global_Cluster(Cell Type)
png(file.path(output_dir, "UMAP - Cells Colored by Global Cluster.png"), width =1800, height=1400, res = 300)
DimPlot(seurat, reduction = "umap", group.by = "Global_Cluster", label = TRUE) +
  ggtitle("UMAP - Cells Colored by Global Cluster")
dev.off()

#Plot UMAP colored by Tissue Type(N,P,T)
png(file.path(output_dir, "UMAP - Cells Colored by Tissue Type.png"), width =1800, height=1400, res = 300)
DimPlot(seurat, reduction = "umap", group.by = "Tissue") +
  ggtitle("UMAP - Cells Colored by Tissue Type")
dev.off()

#UMAP split by Tissue type and grouped by Global_Cluster
png(file.path(output_dir, "UMAP split by tissue type and grouped by Global_Cluster.png"), width =2000, height=1600, res = 300)
DimPlot(seurat, reduction = "umap", split.by = "Tissue", group.by = "Global_Cluster", ncol =3)
dev.off()

#7.Violin Plots of QC Metrics
png(file.path(output_dir, "QC Metrics per cell.png"), width = 2600, height=2000, res = 300)
VlnPlot(seurat, features = c("filter.nUMI", "filter.nGene", "ribo.per"), pt.size = 0.1, ncol = 3) + 
  ggtitle("QC Metrics per Cell")
dev.off()

#6.Cell Type by Tissue
#Do certain cell types occur more in Tumour or normal?
celltype_tissue_counts <- table(seurat$Tissue, seurat$Global_Cluster)
celltype_tissue_counts
#convert to dataframe for ggplot2.
df <- as.data.frame(celltype_tissue_counts)
colnames(df) <- c("Tissue", "CellType", "Count")
#Plot heatmap using ggplot
png(file.path(output_dir, "Cell Type Distribution by Tissue.png"), width =1800, height=1400, res = 300)
ggplot(df, aes(x = CellType, y = Tissue, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd"), name = "Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Distribution by Tissue", x = "Cell Type", y = "Tissue")
dev.off()

#7.Customized interaction plots for selected signaling pathways
# Define pathways of interest
pathways.of.interest <- c("MIF", "GALECTIN", "CD99", "LCK", "SPP1")

# Loop through and generate circle plots for selected pathways
for (pathway in pathways.of.interest) {
  file_path <- file.path(output_dir, paste0("CirclePlot_", pathway, ".png"))
  
  png(file_path, width = 16, height = 14, units = "in", res = 300)
  
  # Title is automatically generated based on pathway
  netVisual_aggregate(
    object = cellchat,
    signaling = pathway,
    layout = "circle",
    vertex.label.cex = 0.8,
    #vertex.label.dist = 2,
    edge.width.max = 6,
    signaling.name = pathway          # avoid title.name conflict
  )
  
  dev.off()
}

