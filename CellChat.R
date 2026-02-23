#Install and load required packages
library(CellChat)
library(patchwork)
library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)

#Set Output Path
output_dir <- "E:/Desktop/Dissertation/Data/Colon_zhang/output/CellChat"
dir.create(output_dir, showWarnings = FALSE)

#Load the Subset & Prepare CellChat Object
#Load your small Seurat object
load("E:/Desktop/Dissertation/Data/Colon_zhang/output/seurat_10x_logTPM.RData")

#Set cell identities based on Global_Cluster
Idents(seurat) <- seurat$Global_Cluster

#Create CellChat Object
cellchat <- createCellChat(object = seurat, group.by = "Global_Cluster", assay = "scRNA")

#Assign CellChatDB human ligand-receptor interaction database
cellchat@DB <- CellChatDB.human

#Preprocess CellChat Data
#Subset Data to Ligand-Receptor Genes & Preprocess
#Only keep genes involved in signalling
cellchat <- subsetData(cellchat)

#Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)

#Identify overexpressed ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)

#Project data onto the protein-protein interaction network
cellchat <- projectData(cellchat, PPI.human)

#Compute Cell-Cell Communication Probabilities
#Compute communication probability
cellchat <- computeCommunProb(cellchat)

#Ligand-Receptor table
df.net <- subsetCommunication(cellchat)
head(df.net)

#Filter out low-confidence interactions
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Infer communication at pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Aggregate network info
cellchat <- aggregateNet(cellchat)

#Get group size
groupSize <- as.numeric(table(cellchat@idents))

#Visualize the Communication Network

#1. Circle Plots (Interaction Count & Strength)

png(file.path(output_dir, "Interaction_count_circle.png"), width = 3000, height=3000, res = 300)
# Expand plot margins and allow drawing outside the plot area
par(mar = c(1, 1, 5, 1))     # bottom, left, top, right
par(xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of Interactions")
dev.off()

#2. Circle Plots (Interaction Strength)

png(file.path(output_dir, "Interaction_Strength.png"), width =3000, height = 3000, res = 300)
# Expand plot margins and allow drawing outside the plot area
par(mar = c(1, 1, 5, 1))     # bottom, left, top, right
par(xpd = TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction Strength")
dev.off()

#3. Number of interactions
# Expand plot margins and allow drawing outside the plot area
png(file.path(output_dir, "Heatmap_count.png"), width = 3000, height = 2000, res = 300)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

#4. Strength of interactions
png(file.path(output_dir, "Heatmap_strength.png"), width = 3000, height = 2000, res = 300)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

#5. Bubble Plot - Signaling Pathway between clusters
png(file.path(output_dir, "Bubbleplot_all_pathways.png"), width=3000, height=2000, res = 300)
netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
dev.off()

#6. Chord Diagram for a specific pathway0
png(file.path(output_dir, "Chord_MIF.png"), width = 4000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = "MIF", layout = "chord")
dev.off()

png(file.path(output_dir, "Chord_GALECTIN.png"), width = 8, height = 8, units = "in", res = 300)
netVisual_aggregate(cellchat, signaling = "GALECTIN", layout = "chord")
dev.off()

#png(file.path(output_dir, "Chord_CCL.png"), width = 4000, height = 3000, res = 300)
#netVisual_aggregate(cellchat, signaling = "CCL", layout = "chord")
#dev.off()

png(file.path(output_dir, "Chord_CD99.png"), width = 4000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = "CD99", layout = "chord")
dev.off()

#png(file.path(output_dir, "Chord_ICAM.png"), width = 4000, height = 3000, res = 300)
#netVisual_aggregate(cellchat, signaling = "ICAM", layout = "chord")
#dev.off()

png(file.path(output_dir, "Chord_LCK.png"), width = 4000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = "LCK", layout = "chord")
dev.off()

#png(file.path(output_dir, "Chord_CD40.png"), width = 4000, height = 3000, res = 300)
#netVisual_aggregate(cellchat, signaling = "CD40", layout = "chord")
#dev.off()

#7. Annotated Chord diagram based on CES1 expression

#Get average CES1 expression per cluster(per celltype)
ces1_avg <- AverageExpression(seurat, features = "CES1", group.by = "Global_Cluster")$scRNA
print(ces1_avg)

#2.Set CES1-Based Node Colors
#Normalize CES1 values between 0-1
ces1_norm <- (ces1_avg - min(ces1_avg)) / (max(ces1_avg) - min(ces1_avg))

#Map to color gradient
library(RColorBrewer)
library(scales)
ces1_colors <- scales::col_numeric(palette = c("#4575b4", "#ffdd57", "#d73027" ), domain = NULL)(as.vector(ces1_norm))

# Create named color vector for nodes
names(ces1_colors) <- names(ces1_avg)

#Plot annotated chord diagrams with CES1 node

png(file.path(output_dir, "MIF_Annotated.png"), width = 4000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = "MIF", layout = "chord", color.use = ces1_colors)
dev.off()

png(file.path(output_dir, "GALECTIN_Annotated.png"), width = 4000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = "GALECTIN", layout = "chord", color.use = ces1_colors)
dev.off()

#png(file.path(output_dir, "CCL_Annotated.png"), width = 4000, height = 3000, res = 300)
#netVisual_aggregate(cellchat, signaling = "CCL", layout = "chord", color.use = ces1_colors)
#dev.off()

png(file.path(output_dir, "CD99_Annotated.png"), width = 4000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = "CD99", layout = "chord", color.use = ces1_colors)
dev.off()

#png(file.path(output_dir, "ICAM_Annotated.png"), width = 4000, height = 3000, res = 300)
#netVisual_aggregate(cellchat, signaling = "ICAM", layout = "chord", color.use = ces1_colors)
#dev.off()

png(file.path(output_dir, "LCK_Annotated.png"), width = 4000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = "LCK", layout = "chord", color.use = ces1_colors)
dev.off()

#png(file.path(output_dir, "CD40_Annotated.png"), width = 4000, height = 3000, res = 300)
#netVisual_aggregate(cellchat, signaling = "CD40", layout = "chord", color.use = ces1_colors)
#dev.off()

#Overlay CES1 Expression Across Cell Types
#Check if CES1 exists
"CES1" %in% rownames(seurat)

#8. Plot CES1 expression by cluster
#Set identities to match group.by
#Create expression bins
#Prepare expression and coordinates
ces1_data <- FetchData(seurat, vars = c("CES1", "gTSNE_1", "gTSNE_2"))

#Create expression bin categories
ces1_data$bin <- cut(
  ces1_data$CES1,
  breaks = c(0, 0.5, 1, 1.5, 2, Inf),
  labels = c("0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2")
)

#Color palette with NA (grey)
ces1_colors <- c(
  "0–0.5" = "#20908d",
  "0.5–1" = "#feb24c",
  "1–1.5" = "#3a528b",
  "1.5–2" = "#bd0026",
  ">2" = "#005824",
  "NA" = "#d9d9d9"
)

#Set factor levels
ces1_data$bin <- as.character(ces1_data$bin)
ces1_data$bin[is.na(ces1_data$bin)] <- "NA"
ces1_data$bin <- factor(ces1_data$bin, levels = names(ces1_colors))

# Plot with two layers
library(ggplot2)

p <- ggplot() +
  # Plot all cells in grey (background)
  geom_point(data = ces1_data, aes(x = gTSNE_1, y = gTSNE_2), color = "#d9d9d9", size = 0.4) +
  # Overlay cells with CES1 > 0 (excluding bin == "NA")
  geom_point(data = ces1_data[ces1_data$bin != "NA", ], aes(x = gTSNE_1, y = gTSNE_2, color = bin), size = 0.7) +
  scale_color_manual(values = ces1_colors) +
  ggtitle("t-SNE - CES1 Expression ") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

#Save plot
ggsave(
  filename = file.path(output_dir, "CES1_Expression_Across_Cell_Types.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

#9. PCA with CES1 Expression
# Run PCA on your Seurat object
#Identify the most variable genes
seurat <- FindVariableFeatures(seurat)

#Scale the data (recommended to only scale variable features)
seurat <- ScaleData(seurat, features = VariableFeatures(seurat))

#Run PCA using only variable features
seurat <- RunPCA(seurat, features = VariableFeatures(seurat))

#Fetch CES1 expression and PCA coordinates
ces1_pca <- FetchData(seurat, vars = c("CES1", "PC_1", "PC_2"))

#Bin CES1 expression values
ces1_pca$bin <- cut(
  ces1_pca$CES1,
  breaks = c(0, 0.5, 1, 1.5, 2, Inf),
  labels = c("0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2")
)

#Assign NA label to zero-expressing cells
ces1_pca$bin <- as.character(ces1_pca$bin)
ces1_pca$bin[is.na(ces1_pca$bin)] <- "NA"
ces1_pca$bin <- factor(ces1_pca$bin, levels = c("0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2", "NA"))

#Define color palette
ces1_colors <- c(
  "0–0.5" = "#4292c6",
  "0.5–1" = "#feb24c",
  "1–1.5" = "#f03b20",
  "1.5–2" = "#bd0026",
  ">2" = "#005824",
  "NA" = "#d9d9d9"
)

#Plot with grey background and color overlay
library(ggplot2)

f_pca <- ggplot() +
  # Background layer: all cells as grey
  geom_point(data = ces1_pca, aes(x = PC_1, y = PC_2), color = "#d9d9d9", size = 0.5) +
  # Overlay CES1-expressing cells
  geom_point(data = ces1_pca[ces1_pca$bin != "NA", ], aes(x = PC_1, y = PC_2, color = bin), size = 0.7) +
  scale_color_manual(values = ces1_colors) +
  ggtitle("PCA - CES1 Expression") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

#Save the plot
ggsave(filename = file.path(output_dir, "PCA_CES1_Expression.png"),
       plot = f_pca,
       width = 8,
       height = 6,
       dpi = 300,
       bg = "white")

#10 UMAP with CES1 Expression
seurat <- RunUMAP(seurat, dims = 1:20)

library(ggplot2)

# Step 1: Calculate CES1 expression bins
ces1_expr <- FetchData(seurat, vars = "CES1")
ces1_expr$CES1_bin <- cut(
  ces1_expr$CES1,
  breaks = c(0, 0.5, 1, 1.5, 2, Inf),
  labels = c("0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2")
)
ces1_expr$CES1_bin <- as.character(ces1_expr$CES1_bin)
ces1_expr$CES1_bin[is.na(ces1_expr$CES1_bin)] <- "NA"
ces1_expr$CES1_bin <- factor(ces1_expr$CES1_bin, levels = c("0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2", "NA"))

# Step 2: Add the CES1_bin column to Seurat metadata
seurat$CES1_bin <- ces1_expr$CES1_bin


# Extract UMAP coordinates
umap_coords <- Embeddings(seurat, "umap") %>% as.data.frame()
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")  # Ensure correct names

# Add CES1 bin info
umap_coords$CES1_bin <- seurat$CES1_bin

# Set plotting order (NA first so color points go on top)
umap_coords$CES1_bin <- factor(
  umap_coords$CES1_bin,
  levels = c(NA, "0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2")
)

# Create plot with grey NA background
f_umap_layered <- ggplot() +
  geom_point(data = umap_coords[is.na(umap_coords$CES1_bin), ],
             aes(x = UMAP_1, y = UMAP_2),
             color = "#d9d9d9", size = 1, alpha = 0.5) +
  geom_point(data = umap_coords[!is.na(umap_coords$CES1_bin), ],
             aes(x = UMAP_1, y = UMAP_2, color = CES1_bin),
             size = 1) +
  scale_color_manual(
    values = c(
      "0–0.5" = "#4292c6",    # deeper blue
      "0.5–1" = "#feb24c",
      "1–1.5" = "#f03b20",
      "1.5–2" = "#bd0026",
      ">2" = "#005824"
    ),
    name = "CES1 Level"
  ) +
  ggtitle("UMAP - CES1 Expression") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Save to file
ggsave(
  filename = file.path(output_dir, "UMAP - CES1 Expression.png"),
  plot = f_umap_layered,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


#11. Boxplot CES1 expression in Tumor vs Normal tissue
# Fetch expression and metadata
ces1_expr <- FetchData(seurat, vars = "CES1")
ces1_expr$Tissue <- seurat$Tissue

# Generate boxplot
b1 <- ggplot(ces1_expr, aes(x = Tissue, y = CES1, fill = Tissue)) +
  geom_boxplot(outlier.colour =  "black", outlier.shape = 16, outlier.size = 1.5) +
  theme_minimal() +
  labs(title = "CES1 Expression - Tumour vs Normal", y = "Log(TPM+1)", x = "") +
  scale_fill_manual(values = c("T" = "firebrick", "N" = "steelblue"))

b1 <- ggplot(ces1_expr, aes(x = Tissue, y = CES1, fill = Tissue)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # avoid double dots
  geom_jitter(shape = 16, size = 1, width = 0.1, aes(color = Tissue)) +
  theme_minimal() +
  labs(title = "CES1 Expression: Tumour vs Normal", y = "Log(TPM+1)", x = "") +
  scale_fill_manual(values = c("T" = "firebrick", "N" = "steelblue")) +
  scale_color_manual(values = c("T" = "firebrick", "N" = "steelblue"))

# Save the plot
ggsave(
  filename = file.path(output_dir, "CES1 Expression - Tumour vs Normal.png"),
  plot = b1,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

#Link CES1 Expression to Communication Network
#Get full interaction dataframe
comm_df <- subsetCommunication(cellchat)

#Average CES1 expression by cell type
ces1_avg <- AverageExpression(seurat, features = "CES1", group.by = "Global_Cluster")$scRNA

#Add CES1 expression in sender and receiver cells
# Extract numeric vector for each source (column name)
comm_df$CES1_sender <- ces1_avg[1, as.character(comm_df$source)]

comm_df$CES1_receiver <- ces1_avg[1, as.character(comm_df$target)]

#12.Visualize CES1 in senders and receiver cells
p1 <- ggplot(comm_df, aes(x = CES1_sender, y = prob)) + 
  geom_point(aes(color = source)) +
  theme_minimal() +
  labs(title = "CES1 in sender Cells vs Interaction Strength", x = "Avg CES1 in Sender", y = "Interaction Probability")
ggsave(filename = file.path(output_dir, "CES1 in sender Cells vs Interaction Strength.png"), plot = p1, width = 6, height = 4, dpi = 300, bg = "white")

p2 <- ggplot(comm_df, aes(x = CES1_receiver, y = prob)) + 
  geom_point(aes(color = target)) +
  theme_minimal() +
  labs(title = "CES1 in Receiver Cells vs Interaction Strength", x = "Avg CES1 in Receiver", y = "Interaction Probability")
ggsave(filename = file.path(output_dir, "CES1 in Receiver Cells vs Interaction Strength.png"), plot = p2, width = 6, height = 4, dpi = 300, bg = "white")

#Save Processes CellChat Object
save(cellchat, file = "cellchat_CRC_CES1_analysis.RData")

#13. Count of each celltype as source and target
#Extract top 100 interactions(based on communication probability)
library(dplyr)
library(ggplot2)
library(patchwork)

# Get all interactions
all_communications <- subsetCommunication(cellchat)

# Top 100 by communication probability
top_100_pathways <- all_communications %>% 
  arrange(desc(prob)) %>%
  slice(1:100)

#Count how often each cell type appears as source and target
# Count source frequency
source_counts <- as.data.frame(table(top_100_pathways$source))
colnames(source_counts) <- c("cell_type", "count")

# Count target frequency
target_counts <- as.data.frame(table(top_100_pathways$target))
colnames(target_counts) <- c("cell_type", "count")

#create bar plots
# Source plot
p1 <- ggplot(source_counts, aes(x = reorder(cell_type, -count), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(title = "Count of Each Cell Type as Source",
       x = "Source Cell Type", y = "Count")
ggsave(filename = file.path(output_dir, "Count of Each Cell Type as Source.png"), plot = p1, width = 6, height = 4, dpi = 300, bg = "white")

# Target plot
p2 <- ggplot(target_counts, aes(x = reorder(cell_type, -count), y = count)) +
  geom_bar(stat = "identity", fill = "salmon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(title = "Count of Each Cell Type as Target",
       x = "Target Cell Type", y = "Count")
ggsave(filename = file.path(output_dir, "Count of Each Cell Type as Target.png"), plot = p2, width = 6, height = 4, dpi = 300, bg = "white")


