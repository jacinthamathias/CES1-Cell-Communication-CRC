#Niechenet on celltype in detail like Myeloid etc

# Load required packages
library(Seurat)
library(nichenetr)
library(tidyverse)
library(future)
library(Matrix)
library(ComplexHeatmap)
library(gridExtra)
library(ggplot2)

# Setup parallel processing for speed
plan(sequential)
options(future.globals.maxSize = 10 * 1024^3)

#Load Data
load("C:/Users/Hp/Desktop/Dissertation/Data/Colon_zhang/output/seurat_10x_logTPM.RData")
output_dir <- "C:/Users/Hp/Desktop/Dissertation/Data/Colon_zhang/output/NiecheNet_eachCellType"
dir.create(output_dir, showWarnings = FALSE)

# Set correct assay
DefaultAssay(seurat) <- "scRNA"
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
seurat <- ScaleData(seurat, features = VariableFeatures(seurat), verbose = FALSE)
seurat <- RunPCA(seurat, features = VariableFeatures(seurat), verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30)

# Define cell types using broader matching
sender_patterns <- c("TAM", "Macro", "Mono")
receiver_patterns <- c("CD8", "B", "CD4", "ILC")
celltype_col <- "Sub_Cluster"

#sender_patterns <- "Tumor"
#receiver_patterns <- c("CD8 T cell", "B cell", "Myeloid cell", "CD4 T cell", "ILC")

# Set cell identities
Idents(seurat) <- factor(seurat@meta.data[[celltype_col]])
all_clusters <- levels(Idents(seurat))
sender_celltype <- grep(paste(sender_patterns, collapse="|"), all_clusters, value = TRUE)
receiver_celltypes <- grep(paste(receiver_patterns, collapse="|"), all_clusters, value = TRUE)

if (length(sender_celltype) == 0) stop("No sender clusters matched the pattern")
if (length(receiver_celltypes) == 0) stop("No receiver clusters matched the pattern")

# Safe wrapper for get_expressed_genes
safe_get_expr_genes <- function(seurat_obj, cluster) {
  tryCatch({
    cells <- WhichCells(seurat_obj, idents = cluster)
    cat("Cluster:", cluster, "| Cells:", length(cells), "\n")
    if (length(cells) == 0) return(character())
    expr_data <- GetAssayData(seurat_obj, slot = "data")
    expr_data <- expr_data[, cells, drop = FALSE]
    gene_pct <- Matrix::rowMeans(expr_data > 0)
    names(gene_pct)[gene_pct >= 0.10]
  }, error = function(e) {
    message(paste("Failed to compute expressed genes for:", cluster, "Error:", e$message))
    return(character())
  })
}

# Get expressed genes
expressed_genes_sender_list <- map(sender_celltype, ~safe_get_expr_genes(seurat, .x)) %>% discard(~length(.x) == 0)
expressed_genes_sender <- if (length(expressed_genes_sender_list) > 0) Reduce(union, expressed_genes_sender_list) else character()

expressed_genes_receiver_list <- map(receiver_celltypes, ~safe_get_expr_genes(seurat, .x)) %>% discard(~length(.x) == 0)
expressed_genes_receiver <- if (length(expressed_genes_receiver_list) > 0) Reduce(intersect, expressed_genes_receiver_list) else character()

if (length(expressed_genes_receiver) == 0 || length(expressed_genes_sender) == 0) stop("No expressed genes found in receiver or sender clusters")

# DEGs in first receiver cluster (as example)
example_receiver <- receiver_celltypes[1]
deg_receiver <- FindMarkers(seurat, ident.1 = example_receiver, group.by = celltype_col, 
                            only.pos = TRUE, logfc.threshold = 0.25) %>%
  rownames_to_column("gene") %>%
  top_n(200, wt = avg_log2FC) %>%
  pull(gene)

# Load NicheNet data
options(timeout = 300)
ligand_target_matrix <- readRDS(url("https://zenodo.org/records/5884439/files/ligand_target_matrix_nsga2r_final.rds"))
ligand_receptor_network <- readRDS(url("https://zenodo.org/records/5884439/files/lr_network_human_21122021.rds"))

ligands <- unique(ligand_receptor_network$from)
ligands_sender <- intersect(ligands, expressed_genes_sender)

# Predict ligand activities
ligand_activities <- predict_ligand_activities(
  geneset = deg_receiver,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligands_sender
)

#Heatmap for Niechenet
top_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(desc(pearson)) %>% pull(test_ligand)
valid_top_ligands <- intersect(top_ligands, rownames(ligand_target_matrix))
valid_deg_receiver <- intersect(deg_receiver, colnames(ligand_target_matrix))

if (length(valid_top_ligands) == 0 || length(valid_deg_receiver) == 0) {
  stop("No overlap between top ligands or DEGs and ligand_target_matrix")
}

heatmap_targets <- ligand_target_matrix[valid_top_ligands, valid_deg_receiver, drop = FALSE]

# Save heatmap
png(file.path(output_dir, "nichenet_heatmap.png"), width = 1200, height = 1000)
if (nrow(heatmap_targets) > 0 && ncol(heatmap_targets) > 0) {
  draw(Heatmap(heatmap_targets, name = "Regulation"))
} else {
  message("No data for heatmap after filtering")
}
dev.off()

#Cell-type Specific Ligand Expression DotPlot
#Ligand Expression DotPlot
top_ligands_plot <- top_ligands[1:10]

png(file.path(output_dir, "LigandExpression_DotPlot.png"), width = 2000, height = 1000, res = 150)
DotPlot(seurat, features = top_ligands_plot, group.by = celltype_col) +
  RotatedAxis() +
  ggtitle("Ligand Expression Across Cell Types")
dev.off()

#Receptor Expression DotPlot in Receiver Cells
#Receptor Expression DotPlot
receptors <- lr_network$to %>% unique() %>% intersect(rownames(seurat))
receptor_subset <- head(receptors, 20)

png(file.path(output_dir, "ReceptorExpression_DotPlot.png"), width = 2000, height = 1000, res = 150)
DotPlot(seurat, features = receptor_subset, group.by = celltype_col) +
  RotatedAxis() +
  ggtitle("Receptor Expression in Receiver Clusters")
dev.off()

#Violin or Ridge Plots for CES1 Expression
library(ggplot2)
png(file.path(output_dir, "CES1_ViolinPlot.png"), width = 2000, height = 1000, res = 150)
VlnPlot(seurat, features = "CES1", group.by = celltype_col, pt.size = 0.1) +
  ggtitle("CES1 Expression Across Cell Types") +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 2), "cm")  # top, right, bottom, left
  )
dev.off()

# Correlation with CES1
ces1_expr <- FetchData(seurat, vars = "CES1")$CES1
expr_matrix <- GetAssayData(seurat, slot = "data")
expr_matrix <- expr_matrix[Matrix::rowMeans(expr_matrix > 0) >= 0.1, ]
expr_scaled <- t(scale(t(expr_matrix), center = TRUE, scale = TRUE))
ces1_expr <- scale(ces1_expr, center = TRUE, scale = TRUE)
cor_scores <- Matrix::rowMeans(expr_scaled * rep(ces1_expr, each = nrow(expr_scaled)))

# Plot top correlated genes
top_pos <- head(sort(cor_scores, decreasing = TRUE), 20)
top_neg <- head(sort(cor_scores, decreasing = FALSE), 20)

# Set consistent theme
base_theme <- theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"))

# Positive correlation plot
g1 <- ggplot(tibble(Gene = names(top_pos), Correlation = top_pos),
             aes(x = reorder(Gene, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  ggtitle("Top Positive Correlated Genes with CES1") +
  base_theme

# Negative correlation plot
g2 <- ggplot(tibble(Gene = names(top_neg), Correlation = top_neg),
             aes(x = reorder(Gene, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  ggtitle("Top Negative Correlated Genes with CES1") +
  base_theme

# Save high-resolution vertically stacked PNG
png(file.path(output_dir, "CES1_Correlation_Bars.png"),
    width = 1600, height = 2000, res = 150)
grid.arrange(g1, g2, ncol = 1)
dev.off()

# Save FeaturePlot
seurat <- RunUMAP(seurat, dims = 1:30)
png(file.path(output_dir, "CES1_FeaturePlot.png"), width = 1000, height = 800)
FeaturePlot(
  seurat,
  features = "CES1",
  reduction = "umap",
  cols = c("grey90", "lightblue", "blue", "darkblue"),
  order = TRUE  # ensures high expression cells plotted on top
)

dev.off()

#DotPlot for Top Ligands in Sender Cells
top10_ligands <- top_ligands[1:10]
png(file.path(output_dir, "TopLigands_DotPlot.png"), width = 2000, height = 1000)
DotPlot(seurat, features = top10_ligands, group.by = celltype_col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#Top ligands by activity
library(ggplot2)

# Prepare top ligands data (if not already)
ligand_plot_df <- ligand_activities %>%
  filter(test_ligand %in% top_ligands) %>%
  mutate(rank = rank(-pearson)) %>%
  arrange(rank)

# Save plot as PNG
png(file.path(output_dir, "TopLigand_Activity_Plot.png"), width = 1600, height = 800, res = 150)

ggplot(ligand_plot_df, aes(x = reorder(test_ligand, -pearson), y = pearson)) +
  geom_point(size = 6, color = "steelblue") +
  theme_minimal() +
  labs(
    x = "Ligand",
    y = "Pearson Correlation",
    title = "Top Ligands by Activity (NicheNet)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Receptor Expression in Receiver Cells
top_ligand_receptors <- lr_network %>% filter(from %in% top_ligands) %>% pull(to) %>% unique()
top_receptors <- head(top_ligand_receptors, 10)

png(file.path(output_dir, "TopReceptors_DotPlot.png"), width = 2000, height = 1000)
DotPlot(seurat, features = top_receptors, group.by = celltype_col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#UMAP with CES1 expression by cluster
png(file.path(output_dir, "CES1_by_SubCluster_UMAP.png"), width = 1000, height = 800)
DimPlot(seurat, group.by = "Sub_Cluster", label = TRUE, reduction = "umap") +
  ggtitle("UMAP colored by Sub_Cluster")
dev.off()

#Sender-Receiver Interaction Heatmap
interaction_matrix <- matrix(0, nrow = length(sender_celltype), ncol = length(receiver_celltypes),
                             dimnames = list(sender_celltype, receiver_celltypes))

for (sender in sender_celltype) {
  for (receiver in receiver_celltypes) {
    expressed_ligands <- safe_get_expr_genes(seurat, sender)
    expressed_receptors <- safe_get_expr_genes(seurat, receiver)
    lr_pairs <- lr_network %>% 
      filter(from %in% expressed_ligands, to %in% expressed_receptors)
    interaction_matrix[sender, receiver] <- nrow(lr_pairs)
  }
}

png(file.path(output_dir, "Sender_Receiver_Interaction_Heatmap.png"), width = 1200, height = 1000)
Heatmap(interaction_matrix, name = "LR Pairs", col = colorRamp2(c(0, max(interaction_matrix)), c("white", "darkred")))
dev.off()

#Ligand-Target Bipartite Graph
ligand_target_edges <- reshape2::melt(as.matrix(heatmap_targets)) %>%
  filter(value > 0.1) %>%  # threshold to declutter
  rename(from = Var1, to = Var2)

if (nrow(ligand_target_edges) > 0) {
  g_bipartite <- graph_from_data_frame(ligand_target_edges, directed = TRUE)
  V(g_bipartite)$type <- ifelse(V(g_bipartite)$name %in% top_ligands, "Ligand", "Target")
  
  png(file.path(output_dir, "Ligand_Target_BipartiteGraph.png"), width = 1600, height = 1200, res = 150)
  ggraph(g_bipartite, layout = "fr") +
    geom_edge_link(aes(alpha = 0.5), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) +
    geom_node_point(aes(color = type), size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
    scale_color_manual(values = c("Ligand" = "dodgerblue", "Target" = "orange")) +
    theme_void() +
    ggtitle("Ligand to Target Network")
  dev.off()
}

#Target Gene Module Score in Receiver Cells
#Target Gene Module Score in Receiver Cells
# AddModuleScore can be used if genes > 5
if (length(valid_deg_receiver) >= 5) {
  seurat <- AddModuleScore(seurat, features = list(valid_deg_receiver), name = "TargetModule")
  
  png(file.path(output_dir, "Receiver_Target_ModuleScore_UMAP.png"), width = 1000, height = 800)
  FeaturePlot(seurat, features = "TargetModule1", reduction = "umap", cols = c("lightgrey", "blue", "darkblue"), order = TRUE) +
    ggtitle("Target Gene Module Score in UMAP")
  dev.off()
}

#Barplot: Ligand Regulates How Many DEGs
ligand_target_counts <- rowSums(heatmap_targets > 0.1)
ligand_target_df <- tibble(Ligand = names(ligand_target_counts),
                           Target_Count = ligand_target_counts)

png(file.path(output_dir, "Ligand_Target_Counts_Barplot.png"), width = 1600, height = 1000, res = 150)
ggplot(ligand_target_df, aes(x = reorder(Ligand, Target_Count), y = Target_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  ggtitle("Number of Regulated Target Genes per Ligand") +
  labs(x = "Ligand", y = "# of Target Genes")
dev.off()

# Export outputs
lr_network <- ligand_receptor_network %>% filter(from %in% top_ligands & to %in% expressed_genes_receiver)
write_csv(ligand_activities, file.path(output_dir, "ligand_activities.csv"))
write_csv(lr_network, file.path(output_dir, "top_ligand_receptors.csv"))

message("NicheNet analysis complete on each celltype. Output saved in output_dir.")

