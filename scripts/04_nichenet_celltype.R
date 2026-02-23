# Niechenet analysis on celltype

# Required libraries
library(Seurat)
library(nichenetr)
library(tidyverse)
library(future)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

plan(sequential)  # safer for large objects
options(future.globals.maxSize = 10 * 1024^3)

# Load Seurat object
load("C:/Users/Hp/Desktop/Dissertation/Data/Colon_zhang/output/seurat_10x_logTPM.RData")

# Set Output Directory
output_dir <- "C:/Users/Hp/Desktop/Dissertation/Data/Colon_zhang/output/NiecheNet_celltype"
dir.create(output_dir, showWarnings = FALSE)

# Set annotations
seurat$celltype <- case_when(
  seurat$Tissue == "T" ~ "Tumor",
  seurat$Global_Cluster == "CD4 T cell" ~ "Treg",
  seurat$Global_Cluster == "Macrophage" ~ "TAM",
  seurat$Global_Cluster == "Fibroblast" ~ "Fibroblast",
  seurat$Global_Cluster == "Myeloid cell" ~ "TAM",
  TRUE ~ "Other"
)
Idents(seurat) <- seurat$celltype

# Load NicheNet ligand-target matrix and network
options(timeout = 300)
ligand_target_matrix <- readRDS(url("https://zenodo.org/records/5884439/files/ligand_target_matrix_nsga2r_final.rds"))
ligand_target_matrix <- t(ligand_target_matrix)  # transpose so genes are cols
stopifnot("CES1" %in% colnames(ligand_target_matrix))  # validation

ligand_receptor_network <- readRDS(url("https://zenodo.org/records/5884439/files/lr_network_human_21122021.rds"))
lr_network <- inner_join(ligand_receptor_network, tibble(from = rownames(t(ligand_target_matrix))), by = "from")

# Set sender and receiver
sender_celltypes <- c("Tumor")
receiver_celltypes <- c("Treg", "TAM", "Fibroblast")

# Get expressed ligands
expressed_sender_genes <- get_expressed_genes(seurat, celltype = sender_celltypes, pct = 0.10)
ligands <- unique(lr_network$from)
expressed_ligands <- union(intersect(ligands, expressed_sender_genes), "CES1")

# Start CES1 ligand activity analysis
ces1_summary_df <- data.frame()
all_ligand_activities <- list()

for (receiver in receiver_celltypes) {
  message("Checking CES1 for receiver: ", receiver)
  
  if (!receiver %in% seurat$celltype) next
  
  receiver_cells <- subset(seurat, subset = celltype == receiver)
  if (ncol(receiver_cells) < 10) next
  
  ces1_expr_all <- FetchData(receiver_cells, vars = "CES1")$CES1
  ces1_positive_idx <- ces1_expr_all > 0.1
  if (sum(ces1_positive_idx) < 10) {
    message("Skipping ", receiver, ": too few CES1-positive cells")
    next
  }
  
  expr_mat <- GetAssayData(receiver_cells, slot = "data")[, ces1_positive_idx]
  ces1_expr <- ces1_expr_all[ces1_positive_idx]
  
  cor_results <- apply(expr_mat, 1, function(g) {
    if (sd(g) == 0 || sd(ces1_expr) == 0) return(NA)
    cor(g, ces1_expr, method = "pearson", use = "complete.obs")
  })
  
  cor_genes <- names(cor_results)[!is.na(cor_results) & abs(cor_results) > 0.1]
  background_genes <- get_expressed_genes(seurat, celltype = receiver, pct = 0.10)
  ces1_targets <- names(which(ligand_target_matrix[, "CES1"] > 0))
  
  overlap <- intersect(cor_genes, background_genes)
  ces1_hits <- intersect(overlap, ces1_targets)
  
  ces1_summary_df <- bind_rows(ces1_summary_df, tibble(
    receiver = receiver,
    ces1_target_hits = length(ces1_hits),
    total_targets = length(ces1_targets),
    expressed_targets = length(intersect(ces1_targets, background_genes)),
    cor_targets = length(overlap)
  ))
  
  # Call predict_ligand_activities
  valid_ligands <- intersect(colnames(ligand_target_matrix), expressed_ligands)
  if (!"CES1" %in% colnames(ligand_target_matrix)) {
    stop("CES1 not found in ligand_target_matrix column names. Ensure matrix is transposed correctly.")
  }
  
  
  ligand_activities <- tryCatch({
    predict_ligand_activities(
      geneset = ces1_hits,
      background_expressed_genes = background_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = valid_ligands
    )
  }, error = function(e) {
    message("Failed to compute ligand activities for ", receiver, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(ligand_activities)) {
    ligand_activities$receiver <- receiver
    all_ligand_activities[[receiver]] <- ligand_activities
    
    # Save to CSV
    write_csv(ligand_activities, file.path(output_dir, paste0("ligand_activities_", receiver, ".csv")))
  }
}

# Combine all ligand activities
ligand_activities_combined <- bind_rows(all_ligand_activities)

# Save combined table
write_csv(ligand_activities_combined, file.path(output_dir, "ligand_activities_combined_all_receivers.csv"))

# Plot overlap count
n1<- ggplot(ces1_summary_df, aes(x = receiver, y = ces1_target_hits)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "CES1 Correlated Targets in Each Receiver",
    x = "Receiver Cell Type",
    y = "# CES1 Correlated Target Genes"
  ) +
  theme_minimal(base_size = 14)
ggsave(file.path(output_dir, "CES1_Correlated target genes.png"), n1, width = 8, height=6, bg = "white")

#Required libraries
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(Seurat)

#Violin: CES1 expression in Treg and Fibroblast
n2 <- 
  VlnPlot(seurat, features = "CES1", group.by = "celltype") +
  theme_minimal() +
  labs(title = "CES1 Expression in Treg and Fibroblast") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "CES1_expression_treg_fibroblast.png"), n2, width = 5, height = 5, bg = "white")


#CES1 GO Enrichment dotplot
# Get ligands (rows) that target CES1
ces1_targets <- rownames(ligand_target_matrix)[ligand_target_matrix[, "CES1"] > 0]

#ces1_targets <- names(which(ligand_target_matrix["CES1", ] > 0))

ces1_entrez <- bitr(
  ces1_targets,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_results <- enrichGO(
  gene = ces1_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  readable = TRUE
)

n3 <- dotplot(go_results, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment for CES1 Targets") +
  theme_minimal(base_size = 13)

ggsave(file.path(output_dir, "CES1_GO_enrichment.png"), n3, width = 10, height = 7, bg = "white")


#CES1 Source Expression (Tumor)
n5 <- VlnPlot(
  seurat,
  features = "CES1",
  group.by = "celltype",
  pt.size = 0.1
) +
  ggtitle("CES1 Expression Across All Cell Types") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "CES1_source_violin.png"), n5, width = 6, height = 5, bg = "white")


#Plot the top targets with expresiion levels
top_targets <- c("IL1B", "CSF1", "CXCL2", "PDGFB", "IL17A", "TGFB1", "TNFSF15")

data_vln <- FetchData(seurat, vars = c("celltype", top_targets)) %>%
  pivot_longer(-celltype, names_to = "gene", values_to = "expression")

n3 <- ggplot(data_vln, aes(x = celltype, y = expression)) +
  geom_violin(scale = "width", fill = "steelblue") +
  facet_wrap(~ gene, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    strip.placement = "outside"
  )

ggsave(file.path(output_dir, "CES1_target_voilin_plot.png"), n3, width = 7, height=6, bg = "white")

# CES1 UMAP Expression Plot
# Normalize and find variable features (if not already done)
library(future)
plan(sequential)  # switch to safe mode
options(future.globals.maxSize = 10 * 1024^3)  # 10GB

seurat <- NormalizeData(seurat)
seurat<- FindVariableFeatures(seurat)

# Run PCA (if not already done)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

# Run UMAP (this is what you're missing)
seurat <- RunUMAP(seurat, dims = 1:10)  # You can adjust dims as needed
library(viridis)

# Rebuild FeaturePlot to separate low vs high expression
n4 <- FeaturePlot(
  seurat,
  features = "CES1",
  reduction = "umap",
  pt.size = 0.2,
  order = TRUE
) +
  scale_color_gradientn(
    colours = c("grey90", "yellow", "red", "blue"),
    values = scales::rescale(c(0, 0.1, 1, 3)),
    name = "CES1 Expression",
    na.value = "grey90"
  ) +
  labs(title = "UMAP Expression of CES1") +
  theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "CES1_UMAP_expression.png"), n4, width = 6, height = 5, bg = "white")

#Heatmap of CES1 Top Target Genes
# Plot heatmap
# Select top 10 ligands that target CES1
top_ces1_ligands <- sort(ligand_target_matrix[, "CES1"], decreasing = TRUE)
top_ces1_ligands <- names(top_ces1_ligands[top_ces1_ligands > 0])[1:10]  # top 10 non-zero
# Scale before heatmap
# Scale before heatmap
seurat <- ScaleData(seurat, features = top_ces1_ligands)

# Heatmap plot
n_heat <- DoHeatmap(seurat, features = top_ces1_ligands, group.by = "celltype", slot = "data") +
  scale_fill_gradientn(
    colours = c("grey90", "skyblue", "darkblue"),
    values = scales::rescale(c(0, 1, 2.5)),
    name = "Expression"
  ) +
  ggtitle("Top Ligands Regulating CES1 Expression") +
  theme_minimal(base_size = 13)

ggsave(file.path(output_dir, "CES1_top_ligands_heatmap.png"), n_heat, width = 8, height = 6, bg = "white")


# Visualize top ligands ranked
n6 <- ggplot(ligand_activities_combined, aes(x = reorder(test_ligand, pearson), y = pearson, fill = receiver)) +
  geom_col(position = "dodge") +
  labs(
    title = "Ranked Ligand Activities per Receiver",
    x = "Ligand", y = "Pearson Score"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8)) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(output_dir, "CES1_ligand_rank_plot.png"), n6, width = 13, height = 8, bg = "white")

message("NicheNet analysis complete on celltype. Output saved in output_dir.")


