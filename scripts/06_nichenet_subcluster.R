# ============================================================================
# File: scripts/NicheNet_detailed_subcluster.R
# ============================================================================

# ------------------------------ REQUIRED LIBS -------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(nichenetr)
  library(tidyverse)
  library(future)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(Matrix)          # rowMeans on sparse
  library(ComplexHeatmap)  # Heatmap(), draw()
  library(circlize)        # colorRamp2()
  library(gridExtra)       # grid.arrange()
  library(reshape2)        # melt()
  library(igraph)          # graph_from_data_frame
  library(ggraph)          # graph plotting
  library(scales)
  library(grid)            # arrow(), unit(), etc.
})

plan(sequential)
options(future.globals.maxSize = 10 * 1024^3)
options(timeout = 300)

# ------------------------------ TUNABLE PARAMS ------------------------------
# Increase network density per your request
n_top_degs      <- 500   # was 200
n_top_ligands   <- 50    # was 20
edge_threshold  <- 0.05  # was 0.1

# ------------------------------- LOAD OBJECT --------------------------------
load("E:/Desktop/Dissertation/Data/Colon_zhang/output/seurat_10x_logTPM.RData")

# --------------------------- OUTPUT DIRECTORY --------------------------------
output_dir <- "E:/Desktop/Dissertation/Data/Colon_zhang/output/NiecheNet_subclusters"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------- SUBCLUSTER RESOLUTION ---------------------------
celltype_col <- "Sub_Cluster"             # key column used throughout
Idents(seurat) <- factor(seurat@meta.data[[celltype_col]])

# --------------------------- BASIC PREPROCESSING ----------------------------
DefaultAssay(seurat) <- DefaultAssay(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000, verbose = FALSE)
seurat <- ScaleData(seurat, features = VariableFeatures(seurat), verbose = FALSE)
seurat <- RunPCA(seurat, features = VariableFeatures(seurat), verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30, verbose = FALSE)

# ---------------------------- NICHE-NET DATA --------------------------------
ligand_target_matrix <- readRDS(
  url("https://zenodo.org/records/5884439/files/ligand_target_matrix_nsga2r_final.rds")
)
# Ensure columns = ligands, rows = target genes
ligand_target_matrix <- t(ligand_target_matrix)
stopifnot("CES1" %in% colnames(ligand_target_matrix))

ligand_receptor_network <- readRDS(
  url("https://zenodo.org/records/5884439/files/lr_network_human_21122021.rds")
)
stopifnot(nrow(ligand_receptor_network) > 0)

# Convenience alias used later
lr_network <- ligand_receptor_network

# ------------------------- SENDER / RECEIVER PICKERS ------------------------
sender_patterns   <- c("TAM", "Macro", "Mono")
receiver_patterns <- c("CD8", "B", "CD4", "ILC")

all_clusters <- levels(Idents(seurat))
sender_celltype2    <- grep(paste(sender_patterns, collapse = "|"), all_clusters, value = TRUE)
receiver_celltypes2 <- grep(paste(receiver_patterns, collapse = "|"), all_clusters, value = TRUE)

if (length(sender_celltype2) == 0) stop("No sender subclusters matched the pattern")
if (length(receiver_celltypes2) == 0) stop("No receiver subclusters matched the pattern")

# ------------------------ EXPRESSED GENES HELPERS ---------------------------
safe_get_expr_genes <- function(seurat_obj, cluster, pct = 0.10) {
  tryCatch({
    cells <- WhichCells(seurat_obj, idents = cluster)
    if (length(cells) == 0) return(character())
    expr_data <- GetAssayData(seurat_obj, slot = "data")[, cells, drop = FALSE]
    gene_pct <- Matrix::rowMeans(expr_data > 0)
    names(gene_pct)[gene_pct >= pct]
  }, error = function(e) character())
}

expressed_genes_sender_list <- purrr::map(sender_celltype2, ~safe_get_expr_genes(seurat, .x)) %>%
  purrr::discard(~length(.x) == 0)
expressed_genes_sender <- if (length(expressed_genes_sender_list) > 0) Reduce(union, expressed_genes_sender_list) else character()

expressed_genes_receiver_list <- purrr::map(receiver_celltypes2, ~safe_get_expr_genes(seurat, .x)) %>%
  purrr::discard(~length(.x) == 0)
expressed_genes_receiver <- if (length(expressed_genes_receiver_list) > 0) Reduce(intersect, expressed_genes_receiver_list) else character()

if (length(expressed_genes_receiver) == 0 || length(expressed_genes_sender) == 0)
  stop("No expressed genes found in receiver or sender subclusters")

# ------------------------------ DEGS (example) ------------------------------
example_receiver <- receiver_celltypes2[1]
deg_receiver <- FindMarkers(
  seurat, ident.1 = example_receiver, group.by = celltype_col,
  only.pos = TRUE, logfc.threshold = 0.25
) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice_head(n = n_top_degs) %>%
  dplyr::pull(gene)

# ------------------------- LIGAND ACTIVITIES (NN) ---------------------------
ligands2 <- unique(ligand_receptor_network$from)
ligands_sender2 <- intersect(ligands2, expressed_genes_sender)

#This is the core NicheNet function. It scores ligands based on how well their 
#predicted target genes overlap with my CES1-related genes in receiver cells
ligand_activities2 <- predict_ligand_activities(
  geneset = deg_receiver,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligands_sender2
)

# ------------------------------- HEATMAP ------------------------------------
top_ligands <- ligand_activities2 %>%
  dplyr::arrange(desc(pearson)) %>%
  dplyr::slice_head(n = n_top_ligands) %>%
  dplyr::pull(test_ligand)
valid_top_ligands  <- intersect(top_ligands, colnames(ligand_target_matrix))   # ligands are COLUMNS
valid_deg_receiver <- intersect(deg_receiver, rownames(ligand_target_matrix))  # targets are ROWS
if (length(valid_top_ligands) == 0 || length(valid_deg_receiver) == 0) stop("No overlap between top ligands/DEGs and matrix")

# For visualization: ligands x targets
heatmap_targets <- t(ligand_target_matrix[valid_deg_receiver, valid_top_ligands, drop = FALSE])

png(file.path(output_dir, "nichenet_heatmap.png"), width = 1200, height = 1000)
if (nrow(heatmap_targets) > 0 && ncol(heatmap_targets) > 0) {
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(heatmap_targets, name = "Regulation"))
} else {
  message("No data for heatmap after filtering")
}
dev.off()

# ------------------------------- DOTPLOTS -----------------------------------
# Ligand Expression DotPlot (top 10 ligands)
top_ligands_plot <- top_ligands[1:min(10, length(top_ligands))]

png(file.path(output_dir, "LigandExpression_DotPlot.png"), width = 2000, height = 1000, res = 150)
DotPlot(seurat, features = top_ligands_plot, group.by = celltype_col) +
  RotatedAxis() + ggtitle("Ligand Expression Across Subclusters")
dev.off()

# Receptor Expression DotPlot
receptors <- lr_network$to %>% unique() %>% intersect(rownames(seurat))
receptor_subset <- head(receptors, 20)

png(file.path(output_dir, "ReceptorExpression_DotPlot.png"), width = 2000, height = 1000, res = 150)
DotPlot(seurat, features = receptor_subset, group.by = celltype_col) +
  RotatedAxis() + ggtitle("Receptor Expression Across Subclusters")
dev.off()

# -------------------------- CES1-RELATED PLOTS ------------------------------
# CES1 violin across subclusters
png(file.path(output_dir, "CES1_ViolinPlot.png"), width = 2000, height = 1000, res = 150)
VlnPlot(seurat, features = "CES1", group.by = celltype_col, pt.size = 0.1) +
  ggtitle("CES1 Expression Across Subclusters") +
  theme(legend.position = "none", plot.margin = unit(c(1,1,1,2), "cm"))
dev.off()

# CES1 correlation bars (global, like your original)
ces1_expr <- FetchData(seurat, vars = "CES1")$CES1
expr_matrix <- GetAssayData(seurat, slot = "data")
expr_matrix <- expr_matrix[Matrix::rowMeans(expr_matrix > 0) >= 0.1, ]
expr_scaled <- t(scale(t(expr_matrix), center = TRUE, scale = TRUE))
ces1_expr <- scale(ces1_expr, center = TRUE, scale = TRUE)
cor_scores <- Matrix::rowMeans(expr_scaled * rep(ces1_expr, each = nrow(expr_scaled)))

top_pos <- head(sort(cor_scores, decreasing = TRUE), 20)
top_neg <- head(sort(cor_scores, decreasing = FALSE), 20)

base_theme <- theme_minimal(base_size = 18) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    plot.title  = element_text(size = 20, face = "bold"),
    plot.margin = ggplot2::margin(10, 30, 10, 10, unit = "pt")
  )

g1 <- ggplot(tibble(Gene = names(top_pos), Correlation = top_pos), aes(x = reorder(Gene, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() + ggtitle("Top Positive Correlated Genes with CES1") + base_theme

g2 <- ggplot(tibble(Gene = names(top_neg), Correlation = top_neg), aes(x = reorder(Gene, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "firebrick") + coord_flip() + ggtitle("Top Negative Correlated Genes with CES1") + base_theme

png(file.path(output_dir, "CES1_Correlation_Bars.png"), width = 2800, height = 2200, res = 220)
grid.arrange(g1, g2, ncol = 1, heights = c(1, 1))
dev.off()

# CES1 UMAP feature plot (dims 1:30 already available)
png(file.path(output_dir, "CES1_FeaturePlot.png"), width = 1000, height = 800)
FeaturePlot(seurat, features = "CES1", reduction = "umap",
            cols = c("grey90", "lightblue", "blue", "darkblue"), order = TRUE)
dev.off()

# -------------------------- MORE NICHE-NET PLOTS ----------------------------
# DotPlot for top ligands again (sender-focused view)
top10_ligands <- top_ligands[1:min(10, length(top_ligands))]
png(file.path(output_dir, "TopLigands_DotPlot.png"), width = 2000, height = 1000, res = 150)
DotPlot(seurat, features = top10_ligands, group.by = celltype_col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Activity scatter for top ligands
ligand_plot_df <- ligand_activities2 %>% dplyr::filter(test_ligand %in% top_ligands) %>%
  dplyr::mutate(rank = rank(-pearson)) %>% dplyr::arrange(rank)

png(file.path(output_dir, "TopLigand_Activity_Plot.png"), width = 1600, height = 800, res = 150)
print(ggplot(ligand_plot_df, aes(x = reorder(test_ligand, -pearson), y = pearson)) +
        geom_point(size = 6, color = "steelblue") + theme_minimal() +
        labs(x = "Ligand", y = "Pearson Correlation", title = "Top Ligands by Activity (NicheNet)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# Receptor expression in receiver cells
top_ligand_receptors <- lr_network %>% dplyr::filter(from %in% top_ligands) %>% dplyr::pull(to) %>% unique()
top_receptors <- head(intersect(top_ligand_receptors, rownames(seurat)), 10)

png(file.path(output_dir, "TopReceptors_DotPlot.png"), width = 2000, height = 1000, res = 150)
DotPlot(seurat, features = top_receptors, group.by = celltype_col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# UMAP labeled by Sub_Cluster
png(file.path(output_dir, "CES1_by_SubCluster_UMAP.png"), width = 1000, height = 800)
DimPlot(seurat, group.by = celltype_col, label = TRUE, reduction = "umap") +
  ggtitle("UMAP colored by Sub_Cluster")
dev.off()

# Sender-Receiver Interaction Heatmap (counts of expressed LR pairs)
interaction_matrix <- matrix(0, nrow = length(sender_celltype2), ncol = length(receiver_celltypes2),
                             dimnames = list(sender_celltype2, receiver_celltypes2))
for (sender in sender_celltype2) {
  for (receiver in receiver_celltypes2) {
    expressed_ligands_sender   <- safe_get_expr_genes(seurat, sender)
    expressed_receptors_receiver <- safe_get_expr_genes(seurat, receiver)
    lr_pairs <- lr_network %>% dplyr::filter(from %in% expressed_ligands_sender, to %in% expressed_receptors_receiver)
    interaction_matrix[sender, receiver] <- nrow(lr_pairs)
  }
}

png(file.path(output_dir, "Sender_Receiver_Interaction_Heatmap.png"), width = 1200, height = 1000)
ComplexHeatmap::Heatmap(interaction_matrix, name = "LR Pairs",
                        col = circlize::colorRamp2(c(0, max(interaction_matrix)), c("white", "darkred"))) %>%
  ComplexHeatmap::draw()
dev.off()

# Ligand-Target Bipartite Graph (threshold to declutter)
ligand_target_edges <- reshape2::melt(as.matrix(heatmap_targets)) %>%
  dplyr::filter(value > edge_threshold) %>% dplyr::rename(from = Var1, to = Var2)

if (nrow(ligand_target_edges) > 0) {
  g_bipartite <- igraph::graph_from_data_frame(ligand_target_edges, directed = TRUE)
  igraph::V(g_bipartite)$type <- ifelse(igraph::V(g_bipartite)$name %in% rownames(heatmap_targets), "Ligand", "Target")
  png(file.path(output_dir, "Ligand_Target_BipartiteGraph.png"), width = 1800, height = 1400, res = 150)
  print(
    ggraph(g_bipartite, layout = "fr") +
      geom_edge_link(aes(alpha = 0.5), arrow = grid::arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) +
      geom_node_point(aes(color = type), size = 5) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
      scale_color_manual(values = c("Ligand" = "dodgerblue", "Target" = "orange")) +
      theme_void() + ggtitle("Ligand to Target Network")
  )
  dev.off()
}

# Target Gene Module Score in Receiver Cells
if (length(valid_deg_receiver) >= 5) {
  seurat <- AddModuleScore(seurat, features = list(valid_deg_receiver), name = "TargetModule")
  png(file.path(output_dir, "Receiver_Target_ModuleScore_UMAP.png"), width = 1000, height = 800)
  FeaturePlot(seurat, features = "TargetModule1", reduction = "umap",
              cols = c("lightgrey", "blue", "darkblue"), order = TRUE) +
    ggtitle("Target Gene Module Score in UMAP")
  dev.off()
}

# Barplot: Ligand regulates how many DEGs (per top ligand set)
ligand_target_counts <- rowSums(heatmap_targets > edge_threshold, na.rm = TRUE)

ligand_target_df <- tibble(
  Ligand = rownames(heatmap_targets),
  Target_Count = as.numeric(ligand_target_counts)
) %>%
  dplyr::filter(Target_Count > 0) %>%                      # hide empty bars
  dplyr::arrange(Target_Count)

png(file.path(output_dir, "Ligand_Target_Counts_Barplot.png"), width = 2000, height = 1200, res = 180)
print(
  ggplot(ligand_target_df, aes(x = forcats::fct_reorder(Ligand, Target_Count), y = Target_Count)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal(base_size = 16) +
    labs(title = "Number of Regulated Target Genes per Ligand",
         x = "Ligand", y = "# of Target Genes")
)
dev.off()

# Export example tables
lr_network_export <- ligand_receptor_network %>% dplyr::filter(from %in% top_ligands & to %in% expressed_genes_receiver)
readr::write_csv(ligand_activities2, file.path(output_dir, "ligand_activities_exampleReceiver.csv"))
readr::write_csv(lr_network_export,  file.path(output_dir, "top_ligand_receptors.csv"))

message("Detailed NicheNet (subcluster) complete. Outputs saved in: ", output_dir)


#==============================
  # ================== CES1: counts + normalized (run at the end) ==================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(broom.mixed)
})

sub_col <- "Sub_Cluster"   # change if your column name differs
stopifnot(sub_col %in% colnames(seurat@meta.data))

# Try to find a per-sample column
sample_col <- if ("Sample" %in% names(seurat@meta.data)) "Sample" else
  if ("Patient" %in% names(seurat@meta.data)) "Patient" else NA

# -------- A) Counts and raw violin (all cells) --------
df_all <- FetchData(seurat, vars = c("CES1", sub_col))
colnames(df_all) <- c("CES1","Sub")

counts <- seurat@meta.data %>%
  dplyr::count(Sub = .data[[sub_col]], name = "cells") %>%   # make a Sub column, then count
  dplyr::arrange(desc(cells))

# order subclusters by size and label with n
df_all$Sub <- factor(df_all$Sub, levels = counts$Sub)
lab_map <- setNames(paste0(counts$Sub, " (n=", counts$cells, ")"), counts$Sub)

p_counts <- ggplot(counts, aes(x = Sub, y = cells)) +
  geom_col(fill = "grey60") +
  geom_text(aes(label = cells), vjust = -0.4, size = 3) +
  labs(title = "Cells per Subcluster", x = NULL, y = "Cells") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(5, 15, 5, 5, unit = "pt")
  )

p_violin <- ggplot(df_all, aes(x = Sub, y = CES1)) +
  geom_violin(fill = "steelblue", alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.15, size = 0.15, alpha = 0.2) +
  scale_x_discrete(labels = lab_map) +
  labs(
    title = "CES1 Expression Across Subclusters",
    x = "Subcluster (n cells)", 
    y = "CES1 (normalized)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(5, 15, 5, 5, unit = "pt")
  )

# save side-by-side: counts | violin
(p_counts | p_violin) +
  plot_layout(widths = c(1, 2))
ggsave(file.path(output_dir, "CES1_counts_and_violin_subclusters.png"),
       width = 16, height = 7, dpi = 300, bg = "white")

# -------- B) Per-sample normalized comparison + stats --------
if (!is.na(sample_col) && sample_col %in% names(seurat@meta.data)) {
  df_norm <- tibble::tibble(
    SampleID = seurat@meta.data[[sample_col]],
    Sub      = seurat@meta.data[[sub_col]],
    CES1     = FetchData(seurat, "CES1")[,1]
  ) %>% tidyr::drop_na(SampleID, Sub)
  
  # per-sample mean CES1 within each subcluster
  ces1_by_sample <- df_norm %>%
    group_by(SampleID, Sub) %>%
    summarise(n_cells = n(),
              mean_CES1 = mean(CES1, na.rm = TRUE),
              .groups = "drop")
  
  # plot: per-sample points + mean±CI per subcluster
  p_norm <- ggplot(ces1_by_sample, aes(x = Sub, y = mean_CES1)) +
    geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 3.2, color = "black") +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.25, color = "black") +
    labs(title = paste0("Per‑sample mean CES1 by Subcluster (normalized across ", sample_col, ")"),
         x = "Subcluster", y = "Mean CES1 per sample") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank())
  
  ggsave(file.path(output_dir, "CES1_per_sample_mean_by_subcluster.png"),
         p_norm, width = 14, height = 6, dpi = 300, bg = "white")
  
  # ANOVA & Kruskal on per-sample means
  anova_mean <- aov(mean_CES1 ~ Sub, data = ces1_by_sample)
  capture.output(summary(anova_mean),
                 file = file.path(output_dir, "ANOVA_per_sample_mean_CES1_by_subcluster.txt"))
  
  kw_mean <- kruskal.test(mean_CES1 ~ Sub, data = ces1_by_sample)
  writeLines(capture.output(kw_mean),
             con = file.path(output_dir, "Kruskal_per_sample_mean_CES1_by_subcluster.txt"))
  
  # Mixed-effects (single-cell, random intercept per sample)
  lmm <- lmer(CES1 ~ Sub + (1|SampleID), data = df_norm)
  tidy_lmm <- broom.mixed::tidy(lmm, effects = "fixed", conf.int = TRUE)
  readr::write_csv(tidy_lmm,
                   file.path(output_dir, "LMM_single_cell_CES1_by_subcluster.csv"))
} else {
  message("Skipped normalized plot: no 'Sample' or 'Patient' column found in meta.data.")
}
# Ensure output folder exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Save counts and violin separately (optional) ---
ggsave(file.path(output_dir, "CES1_counts_subclusters.png"),
       plot = p_counts, width = 8, height = 4, dpi = 300, bg = "white")

ggsave(file.path(output_dir, "CES1_violin_subclusters.png"),
       plot = p_violin, width = 12, height = 6, dpi = 300, bg = "white")

# --- Save the combined panel explicitly ---
panel_counts_violin <- (p_counts | p_violin) + patchwork::plot_layout(widths = c(1, 2))

ggsave(file.path(output_dir, "CES1_counts_and_violin_subclusters.png"),
       plot = panel_counts_violin, width = 16, height = 7, dpi = 300, bg = "white")

# --- Save the per‑sample normalized plot explicitly (only if it exists) ---
if (exists("p_norm")) {
  ggsave(file.path(output_dir, "CES1_per_sample_mean_by_subcluster.png"),
         plot = p_norm, width = 14, height = 6, dpi = 300, bg = "white")
}

#new normalization CES1 barplot and voilin plot with CES1 count
ces1_threshold <- 0.1  # log(TPM+1) threshold

# Make sure 'Tissue' or equivalent column exists
# Replace "Tissue" with your grouping variable (e.g., Sub_Cluster, Condition, etc.)
group_col <- "Tissue"  
stopifnot(group_col %in% colnames(seurat@meta.data))

# Fetch CES1 and group info
ces1_all <- FetchData(seurat, vars = c("CES1", group_col))
colnames(ces1_all) <- c("CES1", "Group")

# Summary: counts and % CES1+
ces1_summary <- ces1_all %>%
  group_by(Group) %>%
  summarise(
    total_cells = n(),
    ces1_pos    = sum(CES1 > ces1_threshold, na.rm = TRUE),
    perc_pos    = 100 * ces1_pos / total_cells,
    .groups     = "drop"
  )

# LEFT plot: CES1+ cell counts with %
p_ces1_counts_pct <- ggplot(ces1_summary, aes(x = Group, y = ces1_pos, fill = Group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", ces1_pos, perc_pos)),
            vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = sprintf("CES1+ cells per %s (> %.2f)", group_col, ces1_threshold),
       x = NULL, y = "CES1-positive cells (count)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# Filter for CES1-positive cells
ces1_filtered <- ces1_all %>% filter(CES1 > ces1_threshold)

# Save removal stats
filter_stats <- ces1_all %>%
  group_by(Group) %>%
  summarise(
    removed       = sum(CES1 <= ces1_threshold, na.rm = TRUE),
    kept          = sum(CES1 > ces1_threshold, na.rm = TRUE),
    total         = n(),
    perc_removed  = 100 * removed / total,
    .groups       = "drop"
  )
write.csv(filter_stats, file.path(output_dir, "CES1_filtering_stats_nichenet.csv"), row.names = FALSE)

# RIGHT plot: CES1 expression (filtered)
p_ces1_expr_filtered <- ggplot(ces1_filtered, aes(Group, CES1, fill = Group)) +
  geom_violin(width = 0.9, alpha = 0.75, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(aes(color = Group), width = 0.15, size = 0.3, alpha = 0.18, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(title = sprintf("CES1 expression per %s (> %.2f)", group_col, ces1_threshold),
       x = NULL, y = "CES1 (log(TPM+1))") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# Combine panels
panel_counts_pct_expr <- p_ces1_counts_pct | p_ces1_expr_filtered

# Save plot
ggsave(file.path(output_dir, "CES1_counts_pct_and_expression_filtered_nichenet.png"),
       panel_counts_pct_expr, width = 12, height = 5, dpi = 300, bg = "white")


# (Optional) a quick confirmation
message("Saved plots to: ", normalizePath(output_dir, winslash = "/"))

# ============================================================================ 

# === CES1 per Sub_Cluster: separate figures (counts + violin) ===
library(forcats)
library(stringr)

ces1_threshold <- 0.10
sub_col <- "Sub_Cluster"
stopifnot(sub_col %in% colnames(seurat@meta.data))

# Pull CES1 + subcluster
df <- FetchData(seurat, vars = c("CES1", sub_col))
colnames(df) <- c("CES1","Sub")

# Summaries
sum_tbl <- df %>%
  dplyr::group_by(Sub) %>%
  dplyr::summarise(
    total_cells = dplyr::n(),
    ces1_pos    = sum(CES1 > ces1_threshold, na.rm = TRUE),
    perc_pos    = 100 * ces1_pos / total_cells,
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(ces1_pos))    # order by CES1+ count (change to total_cells if you prefer)

# Keep only CES1+ for violin; preserve ordering
df_pos <- df %>%
  dplyr::filter(CES1 > ces1_threshold) %>%
  dplyr::mutate(Sub = factor(Sub, levels = sum_tbl$Sub))

# ---- 1) Horizontal bar: CES1+ counts per subcluster ----
p_counts <- ggplot(sum_tbl, aes(y = forcats::fct_reorder(Sub, ces1_pos), x = ces1_pos)) +
  geom_col(fill = "grey60") +
  geom_text(aes(label = sprintf("%d (%.1f%%)", ces1_pos, perc_pos)),
            hjust = -0.05, size = 4) +
  labs(title = sprintf("CES1+ cells per Subcluster (> %.2f)", ces1_threshold),
       x = "CES1-positive cells (count)", y = NULL) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 12),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(labels = function(x) stringr::str_replace_all(x, "_", "\n"))

# ---- 2) Violin: CES1 expression in CES1+ cells ----
p_violin <- ggplot(df_pos, aes(x = Sub, y = CES1)) +
  geom_violin(fill = "steelblue", alpha = .85, color = NA, trim = FALSE) +
  geom_boxplot(width = .16, outlier.shape = NA, alpha = .95) +
  geom_jitter(width = .12, size = .3, alpha = .15) +
  labs(title = sprintf("CES1 expression per Subcluster (> %.2f)", ces1_threshold),
       x = "Subcluster", y = "CES1 (log(TPM+1))") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = function(x) stringr::str_replace_all(x, "_", "\n"))

# ---- Save each separately; canvas scales with #subclusters ----
n_sub   <- nrow(sum_tbl)
h_count <- max(10, n_sub * 0.35)   # more rows -> taller bar chart
w_count <- 16
w_violin <- max(18, n_sub * 0.7)   # many x labels -> wider violin plot
h_violin <- 10

ggsave(file.path(output_dir, "CES1_counts_by_subcluster.png"),
       p_counts, width = w_count, height = h_count, dpi = 400, bg = "white", limitsize = FALSE)
ggsave(file.path(output_dir, "CES1_counts_by_subcluster.pdf"),
       p_counts, device = cairo_pdf, width = w_count, height = h_count, bg = "white")

ggsave(file.path(output_dir, "CES1_expression_filtered_by_subcluster.png"),
       p_violin, width = w_violin, height = h_violin, dpi = 400, bg = "white", limitsize = FALSE)
ggsave(file.path(output_dir, "CES1_expression_filtered_by_subcluster.pdf"),
       p_violin, device = cairo_pdf, width = w_violin, height = h_violin, bg = "white")



# ============================================================================
# ML comparison on NicheNet: TAM-SPP1 vs TAM-C1QC (robust RF version)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(randomForest)
  library(ggplot2)
  library(tidyr)
})

ml_output_dir <- "E:/Desktop/Dissertation/Data/Colon_zhang/output/Niechenet_subcluster_ML"
dir.create(ml_output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------- 1. Extract interactions --------------------------
tam_senders <- c("hM13_TAM-SPP1", "hM12_TAM-C1QC")

interaction_list <- list()

for (sender in tam_senders) {
  # Get ligands expressed in this sender
  ligs <- safe_get_expr_genes(seurat, sender)
  expressed_ligs <- intersect(ligs, lr_network$from)
  
  for (receiver in receiver_celltypes2) {
    rec_genes <- safe_get_expr_genes(seurat, receiver)
    expressed_recs <- intersect(rec_genes, lr_network$to)
    
    lr_pairs <- lr_network %>%
      filter(from %in% expressed_ligs, to %in% expressed_recs) %>%
      mutate(Sender = sender, Receiver = receiver)
    
    if (nrow(lr_pairs) > 0) {
      # Include pathway if available in lr_network
      if ("pathway_name" %in% colnames(lr_pairs)) {
        lr_pairs$feature <- paste(lr_pairs$from,
                                  lr_pairs$to,
                                  lr_pairs$pathway_name,
                                  lr_pairs$Receiver,
                                  sep = "_")
      } else {
        # fallback: no pathway column
        lr_pairs$feature <- paste(lr_pairs$from,
                                  lr_pairs$to,
                                  lr_pairs$Receiver,
                                  sep = "_")
      }
      interaction_list[[paste(sender, receiver, sep = "_")]] <- lr_pairs
    }
  }
}

interaction_df <- bind_rows(interaction_list)

# --------------------- 2. Define groups (SPP1 vs C1QC) ----------------------
interaction_df$group <- ifelse(grepl("SPP1", interaction_df$Sender), "SPP1", "C1QC")

# --------------------- 3. Unique vs common interactions ---------------------
spp1_feats <- unique(interaction_df$feature[interaction_df$group == "SPP1"])
c1qc_feats <- unique(interaction_df$feature[interaction_df$group == "C1QC"])

common_feats <- intersect(spp1_feats, c1qc_feats)
spp1_unique <- setdiff(spp1_feats, common_feats)
c1qc_unique <- setdiff(c1qc_feats, common_feats)

write.csv(data.frame(Interaction = common_feats, Group = "Common"),
          file.path(ml_output_dir, "Common_interactions.csv"), row.names = FALSE)
write.csv(data.frame(Interaction = spp1_unique, Group = "SPP1"),
          file.path(ml_output_dir, "SPP1_unique_interactions.csv"), row.names = FALSE)
write.csv(data.frame(Interaction = c1qc_unique, Group = "C1QC"),
          file.path(ml_output_dir, "C1QC_unique_interactions.csv"), row.names = FALSE)

# -------------------- 4. Encode into ML matrix ------------------------------
ml_encoded <- dcast(interaction_df, group ~ feature, 
                    value.var = "feature", fun.aggregate = length, fill = 0)

ml_encoded$group <- factor(ml_encoded$group)

# Sanitize column names for RF
colnames(ml_encoded) <- make.names(colnames(ml_encoded))

name_map <- data.frame(
  Safe = make.names(c(spp1_unique, c1qc_unique, common_feats)),
  Original = c(spp1_unique, c1qc_unique, common_feats),
  stringsAsFactors = FALSE
)

# -------------------- 5. Random Forest model --------------------------------
set.seed(123)
rf_model <- randomForest(group ~ ., data = ml_encoded, importance = TRUE, ntree = 500)
print(rf_model)

imp <- importance(rf_model)
imp_df <- data.frame(Feature = rownames(imp), imp) %>%
  left_join(name_map, by = c("Feature" = "Safe")) %>%
  mutate(DisplayName = ifelse(!is.na(Original), Original, Feature)) %>%
  arrange(desc(MeanDecreaseGini))

write.csv(imp_df, file.path(ml_output_dir, "RF_feature_importance.csv"), row.names = FALSE)

# -------------------- 6. Plot unique interactions ---------------------------
# ======================================================
# 6. Plot unique interactions (split into PNGs)
# ======================================================

# 6a. Keep only unique features (SPP1-only or C1QC-only)
unique_feats <- c(spp1_unique, c1qc_unique)

imp_unique <- imp_df %>%
  filter(DisplayName %in% unique_feats) %>%
  mutate(Group = ifelse(DisplayName %in% spp1_unique, "SPP1", "C1QC"))

# 6b. Filter out zero/negative importance
imp_unique <- imp_unique %>% filter(MeanDecreaseGini > 0)

# 6c. Order globally by importance
imp_unique <- imp_unique %>%
  arrange(desc(MeanDecreaseGini))

# 6d. Split into chunks
chunk_size <- 80
n_chunks <- ceiling(nrow(imp_unique) / chunk_size)

for (i in 1:n_chunks) {
  imp_chunk <- imp_unique[((i-1)*chunk_size + 1) : min(i*chunk_size, nrow(imp_unique)), ]
  
  p_chunk <- ggplot(imp_chunk, aes(x = reorder(DisplayName, MeanDecreaseGini),
                                   y = MeanDecreaseGini, fill = Group)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Unique interactions (chunk", i, "of", n_chunks, ")"),
         x = "Interaction (Ligand_Receptor_Pathway_Receiver)",
         y = "Importance (MeanDecreaseGini)") +
    theme_bw(base_size = 12) +
    scale_fill_manual(values = c("SPP1" = "darkorange", "C1QC" = "steelblue"))
  
  ggsave(file.path(ml_output_dir, paste0("RF_unique_feature_importance_chunk", i, ".png")),
         p_chunk, width = 18, height = 13, dpi = 300)
}


# -------------------- 7. Receiver summary (bias) -----------------------------
get_receiver <- function(x) {
  sapply(strsplit(x, "_"), function(parts) tail(parts, 1))
}

imp_unique_nonzero <- imp_unique %>%
  filter(MeanDecreaseGini > 0) %>%
  mutate(Receiver = get_receiver(DisplayName))

receiver_summary <- imp_unique_nonzero %>%
  group_by(Group, Receiver) %>%
  summarise(UniqueInteractions = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = UniqueInteractions, values_fill = 0) %>%
  mutate(Difference = C1QC - SPP1)

write.csv(receiver_summary,
          file.path(ml_output_dir, "Receiver_summary_SPP1_vs_C1QC.csv"),
          row.names = FALSE)

# Receiver bias diverging barplot
p2 <- ggplot(receiver_summary, 
             aes(x = reorder(Receiver, Difference), y = Difference,
                 fill = ifelse(Difference > 0, "C1QC", "SPP1"))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("SPP1" = "darkorange", "C1QC" = "steelblue")) +
  labs(title = "Receiver bias (C1QC vs SPP1, NicheNet unique interactions)",
       x = "Receiver cell type", 
       y = "Difference (C1QC − SPP1)") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(file.path(ml_output_dir, "Receiver_bias_C1QC_vs_SPP1.png"),
       p2, width = 10, height = 6, dpi = 300)


# -------------------- 8. Raw counts of unique interactions by receiver -------
receiver_summary_raw <- imp_unique_nonzero %>%
  group_by(Group, Receiver) %>%
  summarise(UniqueInteractions = n(), .groups = "drop")

# Raw counts of unique interactions by receiver
p3 <- ggplot(receiver_summary_raw, 
             aes(x = Receiver, y = UniqueInteractions, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Unique interactions by receiver cell type (NicheNet, ML-driven)",
       x = "Receiver cell type", 
       y = "Number of unique interactions") +
  theme_bw() +
  scale_fill_manual(values = c("SPP1" = "darkorange", "C1QC" = "steelblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

ggsave(file.path(ml_output_dir, "Unique_interactions_by_receiver.png"),
       plot = p3, width = 10, height = 6, dpi = 300)



message("ML on NicheNet (CellChat-style) complete: results in ", ml_output_dir)

