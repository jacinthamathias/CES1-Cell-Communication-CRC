# file: scripts/cellchat_subcluster_full_pipeline_ces1.R

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(gridExtra)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(scales)
})

# ---- Config ----
output_dir <- "E:/Desktop/Dissertation/Data/Colon_zhang/output/CellChat_subcluster"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
obj_path <- "E:/Desktop/Dissertation/Data/Colon_zhang/output/seurat_10x_logTPM.RData"

# ---- Load Seurat ----
load(obj_path)  # expects object named `seurat`
stopifnot(exists("seurat"))
stopifnot("Sub_Cluster" %in% colnames(seurat@meta.data))
seurat$Sub_Cluster <- droplevels(factor(seurat$Sub_Cluster))
Idents(seurat) <- seurat$Sub_Cluster

# ---- Build CellChat on Sub_Cluster ----
# why: align communications with the same grouping used for CES1 averages
if (!exists("CellChatDB.human")) data(CellChatDB.human, package = "CellChat")
if (!exists("PPI.human")) data(PPI.human, package = "CellChat")

cellchat <- createCellChat(object = seurat, group.by = "Sub_Cluster", assay = "scRNA")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#This function calculates the probability of communication
#between cell groups based on ligand–receptor expression
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#This function takes all the individual ligand–receptor
#pairs and aggregates them into pathway-level or network-level communication.
cellchat <- aggregateNet(cellchat)

# communications table
comm_df <- subsetCommunication(cellchat)

# group sizes for plotting
groupSize <- as.numeric(table(cellchat@idents))

# ---- 1. Circle plots (count & strength) ----
png(file.path(output_dir, "Interaction_count_circle.png"), width = 5000, height=5000, res = 300)
par(mar = c(1,1,5,1), xpd = TRUE)
#This is a visualization function. It plots circle diagrams where nodes 
#are cell groups and edges are communication links
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of Interactions")
dev.off()

png(file.path(output_dir, "Interaction_Strength.png"), width = 5000, height=5000, res = 300)
par(mar = c(1,1,5,1), xpd = TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Interaction Strength")
dev.off()

# ---- 2. Heatmaps ----
png(file.path(output_dir, "Heatmap_count.png"), width = 3000, height = 2000, res = 300)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

png(file.path(output_dir, "Heatmap_strength.png"), width = 3000, height = 2000, res = 300)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

# ---- 3. Bubble plot (Top 30 pathways) ----
# file: scripts/single_summary_topK_interactions.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

# ---- Params ----
top_k      <- 30         # number of interactions to show
p_cut      <- 0.05       # set to 1 to keep all
wrap_width <- 65         # y-label wrapping
outfile    <- file.path(output_dir, sprintf("Top%d_interactions_single_summary.png", top_k))

stopifnot(exists("cellchat"), dir.exists(output_dir))

# ---- Data ----
comm_df <- subsetCommunication(cellchat)

# detect p column (why: CellChat versions differ)
p_col <- intersect(c("pval.adj", "pval"), colnames(comm_df))[1]
if (is.na(p_col)) message("[Info] No p-value column detected; size will be constant.")

# pick LR label column
int_col <- intersect(c("interaction_name_2", "interaction_name"), colnames(comm_df))[1]
stopifnot(!is.na(int_col))

# filter by p if available
comm_use <- if (!is.na(p_col)) {
  dplyr::filter(comm_df, .data[[p_col]] < p_cut)
} else {
  comm_df
}

# build labels and -log10 p
comm_use <- comm_use %>%
  mutate(
    lr_label   = .data[[int_col]],
    st_label   = paste0(source, " \u2192 ", target),
    row_label  = paste0(pathway_name, " | ", lr_label, " | ", st_label),
    neglog10_p = if (!is.na(p_col)) -log10(pmax(.data[[p_col]], .Machine$double.xmin)) else NA_real_
  )

# aggregate duplicates
comm_agg <- comm_use %>%
  group_by(row_label, pathway_name) %>%
  summarise(
    prob       = mean(prob, na.rm = TRUE),
    neglog10_p = suppressWarnings(max(neglog10_p, na.rm = TRUE)),
    .groups = "drop"
  )

comm_agg$neglog10_p[!is.finite(comm_agg$neglog10_p)] <- NA_real_

# rank and keep Top-K
top_tbl <- comm_agg %>%
  arrange(desc(prob)) %>%
  slice_head(n = top_k) %>%
  mutate(
    row_label_wrapped = stringr::str_wrap(row_label, width = wrap_width),
    row_label_wrapped = forcats::fct_reorder(row_label_wrapped, prob),
    size_plot = ifelse(is.na(neglog10_p), 3, neglog10_p)  # older ggplot2 compatibility
  )

if (nrow(top_tbl) == 0) stop("No interactions after filtering; relax p_cut or increase top_k.")

# plot
p_single <- ggplot(top_tbl, aes(x = prob, y = row_label_wrapped,
                                color = pathway_name, size = size_plot)) +
  geom_point() +
  scale_size_continuous(name = expression(-log[10](p)), range = c(2.2, 6)) +
  labs(
    title = sprintf("Top %d Cell–Cell Interactions (ranked by probability)", nrow(top_tbl)),
    x = "Communication probability",
    y = NULL,
    color = "Pathway"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(outfile, p_single,
       width = 15, height = max(6, 0.38 * nrow(top_tbl)), dpi = 300, bg = "white")

cat("\n[Single summary] Saved:", outfile, "\n")

#4. Bubbleplot for SPP1, C1Q and related pathways.
# ---- Single bubble: supervisor-focused pathways (SPP1, C1Q, related ECM) ----
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})

#table of all interactions
comm_df <- subsetCommunication(cellchat)

#fuzzy match pathways of interest
interest_patterns <- c("SPP1", "C1Q", "FN1", "LAMININ", "COLLAGEN", "^ITG[AB]")  # add/remove as needed
avail_paths <- sort(unique(comm_df$pathway_name))
matched_paths <- avail_paths[Reduce(`|`, lapply(interest_patterns, function(p) grepl(p, avail_paths, ignore.case = TRUE)))]

if (length(matched_paths) == 0) {
  stop("None of the requested pathways (SPP1/C1Q/ECM-related) are present in your CellChat results.")
}

cat("[Matched pathways]:", paste(matched_paths, collapse = ", "), "\n")

#restrict to those pathways and pick top senders/targets for readability
subset_comm <- comm_df %>% filter(pathway_name %in% matched_paths)

top_n_sources <- 10   # tweak
top_n_targets <- 10

rank_tbl <- subset_comm %>%
  group_by(source, target) %>%
  summarise(score = sum(prob), .groups = "drop")

top_sources <- rank_tbl %>%
  group_by(source) %>% summarise(tot = sum(score), .groups = "drop") %>%
  arrange(desc(tot)) %>% slice_head(n = top_n_sources) %>% pull(source)

top_targets <- rank_tbl %>%
  group_by(target) %>% summarise(tot = sum(score), .groups = "drop") %>%
  arrange(desc(tot)) %>% slice_head(n = top_n_targets) %>% pull(target)

# guard: ensure something remains
if (nrow(subset_comm %>% filter(source %in% top_sources, target %in% top_targets)) == 0) {
  stop("After restricting to top sources/targets, no interactions remain. Increase top_n_*.")
}

#single combined bubble
p_focus <- netVisual_bubble(
  cellchat,
  signaling      = matched_paths,       # <- we use signaling here
  sources.use    = top_sources,
  targets.use    = top_targets,
  remove.isolate = TRUE
) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  ggtitle("Bubbleplot")

ggsave(
  filename = file.path(output_dir, "Bubbleplot.png"),
  plot = p_focus,
  width = 16, height = 7, dpi = 300, bg = "white"
)

#quick console summary
cat("\n[Focus bubble saved]\nPathways:", paste(matched_paths, collapse = ", "),
    "\nSources:", paste(top_sources, collapse = ", "),
    "\nTargets:", paste(top_targets, collapse = ", "), "\n")


# ---- 5. Chord diagrams for selected pathways ----
chord_sig <- c("MIF","GALECTIN","CD99","LCK", "SPP1")
for (sig in chord_sig) {
  png(file.path(output_dir, sprintf("Chord_%s.png", sig)), width = 4000, height = 3000, res = 300)
  netVisual_aggregate(cellchat, signaling = sig, layout = "chord")
  dev.off()
}
# GALECTIN small variant (kept from original)
png(file.path(output_dir, "Chord_GALECTIN_small.png"), width = 8, height = 8, units = "in", res = 300)
netVisual_aggregate(cellchat, signaling = "GALECTIN", layout = "chord")
dev.off()

# ---- 6. CES1 averages by Sub_Cluster + colors ----
ces1_avg_sub <- AverageExpression(seurat, features = "CES1", group.by = "Sub_Cluster")$scRNA
ces1_vec <- as.numeric(ces1_avg_sub[1, ])
names(ces1_vec) <- colnames(ces1_avg_sub)
# scale 0-1 for palette
ces1_min <- min(ces1_vec, na.rm = TRUE)
ces1_max <- max(ces1_vec, na.rm = TRUE)
ces1_norm <- if (isTRUE(all.equal(ces1_min, ces1_max))) rep(0, length(ces1_vec)) else (ces1_vec - ces1_min) / (ces1_max - ces1_min)
ces1_colors <- scales::col_numeric(palette = c("#4575b4", "#ffdd57", "#d73027"), domain = NULL)(as.vector(ces1_norm))
names(ces1_colors) <- names(ces1_vec)

# ---- 7. Annotated chord diagrams with CES1 node colors ----
suppressPackageStartupMessages({library(scales)})

#Diagnostics to ensure we’re aligned
clusters <- names(table(cellchat@idents))              # CellChat vertex order
cat("[CellChat clusters]:", length(clusters), "\n")

#Compute CES1 mean EXACTLY on CellChat groups (avoid name mismatches)
ces1_vec <- tapply(
  X = FetchData(seurat, vars = "CES1")[,1],
  INDEX = cellchat@idents,                            # same grouping CellChat uses
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Ensure order matches vertex order and fill missing with 0
ces1_vec <- ces1_vec[clusters]
ces1_vec[is.na(ces1_vec)] <- 0

#Normalize and create color vector (same length as clusters)
rng <- range(ces1_vec, na.rm = TRUE)
ces1_norm <- if (diff(rng) == 0) rep(0, length(ces1_vec)) else (ces1_vec - rng[1]) / diff(rng)
node_cols <- scales::col_numeric(c("#4575b4", "#ffdd57", "#d73027"), domain = NULL)(ces1_norm)

# Quick check
cat("[Check] length(color.use) =", length(node_cols), " length(clusters) =", length(clusters), "\n")

#Helper: plot only if pathway has edges, use color.use (your CellChat version supports it)
plot_chord_annotated <- function(sig, out_png, thresh = 0) {
  if (nrow(subsetCommunication(cellchat, signaling = sig)) == 0L) {
    message(sprintf("[Skip] No edges for '%s'.", sig)); return(invisible(FALSE))
  }
  png(out_png, width = 4000, height = 3000, res = 300)
  CellChat::netVisual_aggregate(
    object        = cellchat,
    signaling     = sig,
    layout        = "chord",
    color.use     = node_cols,          # MUST be in the same order as 'clusters'
    thresh        = thresh,
    vertex.label.cex = 0.8
  )
  dev.off()
  TRUE
}

#Run for your pathways (adjust list as needed)
for (sig in c("MIF","GALECTIN","CD99","LCK", "SPP1")) {
  out <- file.path(output_dir, sprintf("%s_Annotated.png", sig))
  ok  <- plot_chord_annotated(sig, out_png = out, thresh = 0)
  if (ok) cat("[Saved]", out, "\n")
}



# ---- 8. CES1 overlay on t-SNE reduction ----
# why: names differ across pipelines; normalize to gTSNE_1/2 if present else fall back
fetch_tsne <- function(obj) {
  md <- obj@meta.data
  if (all(c("gTSNE_1","gTSNE_2") %in% colnames(md))) {
    return(md[, c("gTSNE_1","gTSNE_2")])
  }
  if ("tsne" %in% names(Embeddings(obj))) {
    df <- Embeddings(obj, "tsne") %>% as.data.frame()
  } else if ("tsne" %in% tolower(names(seurat@reductions))) {
    df <- Embeddings(obj, names(obj@reductions)[grepl("tsne", names(obj@reductions), ignore.case = TRUE)][1]) %>% as.data.frame()
  } else {
    stop("No t-SNE/gTSNE coordinates found.")
  }
  colnames(df)[1:2] <- c("gTSNE_1","gTSNE_2")
  df
}

ces1_data <- tryCatch({
  tsne_df <- fetch_tsne(seurat)
  x <- FetchData(seurat, vars = c("CES1"))
  cbind(tsne_df, CES1 = x$CES1)
}, error = function(e) {
  # fallback: use PCA 1/2 labels to still generate a plot
  df <- FetchData(seurat, vars = c("CES1","PC_1","PC_2"))
  names(df) <- c("CES1","gTSNE_1","gTSNE_2")
  df
})

ces1_data$bin <- cut(ces1_data$CES1, breaks = c(0, 0.5, 1, 1.5, 2, Inf),
                     labels = c("0–0.5","0.5–1","1–1.5","1.5–2",">2"))
ces1_data$bin <- as.character(ces1_data$bin)
ces1_data$bin[is.na(ces1_data$bin)] <- "NA"
ces1_data$bin <- factor(ces1_data$bin, levels = c("0–0.5","0.5–1","1–1.5","1.5–2",">2","NA"))

cols_tsne <- c("0–0.5"="#20908d","0.5–1"="#feb24c","1–1.5"="#3a528b","1.5–2"="#bd0026",">2"="#005824","NA"="#d9d9d9")

p_tsne <- ggplot() +
  geom_point(data = ces1_data, aes(x = gTSNE_1, y = gTSNE_2), color = "#d9d9d9", size = 0.4) +
  geom_point(data = ces1_data[ces1_data$bin != "NA",], aes(x = gTSNE_1, y = gTSNE_2, color = bin), size = 0.7) +
  scale_color_manual(values = cols_tsne) +
  ggtitle("t-SNE - CES1 Expression") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(output_dir, "CES1_Expression_Across_Cell_Types.png"), p_tsne, width = 8, height = 6, dpi = 300, bg = "white")

# ---- 9. PCA with CES1 ----
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, features = VariableFeatures(seurat))
seurat <- RunPCA(seurat, features = VariableFeatures(seurat))

ces1_pca <- FetchData(seurat, vars = c("CES1","PC_1","PC_2"))
ces1_pca$bin <- cut(ces1_pca$CES1, breaks = c(0, 0.5, 1, 1.5, 2, Inf),
                    labels = c("0–0.5","0.5–1","1–1.5","1.5–2",">2"))
ces1_pca$bin <- as.character(ces1_pca$bin)
ces1_pca$bin[is.na(ces1_pca$bin)] <- "NA"
ces1_pca$bin <- factor(ces1_pca$bin, levels = c("0–0.5","0.5–1","1–1.5","1.5–2",">2","NA"))
cols_pca <- c("0–0.5"="#4292c6","0.5–1"="#feb24c","1–1.5"="#f03b20","1.5–2"="#bd0026",">2"="#005824","NA"="#d9d9d9")

p_pca <- ggplot() +
  geom_point(data = ces1_pca, aes(x = PC_1, y = PC_2), color = "#d9d9d9", size = 0.5) +
  geom_point(data = ces1_pca[ces1_pca$bin != "NA",], aes(x = PC_1, y = PC_2, color = bin), size = 0.7) +
  scale_color_manual(values = cols_pca) +
  ggtitle("PCA - CES1 Expression") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(output_dir, "PCA_CES1_Expression.png"), p_pca, width = 8, height = 6, dpi = 300, bg = "white")

# ---- 10. UMAP with CES1 ----
if (!"umap" %in% names(Embeddings(seurat))) {
  seurat <- RunUMAP(seurat, dims = 1:20)
}

# Extract UMAP coordinates
umap_coords <- Embeddings(seurat, "umap") %>% as.data.frame()
colnames(umap_coords)[1:2] <- c("UMAP_1", "UMAP_2")

# Get CES1 expression and bin it
ces1_expr <- FetchData(seurat, vars = "CES1")
ces1_expr$CES1_bin <- cut(
  ces1_expr$CES1,
  breaks = c(0, 0.5, 1, 1.5, 2, Inf),
  labels = c("0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2"),
  include.lowest = TRUE
)
ces1_expr$CES1_bin <- as.character(ces1_expr$CES1_bin)
ces1_expr$CES1_bin[is.na(ces1_expr$CES1_bin)] <- "NA"

# Combine into one data frame for plotting
umap_df <- cbind(
  umap_coords,
  CES1       = ces1_expr$CES1,
  CES1_Level = factor(ces1_expr$CES1_bin,
                      levels = c("0–0.5", "0.5–1", "1–1.5", "1.5–2", ">2", "NA"))
)

# Define bins to highlight
highlight_bins <- c("0.5–1", "1–1.5", "1.5–2", ">2")

# Make the plot
p_umap <- ggplot() +
  # Grey background for all cells
  geom_point(data = umap_df, aes(UMAP_1, UMAP_2),
             color = "grey80", size = 0.6) +
  # Overlay colored points for CES1-expressing bins
  geom_point(data = subset(umap_df, CES1_Level %in% highlight_bins),
             aes(UMAP_1, UMAP_2, color = CES1_Level),
             size = 0.8) +
  scale_color_manual(
    values = c(
      "0.5–1"  = "orange",
      "1–1.5"  = "red",
      "1.5–2"  = "darkred",
      ">2"     = "darkgreen"
    ),
    name = "CES1 Level"
  ) +
  labs(title = "UMAP - CES1 Expression") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Save
ggsave(file.path(output_dir, "UMAP - CES1 Expression.png"),
       p_umap, width = 8, height = 6, dpi = 300, bg = "white")


# ---- 11. Tumor vs Normal box/violin ----
stopifnot("Tissue" %in% colnames(seurat@meta.data))
ces1_expr2 <- FetchData(seurat, vars = "CES1")
ces1_expr2$Tissue <- seurat$Tissue

b1 <- ggplot(ces1_expr2, aes(x = Tissue, y = CES1, fill = Tissue)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(shape = 16, size = 1, width = 0.1, aes(color = Tissue)) +
  theme_minimal() +
  labs(title = "CES1 Expression: Tumour vs Normal", y = "Log(TPM+1)", x = "") +
  scale_fill_manual(values = c("T" = "firebrick", "N" = "steelblue")) +
  scale_color_manual(values = c("T" = "firebrick", "N" = "steelblue"))

ggsave(file.path(output_dir, "CES1 Expression - Tumour vs Normal.png"), b1, width = 8, height = 6, dpi = 300, bg = "white")

# ---- 12. Link CES1 to communications (sender/receiver) ----
#Cannot get a clean graph for subcluster


# ---- 13. Top-100 interaction counts per role ----
all_communications <- subsetCommunication(cellchat)

top_100_pathways <- all_communications %>%
  arrange(desc(prob)) %>%
  slice_head(n = 100)

# Count WITHOUT zero-levels
source_counts <- top_100_pathways %>%
  mutate(source = as.character(source)) %>%   # remove factor levels
  count(source, sort = TRUE, name = "count")

target_counts <- top_100_pathways %>%
  mutate(target = as.character(target)) %>%
  count(target, sort = TRUE, name = "count")

# (Optional) keep only top 25
source_counts_plot <- source_counts %>% slice_head(n = 25)
target_counts_plot <- target_counts %>% slice_head(n = 25)

# Sources – single color
p_src <- ggplot(source_counts_plot,
                aes(x = reorder(source, count), y = count)) +
  geom_col(fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Count of Each Subcluster as Source",
       x = "Source Subcluster", y = "Count")

# Targets – single color
p_tgt <- ggplot(target_counts_plot,
                aes(x = reorder(target, count), y = count)) +
  geom_col(fill = "salmon") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Count of Each Subcluster as Target",
       x = "Target Subcluster", y = "Count")

ggsave(file.path(output_dir, "Count_of_Source_Subclusters.png"),
       p_src, width = 7, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Count_of_Target_Subclusters.png"),
       p_tgt, width = 7, height = 8, dpi = 300, bg = "white")


# ---- 14. Save processed CellChat object ----
save(cellchat, file = file.path(output_dir, "cellchat_CRC_CES1_analysis_subcluster.RData"))

cat("\n[Done] Outputs written to:", output_dir, "\n")

#-------------------
#-------------------
# =========================
# CES1: tidy plotting & stats
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(broom.mixed)
})

# ---- Checks ----
stopifnot(exists("seurat"))
stopifnot("Tissue" %in% colnames(seurat@meta.data))

# Pick your sample ID column (edit here if needed)
sample_col <- if ("Sample" %in% names(seurat@meta.data)) "Sample" else
  if ("Patient" %in% names(seurat@meta.data)) "Patient" else NA
if (is.na(sample_col)) stop("Please add a per-sample column to meta.data (e.g., Sample or Patient).")

# Shared palette
pal <- c(N = "steelblue", P = "goldenrod", T = "firebrick")

# =========================
# 1) Per‑tissue cell counts
# =========================
counts_by_tissue <- seurat@meta.data %>%
  count(Tissue, name = "cells") %>%
  arrange(Tissue)

p_counts <- ggplot(counts_by_tissue, aes(Tissue, cells, fill = Tissue)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = cells), vjust = -0.5, size = 4) +
  scale_fill_manual(values = pal) +
  labs(title = "Number of cells per tissue", x = NULL, y = "Cells") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "cells_per_tissue_counts.png"),
       p_counts, width = 6, height = 4.2, dpi = 300, bg = "white")

# ==========================================
# 2) Build tidy table for per‑sample analyses
# ==========================================
dat <- data.frame(
  Tissue   = seurat$Tissue,
  SampleID = seurat@meta.data[[sample_col]],
  CES1     = FetchData(seurat, "CES1")[, 1]
)

# Per‑sample summaries
ces1_by_sample <- dat %>%
  group_by(SampleID, Tissue) %>%
  summarise(
    n_cells     = n(),
    mean_CES1   = mean(CES1, na.rm = TRUE),
    median_CES1 = median(CES1, na.rm = TRUE),
    .groups = "drop"
  )

# ==========================================
# 2b) CES1 expression per tissue (all cells)
# (this will be paired with counts in one panel)
# ==========================================
ces1_all <- FetchData(seurat, vars = c("CES1", "Tissue"))

p_expr_tissue <- ggplot(ces1_all, aes(Tissue, CES1, fill = Tissue)) +
  geom_violin(width = 0.9, alpha = 0.75, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(aes(color = Tissue), width = 0.15, size = 0.3, alpha = 0.18, show.legend = FALSE) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(title = "CES1 expression per tissue after normalization", x = NULL, y = "CES1 (log(TPM+1))") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# Save the side‑by‑side panel: counts | expression
panel_counts_expr <- p_counts | p_expr_tissue
ggsave(file.path(output_dir, "panel_counts_and_CES1_expression.png"),
       panel_counts_expr, width = 12, height = 5, dpi = 300, bg = "white")

# =========================
# NEW: CES1+ counts & % + Filtered expression
# =========================
ces1_threshold <- 0.1  # log(TPM+1) threshold

# Get CES1 expression + Tissue
ces1_all <- FetchData(seurat, vars = c("CES1", "Tissue"))

# Summarise counts & percentages
ces1_summary <- ces1_all %>%
  group_by(Tissue) %>%
  summarise(
    total_cells = n(),
    ces1_pos    = sum(CES1 > ces1_threshold, na.rm = TRUE),
    perc_pos    = 100 * ces1_pos / total_cells,
    .groups     = "drop"
  )

# LEFT plot: counts + % label
p_ces1_counts_pct <- ggplot(ces1_summary, aes(x = Tissue, y = ces1_pos, fill = Tissue)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", ces1_pos, perc_pos)),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = pal) +
  labs(title = sprintf("CES1+ cells per tissue (> %.2f)", ces1_threshold),
       x = NULL, y = "CES1-positive cells (count)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# RIGHT plot: filtered violin
ces1_filtered <- ces1_all %>%
  filter(CES1 > ces1_threshold)

# Save removal stats for reporting
filter_stats <- ces1_all %>%
  group_by(Tissue) %>%
  summarise(
    removed       = sum(CES1 <= ces1_threshold, na.rm = TRUE),
    kept          = sum(CES1 > ces1_threshold, na.rm = TRUE),
    total         = n(),
    perc_removed  = 100 * removed / total,
    .groups       = "drop"
  )
write.csv(filter_stats, file.path(output_dir, "CES1_filtering_stats.csv"), row.names = FALSE)

p_ces1_expr_filtered <- ggplot(ces1_filtered, aes(Tissue, CES1, fill = Tissue)) +
  geom_violin(width = 0.9, alpha = 0.75, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(aes(color = Tissue), width = 0.15, size = 0.3, alpha = 0.18, show.legend = FALSE) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(title = sprintf("CES1 expression per tissue (> %.2f)", ces1_threshold),
       x = NULL, y = "CES1 (log(TPM+1))") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# Save new panel
panel_counts_pct_expr <- p_ces1_counts_pct | p_ces1_expr_filtered
ggsave(file.path(output_dir, "panel_CES1_counts_pct_and_expression_filtered.png"),
       panel_counts_pct_expr, width = 12, height = 5, dpi = 300, bg = "white")


# ====================================
# 3) Per‑sample mean CES1 (separate)
# ====================================
p_mean <- ggplot(ces1_by_sample, aes(Tissue, mean_CES1, color = Tissue)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85) +
  stat_summary(fun = mean, geom = "point", size = 4, color = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.28, color = "black") +
  scale_color_manual(values = pal) +
  labs(title = "Per‑sample mean CES1 by tissue",
       x = NULL, y = "Mean CES1 (log(TPM+1))") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "per_sample_mean_CES1_by_tissue.png"),
       p_mean, width = 6.5, height = 4.6, dpi = 300, bg = "white")

# (Removed the previous combined panel of counts | per‑sample mean)

# =====================================
# 4) Counts per sample (stacked by tissue)
# =====================================
counts_by_sample <- dat %>%
  count(SampleID, Tissue, name = "cells")

p_counts_by_sample <- ggplot(counts_by_sample,
                             aes(x = reorder(SampleID, cells), y = cells, fill = Tissue)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(title = "Cell counts per sample (stacked by tissue)",
       x = "Sample", y = "Cells") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "cells_per_sample_and_tissue.png"),
       p_counts_by_sample, width = 8.5, height = 6, dpi = 300, bg = "white")

# =========================
# 5) Statistics to report
# =========================
# ANOVA on per‑sample means
anova_mean <- aov(mean_CES1 ~ Tissue, data = ces1_by_sample)
capture.output(summary(anova_mean),
               file = file.path(output_dir, "ANOVA_per_sample_mean_CES1.txt"))

# Kruskal–Wallis (non‑parametric)
kw_mean <- kruskal.test(mean_CES1 ~ Tissue, data = ces1_by_sample)
writeLines(capture.output(kw_mean),
           con = file.path(output_dir, "Kruskal_per_sample_mean_CES1.txt"))

# Single‑cell LMM (random intercept per sample)
lmm <- lmer(CES1 ~ Tissue + (1|SampleID), data = dat)
tidy_lmm <- broom.mixed::tidy(lmm, effects = "fixed", conf.int = TRUE)
write.csv(tidy_lmm, file.path(output_dir, "LMM_single_cell_CES1_by_tissue.csv"),
          row.names = FALSE)


############################
#### Machine Learning Part #
############################

suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)   # instead of tidyr::pivot_wider
  library(randomForest)
  library(ggplot2)
})

output_dir <- "E:/Desktop/Dissertation/Data/Colon_zhang/output/CellChat_subcluster_ML"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 1. Extract communication data (only TAM senders)
# ----------------------------------------------------------
comm_df <- subsetCommunication(cellchat)

tam_comm <- comm_df %>%
  filter(source %in% c("hM13_TAM-SPP1", "hM12_TAM-C1QC")) %>%
  mutate(group = ifelse(source == "hM13_TAM-SPP1", "SPP1", "C1QC"))

# ----------------------------------------------------------
# 2. Create interaction identifier including RECEIVER
# ----------------------------------------------------------
tam_comm$feature <- paste(tam_comm$ligand, 
                          tam_comm$receptor, 
                          tam_comm$pathway_name, 
                          tam_comm$target,
                          sep = "_")

# ----------------------------------------------------------
# 3. Find unique vs common interactions
# ----------------------------------------------------------
spp1_feats <- unique(tam_comm$feature[tam_comm$group == "SPP1"])
c1qc_feats <- unique(tam_comm$feature[tam_comm$group == "C1QC"])

common_feats <- intersect(spp1_feats, c1qc_feats)
spp1_unique <- setdiff(spp1_feats, common_feats)
c1qc_unique <- setdiff(c1qc_feats, common_feats)

# Save lists with labels
write.csv(data.frame(Interaction = common_feats, Group = "Common"),
          file.path(output_dir, "Common_interactions.csv"), row.names = FALSE)

write.csv(data.frame(Interaction = spp1_unique, Group = "SPP1"),
          file.path(output_dir, "SPP1_unique_interactions.csv"), row.names = FALSE)

write.csv(data.frame(Interaction = c1qc_unique, Group = "C1QC"),
          file.path(output_dir, "C1QC_unique_interactions.csv"), row.names = FALSE)

# ----------------------------------------------------------
# 4. Encode into ML matrix using dcast
# ----------------------------------------------------------
ml_encoded <- dcast(tam_comm, group ~ feature, 
                    value.var = "prob", fun.aggregate = sum, fill = 0)
ml_encoded$group <- factor(ml_encoded$group)

# Sanitize names for RF
colnames(ml_encoded) <- make.names(colnames(ml_encoded))

# Keep mapping of safe names -> original names
name_map <- data.frame(
  Safe = make.names(c(spp1_unique, c1qc_unique, common_feats)),
  Original = c(spp1_unique, c1qc_unique, common_feats),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 5. Random Forest Classification
# ----------------------------------------------------------
set.seed(123)
rf_model <- randomForest(group ~ ., data = ml_encoded, importance = TRUE, ntree = 500)
print(rf_model)

# Feature importance
imp <- importance(rf_model)
imp_df <- data.frame(Feature = rownames(imp), imp)

# Map back to original names
imp_df <- imp_df %>%
  left_join(name_map, by = c("Feature" = "Safe")) %>%
  mutate(DisplayName = ifelse(!is.na(Original), Original, Feature))

# Order by importance
imp_df <- imp_df[order(imp_df$MeanDecreaseGini, decreasing = TRUE), ]

# Save full importance
write.csv(imp_df, file.path(output_dir, "RF_feature_importance.csv"), row.names = FALSE)

# ----------------------------------------------------------
# 6. Plot UNIQUE interactions only
# ----------------------------------------------------------
unique_feats <- c(spp1_unique, c1qc_unique)

imp_unique <- imp_df %>%
  filter(DisplayName %in% unique_feats) %>%
  mutate(Group = ifelse(DisplayName %in% spp1_unique, "SPP1", "C1QC"))

p <- ggplot(imp_unique, aes(x = reorder(DisplayName, MeanDecreaseGini), 
                            y = MeanDecreaseGini, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Unique interactions distinguishing SPP1 vs C1QC",
       x = "Interaction (ligand_receptor_pathway_target)",
       y = "Importance (MeanDecreaseGini)") +
  theme_bw() +
  scale_fill_manual(values = c("SPP1" = "darkorange", "C1QC" = "steelblue"))

ggsave(filename = file.path(output_dir, "RF_unique_feature_importance.png"),
       plot = p,
       width = 12, height = 15, dpi = 300)

# ----------------------------------------------------------
# Plot 2: Only nonzero-importance unique interactions
# ----------------------------------------------------------
imp_unique_nonzero <- imp_unique %>%
  filter(MeanDecreaseGini > 0)

p_nonzero <- ggplot(imp_unique_nonzero, aes(x = reorder(DisplayName, MeanDecreaseGini),
                                            y = MeanDecreaseGini, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Unique interactions (nonzero importance only)",
       x = "Interaction (ligand_receptor_pathway_target)",
       y = "Importance (MeanDecreaseGini)") +
  theme_bw() +
  scale_fill_manual(values = c("SPP1" = "darkorange", "C1QC" = "steelblue"))

ggsave(filename = file.path(output_dir, "RF_unique_feature_importance_nonzero.png"),
       plot = p_nonzero,
       width = 12, height = 15, dpi = 300)

# Also save nonzero table
write.csv(imp_unique_nonzero, 
          file.path(output_dir, "RF_unique_feature_importance_nonzero.csv"), 
          row.names = FALSE)


# Also save unique feature importance as CSV
write.csv(imp_unique, file.path(output_dir, "RF_unique_feature_importance.csv"), row.names = FALSE)

# ----------------------------------------------------------
# 7. Summarize ML-important unique interactions by RECEIVER
# ----------------------------------------------------------

# Function to extract the target (receiver) from interaction string
get_target <- function(x) {
  sapply(strsplit(x, "_"), function(parts) tail(parts, 1))
}

# Keep only important unique features (nonzero importance)
imp_unique_nonzero <- imp_unique %>%
  filter(MeanDecreaseGini > 0) %>%
  mutate(Target = get_target(DisplayName))

# Count important unique interactions per group × target
receiver_summary_full <- imp_unique_nonzero %>%
  group_by(Group, Target) %>%
  summarise(UniqueInteractions = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = UniqueInteractions,
    values_fill = 0
  ) %>%
  mutate(Difference = C1QC - SPP1)

# Save full table
write.csv(receiver_summary_full,
          file.path(output_dir, "Receiver_summary_SPP1_vs_C1QC_ML.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 8. Diverging barplot of difference (C1QC − SPP1)
# ----------------------------------------------------------
p_diff <- ggplot(receiver_summary_full,
                 aes(x = reorder(Target, Difference), y = Difference,
                     fill = Difference > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "darkorange"),
                    labels = c("FALSE" = "SPP1-dominant", "TRUE" = "C1QC-dominant")) +
  labs(title = "Receiver bias (C1QC vs SPP1, ML-driven unique interactions)",
       x = "Receiver cell type",
       y = "Difference (C1QC − SPP1)") +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(filename = file.path(output_dir, "Receiver_bias_C1QC_vs_SPP1.png"),
       plot = p_diff, width = 10, height = 6, dpi = 300)

# ----------------------------------------------------------
# 9. Plot receiver summary (ML-driven counts of unique interactions)
# ----------------------------------------------------------

# Add Target column if missing
imp_unique_nonzero <- imp_unique_nonzero %>%
  mutate(Target = get_target(DisplayName))

# Count ML-driven unique interactions per receiver and group
receiver_summary_long <- imp_unique_nonzero %>%
  group_by(Group, Target) %>%
  summarise(UniqueInteractions = n(), .groups = "drop")

# Plot ML-influenced counts per receiver
p_receiver_ml <- ggplot(receiver_summary_long,
                        aes(x = Target, y = UniqueInteractions, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Unique interactions by receiver cell type (ML-driven)",
       x = "Receiver cell type",
       y = "Number of unique interactions") +
  theme_bw() +
  scale_fill_manual(values = c("SPP1" = "darkorange", "C1QC" = "steelblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save
ggsave(file.path(output_dir, "Unique_interactions_by_receiver_ML.png"),
       p_receiver_ml, width = 10, height = 6, dpi = 300)


#to check the unique interactions count in the file and the plot
library(dplyr)
library(readr)
library(stringr)

# -------------------------------
# Load your unique interaction CSVs
# -------------------------------
c1qc <- read_csv("E:/Desktop/Dissertation/Data/Colon_zhang/output/CellChat_subcluster_ML/C1QC_unique_interactions.csv")
spp1 <- read_csv("E:/Desktop/Dissertation/Data/Colon_zhang/output/CellChat_subcluster_ML/SPP1_unique_interactions.csv")

# -------------------------------
# Helper: extract receiver cell type
# -------------------------------
get_receiver <- function(x) {
  sapply(strsplit(x, "_"), function(parts) tail(parts, 1))
}

# -------------------------------
# Count unique interactions per receiver
# -------------------------------
c1qc_counts <- c1qc %>%
  mutate(Receiver = get_receiver(Interaction)) %>%
  count(Receiver, name = "C1QC_count")

spp1_counts <- spp1 %>%
  mutate(Receiver = get_receiver(Interaction)) %>%
  count(Receiver, name = "SPP1_count")

# -------------------------------
# Merge both counts
# -------------------------------
merged_counts <- full_join(c1qc_counts, spp1_counts, by = "Receiver") %>%
  replace_na(list(C1QC_count = 0, SPP1_count = 0))

# -------------------------------
# OPTIONAL: Cross-check with the data used for plotting
# (this is the file generated in step 9 of your ML code)
# -------------------------------
receiver_summary <- read_csv("E:/Desktop/Dissertation/Data/Colon_zhang/output/CellChat_subcluster_ML/Receiver_summary_SPP1_vs_C1QC_ML.csv")
merged_counts <- merged_counts %>%
left_join(receiver_summary, by = c("Receiver" = "Target"))

# -------------------------------
# Save or view
# -------------------------------
print(merged_counts)
write_csv(merged_counts, "E:/Desktop/Dissertation/Data/Colon_zhang/output/Crosscheck_unique_interactions_per_receiver.csv")

# -------------------------
# Split into 3 equal chunks (Option 1)
# -------------------------
library(dplyr)
library(ggplot2)
library(readr)

# Load feature importance file
# Now when you read the importance file, it will work
imp_unique <- read_csv(file.path(output_dir, "RF_unique_feature_importance.csv"))
# Keep only interactions with non-zero importance
imp_unique <- imp_unique %>%
  filter(MeanDecreaseGini > 0)

# Split into 3 equal parts
n <- nrow(imp_unique)
chunks <- split(imp_unique, cut(seq_len(n), 3, labels = FALSE))

# Function to plot and save each chunk
plot_chunk <- function(df, index) {
  p <- ggplot(df, aes(x = reorder(DisplayName, MeanDecreaseGini),
                      y = MeanDecreaseGini, fill = Group)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Unique interactions (Part", index, ")"),
         x = "Interaction",
         y = "Importance (MeanDecreaseGini)") +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = c("SPP1" = "darkorange", "C1QC" = "steelblue"))
  
  ggsave(filename = file.path(output_dir, paste0("RF_unique_feature_importance_part", index, ".png")),
         plot = p, width = 10, height = 12, dpi = 300)
}

# Loop through each chunk and save
for (i in seq_along(chunks)) {
  plot_chunk(chunks[[i]], i)
}

