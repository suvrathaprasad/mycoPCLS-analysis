#!/usr/bin/env Rscript

# ============================================================================
# modify_plots.R
# Rebuild/modify publication plots using only the .rds files created by
# run_plots.R. Edit the PLOT SETTINGS section only.
# Usage:
# Sheep: Rscript modify_plots.R results/sheep/plots/plots_data sheep results/sheep/plots/modified_plots
# Cattle: Rscript modify_plots.R results/cattle/plots/plots_data cattle results/cattle/plots/modified_plots
# ============================================================================

# ---- PLOT SETTINGS ----------------------------------------------------------

SPECIES_DISPLAY <- list(cattle = "Bovine PCLS transcriptomes", sheep = "Sheep PCLS transcriptomes")
CONDITION_DISPLAY <- list(
  cattle = c(NK_mP_vs_NK_oP = "NC w/ PANTA vs. NC w/o PANTA",
             oP_vs_NK_oP = "M. bovis Infection w/o PANTA vs. NC w/o PANTA",
             mP_vs_NK_mP = "M. bovis Infection w/ PANTA vs. NC w/ PANTA",
             mP_vs_oP = "M. bovis Infection w/ PANTA vs. M. bovis Infection w/o PANTA"),
  sheep = c(infection_vs_NK = "M. bovis Infection vs. NC"))

PLOT_WIDTH <- 8
PLOT_HEIGHT <- 6
BASE_SIZE <- 11
THEME <- ggplot2::theme_bw(base_size = BASE_SIZE) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))

VOLCANO_LABEL_N <- 10
HEATMAP_WIDTH <- 11
HEATMAP_HEIGHT <- 6
EMAP_WIDTH <- 15
EMAP_HEIGHT <- 15

# ---- ARGUMENTS --------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: modify_plots.R <plots_data_dir> <species> <outdir>")

plots_data_dir <- args[1]
species <- tolower(args[2])
outdir <- args[3]

if (!dir.exists(plots_data_dir)) stop("plots_data directory not found: ", plots_data_dir)
if (!species %in% names(SPECIES_DISPLAY)) stop("species must be 'cattle' or 'sheep'")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(clusterProfiler)
  library(enrichplot)
  library(grid)
  if (species == "cattle") library(ggVennDiagram)
})

species_label <- SPECIES_DISPLAY[[species]]

label_comparison <- function(x) {
  map <- CONDITION_DISPLAY[[species]]
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

# ============================================================================
# 1. PCA
# ============================================================================

f <- file.path(plots_data_dir, paste0("1_pca_plot_data_", species, ".rds"))

if (file.exists(f)) {
  x <- readRDS(f)
  p <- ggplot(x$pca_df, aes(PC1, PC2, color = group, shape = dpi)) +
    geom_point(size = 3.5, alpha = 0.9) +
    labs(title = species_label,
         x = paste0("PC1 (", x$pct_var[1], "%)"),
         y = paste0("PC2 (", x$pct_var[2], "%)"),
         color = "Group", shape = "dpi") +
    THEME
  ggsave(file.path(outdir, paste0("1_pca_", species, ".pdf")), p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
}

# ============================================================================
# 2. HEATMAPS
# ============================================================================

files <- list.files(plots_data_dir, "^2_heatmap_data_.*\\.rds$", full.names = TRUE)

for (f in files) {
  x <- readRDS(f)
  tag <- sub("^2_heatmap_data_", "", tools::file_path_sans_ext(basename(f)))
  cd <- sub(paste0("^", species, "_"), "", tag)

  pdf(file.path(outdir, paste0("2_heatmap_", species, "_", cd, ".pdf")), width = HEATMAP_WIDTH, height = HEATMAP_HEIGHT)
  pheatmap(x$mat_z, annotation_col = x$annotation_col, main = paste(species_label, label_comparison(cd), "top candidate genes"),
           show_rownames = TRUE, angle_col = "45", fontsize_row = 7, fontsize_col = 5,
           fontsize = 7, annotation_names_col = TRUE, annotation_legend = TRUE)
  dev.off()
}

# ============================================================================
# 3. VOLCANOES
# ============================================================================

files <- list.files(plots_data_dir, "^3_volcano_data_.*\\.rds$", full.names = TRUE)

for (f in files) {
  x <- readRDS(f)
  df <- if (is.list(x) && "df" %in% names(x)) x$df else x
  labels <- if (is.list(x) && "labels" %in% names(x)) x$labels else df[df$candidate, , drop = FALSE]
  comparison <- if (is.list(x) && "comparison" %in% names(x)) x$comparison else "unknown"
  dpi <- if (is.list(x) && "dpi" %in% names(x)) x$dpi else "unknown"

  labels$label <- if ("label" %in% names(labels)) labels$label else labels$gene_id

  p <- ggplot(df, aes(log2FoldChange, -log10(padj_plot), color = direction)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "black")) +
    ggrepel::geom_text_repel(data = labels, aes(label = label), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    labs(title = paste(species_label, label_comparison(comparison), "dpi", dpi),
         x = "log2 fold change", y = "-log10(padj)", color = NULL) +
    THEME

  tag <- tools::file_path_sans_ext(basename(f))
  tag <- sub("^3_volcano_data_", "", tag)
  ggsave(file.path(outdir, paste0("3_volcano_", tag, ".pdf")), p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
}

# ============================================================================
# 4. TEMPORAL TRAJECTORIES
# ============================================================================

files <- list.files(plots_data_dir, "^4_trajectory_data_.*\\.rds$", full.names = TRUE)

for (f in files) {
  x <- readRDS(f)
  traj <- if (is.list(x) && "traj" %in% names(x)) x$traj else x
  comparison <- if (is.list(x) && "comparison" %in% names(x)) x$comparison else "unknown"

  p <- ggplot(traj, aes(dpi, log2FC, color = gene_label, group = gene_id)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = paste(species_label, label_comparison(comparison), "top candidate trajectories"),
         x = "dpi", y = "log2 fold change", color = "Gene") +
    THEME

  ggsave(file.path(outdir, paste0("4_trajectory_", species, "_", comparison, ".pdf")), p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
}

# ============================================================================
# 5. VENN
# ============================================================================

f <- file.path(plots_data_dir, paste0("5_venn_data_", species, ".rds"))

if (file.exists(f) && species == "cattle") {
  x <- readRDS(f)
  p <- ggVennDiagram(x$candidate_sets, category.names = x$labels, set_size = 4, label_size = 4) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    labs(title = paste(species_label, "candidate gene overlap across comparisons")) +
    THEME +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(file.path(outdir, paste0("5_candidate_overlap_", species, ".pdf")), p, width = 12, height = 9)
}

# ============================================================================
# 6. EMAP
# ============================================================================

files <- list.files(plots_data_dir, "^6_emap_data_.*\\.rds$", full.names = TRUE)

for (f in files) {
  x <- readRDS(f)
  res <- if (is.list(x) && "result" %in% names(x)) x$result else x
  tag <- if (is.list(x) && "tag" %in% names(x)) x$tag else sub("^6_emap_data_|\\.rds$", "", basename(f))
  src <- if (is.list(x) && "source" %in% names(x)) x$source else ""

  p <- tryCatch(
    emapplot(res, showCategory = 20) +
      ggtitle(paste(species_label, tag, src)) +
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 13, face = "bold"),
            legend.key.size = unit(1.3, "cm")),
    error = function(e) NULL)

  if (!is.null(p)) ggsave(file.path(outdir, paste0("6_emap_", tag, "_", src, ".pdf")), p, width = EMAP_WIDTH, height = EMAP_HEIGHT)
}

# ============================================================================
# 7. CROSS-DPI GO
# ============================================================================

files <- list.files(plots_data_dir, "^7_dpi_comparison_data_.*\\.rds$", full.names = TRUE)

for (f in files) {
  x <- readRDS(f)
  res <- if (is.list(x) && "result" %in% names(x)) x$result else x
  comparison <- if (is.list(x) && "comparison" %in% names(x)) x$comparison else "unknown"
  direction <- if (is.list(x) && "direction" %in% names(x)) x$direction else "unknown"

  p <- dotplot(res, showCategory = 5) +
    ggtitle(paste(species_label, label_comparison(comparison), direction, "GO across dpi")) +
    THEME +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(outdir, paste0("7_dpi_comparison_", species, "_", comparison, "_", direction, ".pdf")),
         p, width = 10, height = 8)
}

# ============================================================================
# 8a. CANDIDATES OVER TIME
# ============================================================================

f <- file.path(plots_data_dir, paste0("8a_plot_data_", species, ".rds"))

if (file.exists(f)) {
  cs_long <- readRDS(f)

  p <- ggplot(cs_long, aes(dpi, count, fill = direction)) +
    geom_col(position = "identity") +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
    facet_wrap(~comparison) +
    labs(title = paste(species_label, "candidate genes over time"),
         x = "dpi", y = "Number of candidate genes", fill = NULL) +
    THEME

  ggsave(file.path(outdir, paste0("8a_candidates_overview_", species, ".pdf")),
         p, width = 10, height = 7)
}

# ============================================================================
# 8b. CROSS-COMPARISON TRAJECTORIES
# ============================================================================

f <- file.path(plots_data_dir, paste0("8b_plot_data_", species, ".rds"))

if (file.exists(f)) {
  x <- readRDS(f)
  overlay_long <- x$overlay_long

  p <- ggplot(overlay_long, aes(dpi, log2FC, color = comparison_label, group = comparison)) +
    geom_line(alpha = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~facet_label, scales = "free_y") +
    labs(title = paste(species_label, "gene response across comparisons"),
         x = "dpi", y = "log2 fold change", color = "Comparison") +
    THEME

  ggsave(file.path(outdir, paste0("8b_cross_comparison_trajectories_", species, ".pdf")),
         p, width = 11, height = 8)
}

# ============================================================================
# 8c. RECURRING TERMS
# ============================================================================

f <- file.path(plots_data_dir, paste0("8c_plot_data_", species, ".rds"))

if (file.exists(f)) {
  recurring <- readRDS(f)

  p <- ggplot(recurring, aes(x = reorder(Description, n_significant_results),
                             y = n_significant_results, fill = term_type)) +
    geom_col() +
    coord_flip() +
    labs(title = paste(species_label, "most recurring enriched terms"),
         x = NULL, y = "Number of significant results", fill = "Type") +
    THEME

  ggsave(file.path(outdir, paste0("8c_recurring_terms_", species, ".pdf")),
         p, width = 10, height = 8)
}

# ============================================================================
# 8d. HUB GENES
# ============================================================================

f <- file.path(plots_data_dir, paste0("8d_plot_data_", species, ".rds"))

if (file.exists(f)) {
  hub_genes <- readRDS(f)

  p <- ggplot(hub_genes, aes(x = reorder(label, n_terms), y = n_terms)) +
    geom_col(fill = "darkorange") +
    coord_flip() +
    labs(title = paste(species_label, "candidate genes in the most enriched terms"),
         x = NULL, y = "Number of distinct GO/KEGG terms") +
    THEME

  ggsave(file.path(outdir, paste0("8d_hub_genes_", species, ".pdf")),
         p, width = 9, height = 8)
}

cat("Modified plots written to:", outdir, "\n")

