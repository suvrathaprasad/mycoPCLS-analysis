#!/usr/bin/env Rscript
# ============================================================================
# run_plots.R : publication-ready plots, rebuilt from existing RDS
#
# Reads run_pca_deseq2.R's SECTION_7 dds_*.rds (real DESeqDataSets) and run_enrichment.R's
# GO ORA .rds objects (real enrichResult objects) 
#
# GENE SYMBOLS, NOT ENSEMBL IDs: heatmap rownames, volcano labels, and
# trajectory legends show gene symbols where resolvable, reusing the
# cached annotation map (<species>_annotation_map.rds) rather than doing an
# independent lookup or any new network call. entrez_go is already the
# correct organism's Entrez ID for whichever source_go says (bovine direct,
# or human ortholog) — so this only has to pick the matching OrgDb
# (org.Bt.eg.db for "direct", org.Hs.eg.db for "human_ortholog") and query
# SYMBOL, all local. Genes with no resolvable symbol fall back to their raw
# Ensembl ID & not every gene has a registered symbol.
#
# SECTION 8 (cross-comparison overview plots) READS run_tables.R's OUTPUT —
# TABLES must run before PLOTS. tables_dir is now taken as an explicit 5th CLI
# argument, not derived from deseq2_outdir — run_tables.R writes its *.tsv
# files flat, directly into whatever <outdir> it's given (no "tables/"
# subfolder), so there's no fixed subpath to derive here; the caller must
# point <tables_dir> at wherever those files actually landed. Those tables
# already carry gene_symbol/entrez_id/gene_id from their own join against
# the same annotation cache, so section 8 reuses those columns directly
# rather than re-resolving symbols a second time.
#
# EDITABLE SETTINGS — change labels/titles/colors/counts here, not buried
# in the plotting code below:
# ============================================================================
SPECIES_DISPLAY   <- list(cattle = "Bovine PCLS transcriptomes", sheep = "Sheep PCLS transcriptomes")
COMPARISON_DISPLAY <- list(cattle = "M. bovis Infection ± PANTA Treatment",
                           sheep = "M. bovis Infection vs. Non.Infected")
CONDITION_DISPLAY <- list(cattle = c(NK_mP_vs_NK_oP = "NC w/ PANTA vs. NC w/o PANTA",
                                     oP_vs_NK_oP = "M. bovis Infection w/o PANTA vs. NC w/o PANTA",
                                     mP_vs_NK_mP = "M. bovis Infection  w/ PANTA vs. NC w/ PANTA",
                                     mP_vs_oP = "M. bovis Infection w/ PANTA vs. M. bovis Infection w/o PANTA"), 
                          sheep = c(infection_vs_NK = "M. bovis Infection vs. NC"))
GROUP_DISPLAY <- list(
  sheep = c(NK = "NC", infection = "23MA0307"),
  cattle = c(
    NK_oP = "NC w/o PANTA",
    NK_mP = "NC w/ PANTA",
    oP = "23MA0307 w/o PANTA",
    mP = "23MA0307 w/ PANTA"))

TOP_N_HEATMAP     <- 30     # genes per comparison in the expression heatmap
TOP_N_TRAJECTORY  <- 15     # genes per comparison in the temporal line plot
VOLCANO_LABEL_N   <- 10     # top genes labeled on each volcano plot
TOP_N_RECURRING_TERMS <- 20  # most-recurring enriched terms across all contrasts
TOP_N_HUB_GENES       <- 20  # candidate genes in the most distinct enriched terms
TOP_N_OVERLAY_GENES   <- 8   # genes in the cross-comparison trajectory overlay (cattle only)
TERM_SIG_THRESHOLD    <- 0.05  # p.adjust cutoff for "recurring term" counting
PLOT_WIDTH        <- 8
PLOT_HEIGHT       <- 6
THEME <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
# ============================================================================

# ---- Arguments -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) stop("Usage: run_plots.R <section7_dir> <enrichment_dir> <species> <outdir> <cache_dir> <tables_dir> <pca_rds>")

section7_dir <- args[1]
enrichment_dir <- args[2]
species <- tolower(args[3])
outdir <- args[4]
cache_dir <- args[5]
tables_dir <- args[6]
pca_rds <- args[7]

if (!species %in% names(SPECIES_DISPLAY)) stop("species must be 'cattle' or 'sheep', got: ", species)
if (!dir.exists(section7_dir)) stop("SECTION_7 directory not found: ", section7_dir)
if (!dir.exists(enrichment_dir)) stop("Enrichment directory not found: ", enrichment_dir)
if (!dir.exists(tables_dir)) message("Warning: tables directory not found: ", tables_dir)

species_label <- SPECIES_DISPLAY[[species]]
comparison_label <- COMPARISON_DISPLAY[[species]]

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(clusterProfiler)
  library(enrichplot)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  if (species == "cattle") library(org.Bt.eg.db)
  if (species == "cattle") library(ggVennDiagram)
})

plots_dir <- outdir
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
plots_data_dir <- file.path(plots_dir, "plots_data")
dir.create(plots_data_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Gene ID -> display label lookup ---------------------------------------

annot_map_file <- file.path(cache_dir, paste0(species, "_annotation_map.rds"))
gene_labels <- NULL

if (file.exists(annot_map_file)) {
  annot_map <- readRDS(annot_map_file)

  direct_ids <- unique(na.omit(annot_map$entrez_go[annot_map$source_go == "direct"]))
  ortholog_ids <- unique(na.omit(annot_map$entrez_go[annot_map$source_go == "human_ortholog"]))

  sym_direct <- if (species == "cattle" && length(direct_ids) > 0) AnnotationDbi::select(org.Bt.eg.db, keys = direct_ids, keytype = "ENTREZID", columns = "SYMBOL") else data.frame(ENTREZID = character(0), SYMBOL = character(0))
  sym_ortholog <- if (length(ortholog_ids) > 0) AnnotationDbi::select(org.Hs.eg.db, keys = ortholog_ids, keytype = "ENTREZID", columns = "SYMBOL") else data.frame(ENTREZID = character(0), SYMBOL = character(0))

  sym_lookup <- rbind(sym_direct, sym_ortholog)
  sym_lookup <- sym_lookup[!duplicated(sym_lookup$ENTREZID) & !is.na(sym_lookup$SYMBOL), ]
  symbol_for_entrez <- setNames(sym_lookup$SYMBOL, sym_lookup$ENTREZID)
  resolved <- symbol_for_entrez[as.character(annot_map$entrez_go)]
  gene_labels <- setNames(ifelse(!is.na(resolved), resolved, annot_map$gene_id), annot_map$gene_id)

  cat(sum(!is.na(resolved)), "of", length(gene_labels), "genes resolved to a symbol.\n")
} else {
  message("No annotation cache found at ", annot_map_file, " — plots will show raw Ensembl IDs instead of gene symbols.")
}

label_genes <- function(ids) {
  if (is.null(gene_labels)) return(ids)
  out <- unname(gene_labels[ids])
  out[is.na(out)] <- ids[is.na(out)]
  out
}

label_comparison <- function(x) {
  map <- CONDITION_DISPLAY[[species]]
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

ordered_dpi <- function(x) factor(x, levels = sort(unique(as.numeric(as.character(x)))))

# ---- Discover comparisons --------------------------------------------------

comparison_dirs <- list.dirs(section7_dir, recursive = FALSE, full.names = FALSE)
cat("Found", length(comparison_dirs), "comparison(s) for", species_label, ":", paste(comparison_dirs, collapse = ", "), "\n")

# ============================================================================
# 1. PCA
# ============================================================================

pca_list <- readRDS(pca_rds)
dds_raw <- pca_list$dds_raw

if (is.null(dds_raw)) stop("PCA RDS does not contain dds_raw.")

if (is.null(sizeFactors(dds_raw))) dds_raw <- DESeq2::estimateSizeFactors(dds_raw)

expr_mat <- log2(counts(dds_raw, normalized = TRUE) + 1)
meta <- as.data.frame(colData(dds_raw))

pc <- prcomp(t(expr_mat), scale. = FALSE)
pct_var <- round(100 * (pc$sdev^2 / sum(pc$sdev^2)), 1)

pca_df <- data.frame(
  PC1   = pc$x[, 1], PC2 = pc$x[, 2],
  group = unname(GROUP_DISPLAY[[species]][as.character(meta$condition)]),
  dpi   = ordered_dpi(meta$dpi)
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = dpi)) +
  geom_point(size = 3.5, alpha = 0.9) +
  labs(title = paste(species_label, comparison_label),
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)"),
       color = "Group", shape = "dpi") +
  THEME

saveRDS(list(pca_df = pca_df, pct_var = pct_var), file.path(plots_data_dir, paste0("1_pca_plot_data_", species, ".rds")))
ggsave(file.path(plots_dir, paste0("1_pca_", species, ".pdf")), p_pca, width = PLOT_WIDTH, height = PLOT_HEIGHT)
cat("Wrote PCA plot using dds_raw with", ncol(dds_raw), "individual samples.\n")

# ============================================================================
# 2. Heatmap — top candidate genes
# ============================================================================

for (cd in comparison_dirs) {
  dds_files <- list.files(
    file.path(section7_dir, cd),
    pattern = "^dds_.*\\.rds$",
    full.names = TRUE
  )

  temporal_file <- list.files(
    file.path(section7_dir, cd),
    pattern = "_temporal_summary\\.tsv$",
    full.names = TRUE
  )

  if (length(dds_files) == 0 || length(temporal_file) == 0) {
    message(
      "Skipping heatmap for ", cd,
      " — missing dds RDS or temporal summary."
    )
    next
  }

  dds_pair <- readRDS(dds_files[1])
  temporal <- read.delim(temporal_file[1], stringsAsFactors = FALSE)

  top_genes <- temporal %>%
    dplyr::filter(passed_candidate_filter_any_DPI) %>%
    arrange(desc(max_abs_log2FC)) %>%
    slice_head(n = TOP_N_HEATMAP) %>%
    pull(gene_id)

  if (length(top_genes) < 2) {
    message(
      "Skipping heatmap for ", cd,
      " — fewer than 2 candidate genes."
    )
    next
  }

  norm_counts <- counts(dds_pair, normalized = TRUE)
  top_genes <- intersect(top_genes, rownames(norm_counts))

  if (length(top_genes) < 2) {
    message(
      "Skipping heatmap for ", cd,
      " fewer than 2 candidate genes present in normalized counts."
    )
    next
  }

  mat <- log2(norm_counts[top_genes, , drop = FALSE] + 1)

  mat_z <- t(scale(t(mat)))

  keep <- apply(mat_z, 1, function(x) all(is.finite(x)))

  mat_z <- mat_z[keep, , drop = FALSE]

  if (nrow(mat_z) < 2) {
    message(
      "Skipping heatmap for ", cd,
      " fewer than 2 genes with finite z-scores."
    )
    next
  }

  rownames(mat_z) <- make.unique(label_genes(rownames(mat_z)))

  # Build annotation data frame
  ann_col <- data.frame(
    group = unname(
      GROUP_DISPLAY[[species]][
        as.character(colData(dds_pair)$condition)
      ]
    ),
    dpi = ordered_dpi(colData(dds_pair)$dpi),
    row.names = colnames(dds_pair)
  )

  # Track numerical DPI value for chronological column sorting
  ann_col$dpi_num <- as.numeric(as.character(ann_col$dpi))

  # Order sample columns first by group,
  # then ascending by timepoint
  sample_order <- rownames(ann_col[order(ann_col$group, ann_col$dpi_num), ])

  mat_z <- mat_z[, sample_order, drop = FALSE]
  ann_col <- ann_col[sample_order, , drop = FALSE]

  # Drop helper numeric column
  ann_col$dpi_num <- NULL

  # Drop unused DPI levels and sort levels numerically
  present_dpis_numeric <- sort(unique(as.numeric(as.character(ann_col$dpi))))

  ann_col$dpi <- factor(as.character(ann_col$dpi), levels = as.character(present_dpis_numeric))

  # Drop unused group levels.
  # This ensures the group legend contains only groups
  # actually present in this comparison.
  present_groups <- unique(as.character(ann_col$group))

  ann_col$group <- factor(as.character(ann_col$group), levels = present_groups)

  heatmap_title <- paste(
    strwrap(
      paste(
        species_label,
        label_comparison(cd),
        "top candidate genes"
      ),
      width = 55
    ),
    collapse = "\n"
  )

  pdf(file.path(plots_dir, paste0("2_heatmap_", species, "_", cd, ".pdf")), width = 11, height = PLOT_HEIGHT)

  if (species == "cattle") {

    # Master DPI color lookup
    all_dpi_colors <- c(
      "0" = "#1B9E77",
      "1" = "#D95F02",
      "2" = "#7570B3",
      "3" = "#E7298A",
      "6" = "#66A61E",
      "10" = "#E6AB02"
    )

    # Keep only DPI colors corresponding to DPIs
    # actually present in this heatmap
    active_dpi_colors <- all_dpi_colors[
      names(all_dpi_colors) %in% levels(ann_col$dpi)
    ]

    # Master group color lookup for cattle
    all_group_colors <- c(
      "NC w/o PANTA" = "#4C78A8",
      "NC w/ PANTA" = "#72B7B2",
      "23MA0307 w/o PANTA" = "#F58518",
      "23MA0307 w/ PANTA" = "#E45756"
    )

    # Keep only group colors corresponding to groups
    # actually present in this heatmap
    active_group_colors <- all_group_colors[
      names(all_group_colors) %in% levels(ann_col$group)
    ]

    ann_colors <- list(group = active_group_colors, dpi = active_dpi_colors)

  } else {

    # Master group color lookup for sheep
    all_group_colors <- c(
      "NC" = "#4C78A8",
      "23MA0307" = "#72B7B2"
    )

    # Keep only group colors corresponding to groups
    # actually present in this heatmap
    active_group_colors <- all_group_colors[
      names(all_group_colors) %in% levels(ann_col$group)
    ]

    ann_colors <- list(group = active_group_colors)
  }

  pheatmap(
    mat_z,
    cluster_cols = FALSE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    main = heatmap_title,
    show_rownames = TRUE,
    angle_col = "45",
    fontsize_row = 7,
    fontsize_col = 5,
    fontsize = 7,
    annotation_names_col = TRUE,
    annotation_legend = TRUE
  )

  saveRDS(
    list(
      mat = mat,
      mat_z = mat_z,
      annotation_col = ann_col),
      file.path(plots_data_dir, paste0("2_heatmap_data_", species, "_", cd, ".rds")))

  dev.off()

  cat("Wrote heatmap for", cd, "\n")
}

# ============================================================================
# 3. Volcano plots — per comparison/dpi
# ============================================================================

expl_files <- list.files(section7_dir, pattern = "_DESeq2_exploratory\\.tsv$", recursive = TRUE, full.names = TRUE)

for (ef in expl_files) {
  comparison_dir <- basename(dirname(ef))
  fname <- basename(ef)
  m <- regmatches(fname, regexec("_dpi([^_]+)_DESeq2_exploratory\\.tsv$", fname))
  dpi_label <- if (length(m[[1]]) == 2) m[[1]][2] else "unknown"

  df <- read.delim(ef, stringsAsFactors = FALSE)
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]
  if (nrow(df) == 0) next

  df$direction <- ifelse(df$candidate & df$log2FoldChange > 0, "Up", ifelse(df$candidate & df$log2FoldChange < 0, "Down", "NS"))
  df$padj_plot <- pmax(df$padj, .Machine$double.xmin)

  labels <- df %>% dplyr::filter(candidate) %>% dplyr::arrange(padj) %>% slice_head(n = VOLCANO_LABEL_N)
  labels$label <- label_genes(labels$gene_id)

  p_volc <- ggplot(df, aes(log2FoldChange, -log10(padj_plot), color = direction)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "black")) +
    ggrepel::geom_text_repel(data = labels, aes(label = label), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    labs(title = paste(species_label, label_comparison(comparison_dir), "dpi", dpi_label), x = "log2 fold change", y = "-log10(padj)", color = NULL) +
    THEME

  saveRDS(list(df = df, labels = labels, comparison = comparison_dir, dpi = dpi_label), file.path(plots_data_dir, paste0("3_volcano_data_", species, "_", comparison_dir, "_dpi", dpi_label, ".rds")))
  ggsave(file.path(plots_dir, paste0("3_volcano_", species, "_", comparison_dir, "_dpi", dpi_label, ".pdf")), p_volc, width = PLOT_WIDTH, height = PLOT_HEIGHT)
}

cat("Wrote volcano plots.\n")

# ============================================================================
# 4. Temporal trajectory — top candidate genes
# ============================================================================

for (cd in comparison_dirs) {
  temporal_file <- list.files(file.path(section7_dir, cd), pattern = "_temporal_summary\\.tsv$", full.names = TRUE)
  lfc_file <- list.files(file.path(section7_dir, cd), pattern = "_infection_LFC_matrix\\.tsv$", full.names = TRUE)

  if (length(temporal_file) == 0 || length(lfc_file) == 0) {
    message("Skipping temporal trajectory for ", cd, " — missing summary or LFC matrix.")
    next
  }

  temporal <- read.delim(temporal_file[1], stringsAsFactors = FALSE)
  lfc_mat <- read.delim(lfc_file[1], stringsAsFactors = FALSE)

  top_genes <- temporal %>%
    dplyr::filter(passed_candidate_filter_any_DPI) %>%
    dplyr::arrange(desc(max_abs_log2FC)) %>%
    slice_head(n = TOP_N_TRAJECTORY) %>%
    pull(gene_id)

  if (length(top_genes) == 0) next

  traj <- lfc_mat %>%
    dplyr::filter(gene_id %in% top_genes) %>%
    pivot_longer(cols = starts_with("dpi"), names_to = "dpi", values_to = "log2FC") %>%
    mutate(dpi = as.numeric(gsub("dpi", "", dpi)), gene_label = label_genes(gene_id))
  
  plot_title <- paste(species_label, label_comparison(cd), "top candidate trajectories")

  p_traj <- ggplot(traj, aes(dpi, log2FC, color = gene_label, group = gene_id)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = stringr::str_wrap(plot_title, width = 65), x = "dpi", y = "log2 fold change", color = "Gene") +
    THEME + theme(plot.title = element_text(hjust = 0.5))

  saveRDS(list(traj = traj, comparison = cd, top_genes = top_genes), file.path(plots_data_dir, paste0("4_trajectory_data_", species, "_", cd, ".rds")))
  ggsave(file.path(plots_dir, paste0("4_trajectory_", species, "_", cd, ".pdf")), p_traj, width = PLOT_WIDTH, height = PLOT_HEIGHT)
}

cat("Wrote temporal trajectory plots.\n")

# ============================================================================
# 5. Candidate overlap Venn — cattle only
# ============================================================================

if (species == "cattle" && length(comparison_dirs) > 1) {
  candidate_sets <- lapply(comparison_dirs, function(cd) {
    f <- list.files(file.path(section7_dir, cd), pattern = "_candidate_genes_all_DPI\\.tsv$", full.names = TRUE)
    if (length(f) == 0) return(character(0))
    unique(read.delim(f[1], stringsAsFactors = FALSE)$gene_id)
  })

  keep <- lengths(candidate_sets) > 0
  candidate_sets <- candidate_sets[keep]
  venn_labels <- label_comparison(comparison_dirs)[keep]
  venn_labels <- vapply(venn_labels, function(x) paste(strwrap(x, width = 24), collapse = "\n"), character(1))
  names(candidate_sets) <- venn_labels

  if (length(candidate_sets) >= 2) {
    p_venn <- ggVennDiagram(candidate_sets, category.names = venn_labels, set_size = 4, label_size = 4) +
      scale_fill_gradient(low = "white", high = "firebrick") +
      coord_cartesian(xlim = c(-0.15, 1.2), ylim = c(-0.15, 1.25), expand = FALSE) +
      labs(title = stringr::str_wrap(paste(species_label, "candidate gene overlap across comparisons"), width = 60)) +
      THEME +
      theme(plot.title = element_text(hjust = 0.5))

    set_layer <- which(vapply(p_venn$layers, function(layer) inherits(layer$geom, "GeomText") || inherits(layer$geom, "GeomLabel"), logical(1)))
    if (length(set_layer) > 0) {
      p_venn$layers[[set_layer[1]]]$aes_params$vjust <- 0
      p_venn$layers[[set_layer[1]]]$position <- position_nudge(y = 0.05)
    }
    
    saveRDS(list(candidate_sets = candidate_sets, labels = venn_labels),
        file.path(plots_data_dir, paste0("5_venn_data_", species, ".rds")))

    ggsave(file.path(plots_dir, paste0("5_candidate_overlap_", species, ".pdf")), p_venn, width = 12, height = 9)
    cat("Wrote candidate overlap Venn.\n")
  }
} else {
  cat("Skipping candidate overlap Venn (single comparison for", species_label, ").\n")
}

# ============================================================================
# 6. emapplot — GO ORA term relationships among enriched terms
# ============================================================================

go_result_files <- list.files(enrichment_dir, pattern = "_GO_results\\.rds$", recursive = TRUE, full.names = TRUE)
go_result_files <- go_result_files[!grepl("_GSEA_", go_result_files)]

for (grf in go_result_files) {
  go_results <- readRDS(grf)
  fname <- basename(grf)
  tag <- sub("_GO_results\\.rds$", "", fname)

  m <- regexec("^[^_]+_(.+)_dpi([^_]+)_(up|down)$", tag)
  parts <- regmatches(tag, m)[[1]]
  if (length(parts) != 4) {
    message("Could not parse GO filename: ", fname)
    next
  }

  comparison <- parts[2]
  dpi_label <- parts[3]
  direction <- parts[4]
  comparison_label <- label_comparison(comparison)

  for (src in names(go_results)) {
    res <- go_results[[src]]
    if (nrow(as.data.frame(res)) < 3) next

    res_sim <- tryCatch(pairwise_termsim(res), error = function(e) NULL)
    if (is.null(res_sim)) next

    saveRDS(list(result = res_sim, tag = tag, source = src, comparison = comparison, dpi = dpi_label, direction = direction), file.path(plots_data_dir, paste0("6_emap_data_", tag, "_", src, ".rds")))

    plot_title <- stringr::str_wrap(paste(species_label, comparison_label, "dpi", dpi_label, direction, "GO (", src, ") term map"), width = 80)
    p_emap <- tryCatch(emapplot(res_sim, showCategory = 20) + ggtitle(plot_title) +
                 theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                        legend.text = element_text(size = 12),
                        legend.title = element_text(size = 13, face = "bold"),
                        legend.key.size = unit(1.3, "cm")), error = function(e) NULL)

    if (!is.null(p_emap)) ggsave(file.path(plots_dir, paste0("6_emap_", tag, "_", src, ".pdf")), p_emap, width = 15, height = 15)
  }
}

cat("Wrote emapplots.\n")

# ============================================================================
# 7. Cross-dpi comparison dotplot — GO ORA enrichment patterns across time
# ============================================================================

for (cd in comparison_dirs) {
  for (dir_label in c("up", "down")) {
    go_files <- list.files(file.path(enrichment_dir, cd), pattern=paste0("_", dir_label, "_GO_results\\.rds$"), recursive=TRUE, full.names=TRUE)
    if (!length(go_files)) next

    dpi <- suppressWarnings(as.numeric(sub(".*_dpi([^_]+)_.*", "\\1", basename(go_files))))
    o <- order(dpi); go_files <- go_files[o]; dpi <- dpi[o]

    merged_list <- list()
    for (i in seq_along(go_files)) {
      x <- readRDS(go_files[i])
      for (src in names(x)) {
        if (is.null(x[[src]]) || !nrow(as.data.frame(x[[src]]))) next
        r <- x[[src]]
        r@result <- r@result[!is.na(r@result$GeneRatio) & !is.na(r@result$BgRatio) & !is.na(r@result$p.adjust), ]
        if (!nrow(r@result)) next
        r@result <- r@result[order(r@result$p.adjust), ][seq_len(min(5,nrow(r@result))), ]
        merged_list[[paste0("dpi",dpi[i],"_",src)]] <- r
      }
    }

    if (!length(merged_list)) next
    merged <- tryCatch(merge_result(merged_list), error=function(e) { message("merge_result: ",e$message); NULL })
    if (is.null(merged)) next

    res <- merged@compareClusterResult
    res <- res[!is.na(res$GeneRatio) & !is.na(res$p.adjust), ]
    clusters <- names(merged_list)
    dpi_num <- as.numeric(sub("^dpi([0-9]+)_.*$","\\1",clusters))
    src <- sub("^dpi[0-9]+_","",clusters)
    merged@compareClusterResult <- res
    merged@compareClusterResult$Cluster <- factor(merged@compareClusterResult$Cluster, levels=clusters[order(match(src,c("direct","ortholog")),dpi_num)])

    p_cmp <- dotplot(merged, showCategory=5) +
      ggtitle(paste(species_label,label_comparison(cd),dir_label,"GO across dpi")) +
      THEME +
      theme(plot.title=element_text(size=11),axis.text.x=element_text(angle=45,hjust=1))

    saveRDS(list(result=merged,plot_data=res,comparison=cd,direction=dir_label), file.path(plots_data_dir,paste0("7_dpi_comparison_data_",species,"_",cd,"_",dir_label,".rds")))
    ggsave(file.path(plots_dir,paste0("7_dpi_comparison_",species,"_",cd,"_",dir_label,".pdf")),p_cmp,width=12,height=10)
  }
}

cat("Wrote cross-dpi comparison dotplots.\n")

# ============================================================================
# 8. Cross-comparison overview plots — sourced from TABLES
# ============================================================================

# ---- 8a. Candidates over time ----------------------------------------------

comparison_summary_file <- file.path(tables_dir, paste0("2_comparison_summary_", species, ".tsv"))

if (file.exists(comparison_summary_file)) {
  cs <- read.delim(comparison_summary_file, stringsAsFactors = FALSE)

  cs_long <- cs %>%
    mutate(candidates_down = -candidates_down, comparison = label_comparison(comparison)) %>%
    pivot_longer(cols = c(candidates_up, candidates_down), names_to = "direction", values_to = "count") %>%
    mutate(direction = factor(ifelse(direction == "candidates_up", "Up", "Down"), levels = c("Up", "Down")), dpi = as.numeric(dpi))
  
  saveRDS(cs_long, file.path(plots_data_dir, paste0("8a_plot_data_", species, ".rds")))

  p_overview <- ggplot(cs_long, aes(dpi, count, fill = direction)) +
    geom_col(position = "identity") +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
    facet_wrap(~comparison) +
    labs(title = paste(species_label, "candidate genes over time"), x = "dpi", y = "Number of candidate genes", fill = NULL) +
    THEME

  ggsave(file.path(plots_dir, paste0("8a_candidates_overview_", species, ".pdf")), p_overview, width = 10, height = 7)
  cat("Wrote candidates-over-time overview.\n")
} else {
  message("Skipping candidates overview — missing ", comparison_summary_file)
}

# ---- 8b. Cross-comparison gene trajectory overlay -------------------------

if (species == "cattle") {
  traj_file <- file.path(tables_dir, paste0("5_temporal_trajectories_", species, ".tsv"))

  if (file.exists(traj_file)) {
    tt <- read.delim(traj_file, stringsAsFactors = FALSE)

    top_overlay_genes <- tt %>%
      group_by(gene_id) %>%
      dplyr::mutate(n_comparisons = n_distinct(comparison)) %>%
      ungroup() %>%
      dplyr::filter(n_comparisons >= 2) %>%
      distinct(gene_id, gene_symbol, entrez_id, max_abs_log2FC) %>%
      dplyr::arrange(desc(max_abs_log2FC)) %>%
      slice_head(n = TOP_N_OVERLAY_GENES) %>%
      pull(gene_id)

    if (length(top_overlay_genes) > 0) {
      overlay_long <- tt %>%
        filter(gene_id %in% top_overlay_genes) %>%
        pivot_longer(cols = matches("^dpi[0-9]+$"), names_to = "dpi", values_to = "log2FC")%>%
        mutate(dpi = as.numeric(gsub("dpi", "", dpi)), comparison_label = label_comparison(comparison), 
               facet_label = ifelse(!is.na(gene_symbol), gene_symbol, gene_id))

      saveRDS(list(top_overlay_genes = top_overlay_genes, overlay_long = overlay_long),
              file.path(plots_data_dir, paste0("8b_plot_data_", species, ".rds")))
      
      p_overlay <- ggplot(overlay_long, aes(dpi, log2FC, color = comparison_label, group = comparison)) +
        geom_line(alpha = 0.8) +
        geom_point(size = 1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        facet_wrap(~facet_label, scales = "free_y") +
        labs(title = paste(species_label, "gene response across comparisons"), x = "dpi", y = "log2 fold change", color = "Comparison") +
        THEME

      ggsave(file.path(plots_dir, paste0("8b_cross_comparison_trajectories_", species, ".pdf")), p_overlay, width = 11, height = 8)
      cat("Wrote cross-comparison trajectory overlay.\n")
    } else {
      message("Skipping cross-comparison trajectory overlay — no genes recur across >=2 comparisons.")
    }
  } else {
    message("Skipping cross-comparison trajectory overlay — missing ", traj_file)
  }
} else {
  cat("Skipping cross-comparison trajectory overlay (single comparison for", species_label, ").\n")
}

# ---- 8c. Top recurring enriched terms --------------------------------------

ora_file <- file.path(tables_dir, paste0("3_master_enrichment_ora_", species, ".tsv"))

if (file.exists(ora_file)) {
  ora <- read.delim(ora_file, stringsAsFactors = FALSE)

  recurring <- ora %>%
    filter(p.adjust <= TERM_SIG_THRESHOLD) %>%
    distinct(Description, term_type, comparison, dpi, direction) %>%
    count(Description, term_type, name = "n_significant_results") %>%
    arrange(desc(n_significant_results)) %>%
    slice_head(n = TOP_N_RECURRING_TERMS)

  if (nrow(recurring) > 0) {
    saveRDS(recurring, file.path(plots_data_dir, paste0("8c_plot_data_", species, ".rds")))

    p_recurring <- ggplot(recurring, aes(x = reorder(Description, n_significant_results), y = n_significant_results, fill = term_type)) +
      geom_col() +
      coord_flip() +
      labs(title = paste(species_label, "most recurring enriched terms"), x = NULL, y = "Number of significant (comparison, dpi, direction) results", fill = "Type") +
      THEME

    ggsave(file.path(plots_dir, paste0("8c_recurring_terms_", species, ".pdf")), p_recurring, width = 10, height = 8)
    cat("Wrote recurring terms plot.\n")
  }
} else {
  message("Skipping recurring terms plot — missing ", ora_file)
}

# ---- 8d. Hub-gene plot -----------------------------------------------------

crosswalk_file <- file.path(tables_dir, paste0("4_gene_pathway_crosswalk_", species, ".tsv"))

if (file.exists(crosswalk_file)) {
  cw <- read.delim(crosswalk_file, stringsAsFactors = FALSE)

  hub_genes <- cw %>%
    distinct(gene_id, gene_symbol, entrez_id, ID) %>%
    count(gene_id, gene_symbol, entrez_id, name = "n_terms") %>%
    arrange(desc(n_terms)) %>%
    slice_head(n = TOP_N_HUB_GENES) %>%
    mutate(label = ifelse(!is.na(gene_symbol), gene_symbol, gene_id))

  if (nrow(hub_genes) > 0) {
    saveRDS(hub_genes, file.path(plots_data_dir, paste0("8d_plot_data_", species, ".rds")))
    p_hub <- ggplot(hub_genes, aes(x = reorder(label, n_terms), y = n_terms)) +
      geom_col(fill = "darkorange") +
      coord_flip() +
      labs(title = paste(species_label, "candidate genes in the most enriched terms"), x = NULL, y = "Number of distinct GO/KEGG terms") +
      THEME

    ggsave(file.path(plots_dir, paste0("8d_hub_genes_", species, ".pdf")), p_hub, width = 9, height = 8)
    cat("Wrote hub-gene plot.\n")
  }
} else {
  message("Skipping hub-gene plot — missing ", crosswalk_file)
}

# ============================================================================
# 9. Recurring candidate genes across comparisons
# ============================================================================

recurring_candidate_files <- lapply(
  comparison_dirs,
  function(cd) {
    f <- list.files(
      file.path(section7_dir, cd),
      pattern = "_candidate_genes_all_DPI\\.tsv$",
      full.names = TRUE
    )

    if (!length(f)) return(NULL)

    x <- read.delim(
      f[1],
      stringsAsFactors = FALSE
    )

    if (!"gene_id" %in% names(x)) return(NULL)

    x %>%
      dplyr::select(gene_id) %>%
      dplyr::distinct() %>%
      dplyr::mutate(comparison = cd)
  }
)

recurring_candidate_table <- bind_rows(recurring_candidate_files)

if (nrow(recurring_candidate_table)) {

  recurring_genes <- recurring_candidate_table %>%
    group_by(gene_id) %>%
    summarise(
      n_comparisons = n_distinct(comparison),
      comparisons = paste(
        unique(comparison),
        collapse = ", "
      ),
      .groups = "drop"
    ) %>%
    arrange(
      desc(n_comparisons),
      gene_id
    )

  write.table(
    recurring_genes,
    file.path(
      plots_data_dir,
      paste0(
        "9_recurring_candidate_genes_",
        species,
        ".tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  top_recurring <- recurring_genes %>%
    dplyr::filter(n_comparisons >= 2) %>%
    slice_head(n = 30) %>%
    dplyr::mutate(
      label = label_genes(gene_id)
    )

  if (nrow(top_recurring)) {

    p_recurring_genes <- ggplot(
      top_recurring,
      aes(
        x = reorder(label, n_comparisons),
        y = n_comparisons
      )
    ) +
      geom_col(
        fill = "#0072B2",
        width = 0.75
      ) +
      geom_text(
        aes(label = n_comparisons),
        hjust = -0.2,
        size = 3.2
      ) +
      coord_flip() +
      scale_y_continuous(
        breaks = seq_len(length(comparison_dirs)),
        limits = c(0, max(top_recurring$n_comparisons) + 0.7),
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        title = paste(
          species_label,
          "recurring candidate genes"
        ),
        subtitle = "Genes identified as candidates in multiple comparisons",
        x = NULL,
        y = "Number of comparisons"
      ) +
      THEME

    ggsave(
      file.path(
        plots_dir,
        paste0(
          "9_recurring_candidate_genes_",
          species,
          ".pdf"
        )
      ),
      p_recurring_genes,
      width = 9,
      height = max(6, 0.25 * nrow(top_recurring) + 3)
    )

    cat("Wrote recurring candidate gene plot.\n")
  }
}

# ============================================================================
# 10. Candidate-gene response heatmap
# ============================================================================

candidate_stats_list <- list()

for (ef in expl_files) {

  comparison_dir <- basename(dirname(ef))
  fname <- basename(ef)

  m <- regmatches(
    fname,
    regexec(
      "_dpi([^_]+)_DESeq2_exploratory\\.tsv$",
      fname
    )
  )

  if (length(m[[1]]) != 2) next

  dpi_label <- m[[1]][2]

  df <- read.delim(
    ef,
    stringsAsFactors = FALSE
  )

  if (!all(
    c("gene_id", "log2FoldChange", "padj") %in% names(df)
  )) next

  df <- df %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj)
    ) %>%
    dplyr::mutate(
      comparison = comparison_dir,
      dpi = dpi_label,
      dpi_num = suppressWarnings(as.numeric(dpi_label))
    )

  candidate_stats_list[[length(candidate_stats_list) + 1]] <- df
}

if (length(candidate_stats_list)) {

  candidate_heatmap_data <- bind_rows(candidate_stats_list)

  # Keep genes that are candidates at least once
  if ("candidate" %in% names(candidate_heatmap_data)) {

    candidate_heatmap_data <- candidate_heatmap_data %>%
      dplyr::filter(candidate)

  }

  if (nrow(candidate_heatmap_data)) {

    # Keep the strongest genes overall
    top_candidate_genes <- candidate_heatmap_data %>%
      group_by(gene_id) %>%
      summarise(
        max_abs_log2FC = max(
          abs(log2FoldChange),
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(desc(max_abs_log2FC)) %>%
      slice_head(n = 40) %>%
      pull(gene_id)

    heatmap_long <- candidate_heatmap_data %>%
      dplyr::filter(gene_id %in% top_candidate_genes) %>%
      dplyr::mutate(
        comparison_label = label_comparison(comparison),
        sample_column = paste(
          comparison_label,
          "dpi",
          dpi,
          sep = "_"
        )
      ) %>%
      dplyr::select(
        gene_id,
        comparison,
        dpi,
        dpi_num,
        sample_column,
        log2FoldChange
      )

    # If a gene occurs more than once for the same comparison/dpi,
    # retain the first observation
    heatmap_long <- heatmap_long %>%
      dplyr::distinct(
        gene_id,
        comparison,
        dpi,
        .keep_all = TRUE
      )

    heatmap_matrix <- heatmap_long %>%
      dplyr::select(
        gene_id,
        sample_column,
        log2FoldChange
      ) %>%
      pivot_wider(
        names_from = sample_column,
        values_from = log2FoldChange
      ) %>%
      as.data.frame()

    rownames(heatmap_matrix) <- heatmap_matrix$gene_id
    heatmap_matrix$gene_id <- NULL

    heatmap_matrix <- as.matrix(heatmap_matrix)

    # Reorder genes by strongest overall response
    row_order <- order(
      apply(
        abs(heatmap_matrix),
        1,
        max,
        na.rm = TRUE
      ),
      decreasing = TRUE
    )

    heatmap_matrix <- heatmap_matrix[row_order, , drop = FALSE]

    rownames(heatmap_matrix) <- make.unique(
      label_genes(rownames(heatmap_matrix))
    )

    pdf(
      file.path(
        plots_dir,
        paste0(
          "10_candidate_response_heatmap_",
          species,
          ".pdf"
        )
      ),
      width = 13,
      height = max(
        8,
        0.28 * nrow(heatmap_matrix) + 5
      )
    )

    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(
        c(
          "#2166AC",
          "#67A9CF",
          "#D1E5F0",
          "white",
          "#FDDBC7",
          "#EF8A62",
          "#B2182B"
        )
      )(100),
      breaks = seq(
        -max(abs(heatmap_matrix), na.rm = TRUE),
        max(abs(heatmap_matrix), na.rm = TRUE),
        length.out = 101
      ),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = NA,
      main = paste(
        species_label,
        "candidate gene response across comparisons and dpi"
      ),
      fontsize_row = 8,
      fontsize_col = 8,
    )

    dev.off()

    saveRDS(
      heatmap_matrix,
      file.path(
        plots_data_dir,
        paste0(
          "10_candidate_response_heatmap_",
          species,
          ".rds"
        )
      )
    )

    cat("Wrote candidate response heatmap.\n")
  }
}

# ============================================================================
# 11. GO / KEGG enrichment barplots of most signicant enriched terms
# ============================================================================

go_plot_files <- list.files(
  enrichment_dir,
  pattern = "_GO_results\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

go_plot_files <- go_plot_files[
  !grepl("_GSEA_", go_plot_files)
]

if (length(go_plot_files)) {

  go_bar_list <- list()

  for (grf in go_plot_files) {

    x <- readRDS(grf)

    fname <- basename(grf)

    for (src in names(x)) {

      res <- x[[src]]

      if (is.null(res)) next

      res_df <- as.data.frame(res)

      if (!nrow(res_df)) next

      required <- c(
        "Description",
        "p.adjust",
        "Count"
      )

      if (!all(required %in% names(res_df))) next

      res_df <- res_df %>%
        dplyr::filter(
          !is.na(p.adjust),
          p.adjust <= 0.05
        ) %>%
        dplyr::arrange(p.adjust) %>%
        slice_head(n = 10)

      if (!nrow(res_df)) next

      go_bar_list[[length(go_bar_list) + 1]] <-
        res_df %>%
        dplyr::mutate(
          source = src,
          file = fname
        )
    }
  }

  if (length(go_bar_list)) {

    go_bar_data <- bind_rows(go_bar_list)

    # Plot the most significant terms overall
    go_bar_data <- go_bar_data %>%
      dplyr::arrange(p.adjust) %>%
      slice_head(n = 30) %>%
      dplyr::mutate(
        Description = factor(
          Description,
          levels = rev(unique(Description))
        )
      )

    p_go_bar <- ggplot(
      go_bar_data,
      aes(
        x = -log10(p.adjust),
        y = Description,
        fill = source
      )
    ) +
      geom_col(
        width = 0.7
      ) +
      labs(
        title = paste(
          species_label,
          "most significant enriched GO terms"
        ),
        x = expression(-log[10](adjusted~p)),
        y = NULL,
        fill = "Annotation source"
      ) +
      THEME

    ggsave(
      file.path(
        plots_dir,
        paste0(
          "11_GO_enrichment_barplot_",
          species,
          ".pdf"
        )
      ),
      p_go_bar,
      width = 11,
      height = max(
        7,
        0.25 * nrow(go_bar_data) + 4
      )
    )

    saveRDS(
      go_bar_data,
      file.path(
        plots_data_dir,
        paste0(
          "11_GO_enrichment_barplot_",
          species,
          ".rds"
        )
      )
    )

    cat("Wrote GO enrichment barplot.\n")
  }
}

# ============================================================================
# 12. Enrichment activity over dpi - number of significant terms over dpi
# ============================================================================

if (length(go_plot_files)) {

  enrichment_over_time <- list()

  for (grf in go_plot_files) {

    fname <- basename(grf)

    m <- regmatches(
      fname,
      regexec(
        "^[^_]+_(.+)_dpi([^_]+)_(up|down)_GO_results\\.rds$",
        fname
      )
    )[[1]]

    if (length(m) != 4) next

    comparison <- m[2]
    dpi <- suppressWarnings(as.numeric(m[3]))
    direction <- m[4]

    x <- readRDS(grf)

    for (src in names(x)) {

      res <- x[[src]]

      if (is.null(res)) next

      res_df <- as.data.frame(res)

      if (!"p.adjust" %in% names(res_df)) next

      n_sig <- sum(
        !is.na(res_df$p.adjust) &
        res_df$p.adjust <= TERM_SIG_THRESHOLD
      )

      enrichment_over_time[[length(enrichment_over_time) + 1]] <-
        data.frame(
          comparison = comparison,
          dpi = dpi,
          direction = direction,
          source = src,
          n_significant_terms = n_sig,
          stringsAsFactors = FALSE
        )
    }
  }

  if (length(enrichment_over_time)) {

    enrichment_over_time <- bind_rows(
      enrichment_over_time
    ) %>%
      dplyr::mutate(
        comparison_label = label_comparison(comparison)
      )

    p_enrichment_time <- ggplot(
      enrichment_over_time,
      aes(
        x = dpi,
        y = n_significant_terms,
        color = source,
        group = source
      )
    ) +
      geom_line(
        linewidth = 0.9
      ) +
      geom_point(
        size = 2.5
      ) +
      facet_grid(
        direction ~ comparison_label
      ) +
      scale_x_continuous(
        breaks = sort(
          unique(enrichment_over_time$dpi)
        )
      ) +
      labs(
        title = paste(
          species_label,
          "GO enrichment over time"
        ),
        x = "dpi",
        y = "Number of significant GO terms",
        color = "Annotation source"
      ) +
      THEME + theme(strip.text = element_text(size = 7))

    ggsave(
      file.path(
        plots_dir,
        paste0(
          "12_GO_enrichment_over_time_",
          species,
          ".pdf"
        )
      ),
      p_enrichment_time,
      width = 12,
      height = 7
    )

    saveRDS(
      enrichment_over_time,
      file.path(
        plots_data_dir,
        paste0(
          "12_GO_enrichment_over_time_",
          species,
          ".rds"
        )
      )
    )

    cat("Wrote GO enrichment-over-time plot.\n")
  }
}

cat(
  "All plots written to:",
  plots_dir,
  "\n"
)
