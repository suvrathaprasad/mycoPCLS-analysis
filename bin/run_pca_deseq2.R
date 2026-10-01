#!/usr/bin/env Rscript

# ============================================================================
# run_pca_deseq2.R
#
# Final stage of the combined pipeline:
# featureCounts matrix -> PCA -> per-timepoint DESeq2.
#
# IMPORTANT:
# - Sequencing runs remain separate; collapseReplicates() is NOT used.
# - One biological sample exists per condition x dpi cell.
# - Dispersions are estimated across the full time course using:
#       ~ dpi + condition
# - A saturated cell-means model is then used:
#       ~ 0 + dpi_condition
# - P-values are exploratory because biological replication within each
#   condition x dpi cell is absent.
#
# Usage:
#   run_pca_deseq2.R <counts_tsv> <metadata_csv> <species> <outdir>
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4)
  stop("Usage: run_deseq2_pca.R <counts_tsv> <metadata_csv> <species> <outdir>")

counts_tsv <- args[1]
metadata_csv <- args[2]
species <- tolower(args[3])
outdir <- args[4]

suppressPackageStartupMessages({
  library(DESeq2)
  library(PCAtools)
  library(ggplot2)
  library(dplyr)
})

# ============================================================================
# 0. Output directories
# ============================================================================

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "PCA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "SECTION_6"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "SECTION_7"), recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. Load metadata
# ============================================================================

meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c("sample", "group", "condition", "dpi")
if (!all(required_columns %in% colnames(meta)))
  stop("metadata must contain columns: ", paste(required_columns, collapse = ", "))

if (!"replicate" %in% colnames(meta))
  meta$replicate <- ave(seq_len(nrow(meta)), meta$group, FUN = seq_along)

meta$sample <- as.character(meta$sample)
meta$group <- as.character(meta$group)
meta$condition <- as.character(meta$condition)
meta$dpi <- as.character(meta$dpi)
meta$replicate <- as.character(meta$replicate)
rownames(meta) <- meta$sample

if (anyDuplicated(meta$sample)) {
  dupes <- unique(meta$sample[duplicated(meta$sample)])
  stop("Duplicate sample IDs found: ", paste(dupes, collapse = ", "))
}

# ============================================================================
# 2. Import featureCounts matrix
# ============================================================================

counts_raw <- read.delim(
  counts_tsv,
  header = TRUE,
  comment.char = "#",
  check.names = FALSE,
  row.names = 1
)

colnames(counts_raw) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts_raw))

if (!all(meta$sample %in% colnames(counts_raw))) {
  missing_samples <- setdiff(meta$sample, colnames(counts_raw))
  stop("Missing count columns for: ", paste(missing_samples, collapse = ", "),
       "\nCheck that sample IDs in the metadata match the samplesheet.")
}

counts_mat <- as.matrix(counts_raw[, meta$sample, drop = FALSE])
storage.mode(counts_mat) <- "integer"

if (!identical(colnames(counts_mat), meta$sample))
  stop("Count matrix columns do not match metadata sample order.")

cat("Loaded count matrix:", nrow(counts_mat), "genes x", ncol(counts_mat),
    "sequencing runs\n")

# ============================================================================
# 3. Create run-level DESeqDataSet
#
# IMPORTANT:
# NO collapseReplicates().
# Every sequencing run remains an individual column.
# ============================================================================

dds_raw <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta, design = ~1)

# ============================================================================
# 4. Species-specific condition pairs
# ============================================================================

if (species == "cattle") {
  condition_pairs <- list(
    c(control = "NK_oP", treatment = "oP"),
    c(control = "NK_mP", treatment = "mP"),
    c(control = "oP", treatment = "mP"),
    c(control = "NK_oP", treatment = "NK_mP")
  )
} else if (species == "sheep") {
  condition_pairs <- list(c(control = "NK", treatment = "infection"))
} else {
  stop("species must be 'rind' or 'schaf', got: ", species)
}

# ============================================================================
# 5. PCA
#
# PCA 1:
# Biological-sample PCA. The two runs belonging to the same biological
# sample/group are averaged for visualization only.
#
# PCA 2:
# Every sequencing run is shown individually.
# ============================================================================

# ----------------------------------------------------------------------------
# 5a. Normalized expression
# ----------------------------------------------------------------------------

dds_pca <- dds_raw
dds_pca <- estimateSizeFactors(dds_pca)
norm_counts_pca <- counts(dds_pca, normalized = TRUE)
log_norm_counts <- log2(norm_counts_pca + 1)

# ----------------------------------------------------------------------------
# 5b. Biological-sample PCA
# ----------------------------------------------------------------------------

group_levels <- unique(meta$group)

log_norm_counts_bio <- sapply(group_levels, function(g) {
  runs <- meta$sample[meta$group == g]
  if (length(runs) == 0) stop("No runs found for group: ", g)
  rowMeans(log_norm_counts[, runs, drop = FALSE])
})

colnames(log_norm_counts_bio) <- group_levels

bio_meta <- meta[!duplicated(meta$group), , drop = FALSE]
rownames(bio_meta) <- bio_meta$group
bio_meta <- bio_meta[group_levels, , drop = FALSE]

bio_meta$dpi <- factor(
  bio_meta$dpi,
  levels = sort(unique(bio_meta$dpi), na.last = TRUE)
)

p_bio <- PCAtools::pca(log_norm_counts_bio, metadata = bio_meta)

pca_plot <- PCAtools::biplot(
  p_bio,
  colby = "condition",
  shape = "dpi",
  legendPosition = "right",
  title = paste("PCA", species)
)

ggsave(file.path(outdir, "PCA", paste0("pca_", species, ".pdf")),
       pca_plot, width = 8, height = 6)

saveRDS(
  list(
    log_norm_counts_bio = log_norm_counts_bio,
    bio_meta = bio_meta,
    dds_raw = dds_raw,
    method = "DESeq2_normalized_counts_mean_for_PCA_only"
  ),
  file.path(outdir, "PCA", paste0("dds_", species, "_collapsed.rds"))
)

write.csv(
  bio_meta,
  file.path(outdir, "PCA", paste0("collapsed_metadata_", species, ".csv")),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 5c. Technical-replicate PCA
# ----------------------------------------------------------------------------

pca_meta_raw <- meta
pca_meta_raw$dpi <- factor(
  pca_meta_raw$dpi,
  levels = sort(unique(pca_meta_raw$dpi), na.last = TRUE)
)

p_raw <- PCAtools::pca(log_norm_counts, metadata = pca_meta_raw)

pca_plot_raw <- PCAtools::biplot(
  p_raw,
  colby = "condition",
  shape = "dpi",
  lab = pca_meta_raw$replicate,
  legendPosition = "right",
  title = paste("PCA with replicates of", species)
)

ggsave(
  file.path(outdir, "PCA", paste0("pca_", species, "_technical_replicates.pdf")),
  pca_plot_raw, width = 8, height = 6
)

# ============================================================================
# Helper function:
# Estimate DESeq2 dispersions using the full time-course model.
# ============================================================================

fit_pair_deseq2 <- function(dds_raw, meta, ctrl, trt) {

  # --------------------------------------------------------------------------
  # Select samples belonging to this condition pair
  # --------------------------------------------------------------------------

  pair_samples <- meta$sample[meta$condition %in% c(ctrl, trt)]

  if (length(pair_samples) < 2)
    stop("Fewer than 2 sequencing runs for pair ", trt, " vs ", ctrl)

  dds_pair <- dds_raw[, pair_samples, drop = FALSE]

  # --------------------------------------------------------------------------
  # Align metadata
  # --------------------------------------------------------------------------

  meta_pair <- meta[match(colnames(dds_pair), meta$sample), , drop = FALSE]

  if (anyNA(meta_pair$sample))
    stop("Could not align metadata to DESeq2 columns for ", trt, " vs ", ctrl)

  if (!identical(colnames(dds_pair), meta_pair$sample))
    stop("Metadata/sample alignment failed for ", trt, " vs ", ctrl)

  colData(dds_pair) <- S4Vectors::DataFrame(meta_pair)

  # --------------------------------------------------------------------------
  # Set factors
  # --------------------------------------------------------------------------

  dds_pair$condition <- factor(dds_pair$condition, levels = c(ctrl, trt))

  dds_pair$dpi <- factor(
    dds_pair$dpi,
    levels = sort(unique(as.numeric(as.character(dds_pair$dpi))))
  )

  # --------------------------------------------------------------------------
  # Estimate dispersions using the full time-course model
  # --------------------------------------------------------------------------

  design(dds_pair) <- ~ dpi + condition
  dds_pair <- estimateSizeFactors(dds_pair)
  dds_pair <- estimateDispersions(dds_pair)

  dispersion_function <- dispersionFunction(dds_pair)
  gene_wise_dispersion <- mcols(dds_pair)$dispGeneEst
  final_dispersion <- mcols(dds_pair)$dispersion

  # --------------------------------------------------------------------------
  # Saturated cell-means design
  # --------------------------------------------------------------------------

  dds_pair$dpi_condition <- factor(
    paste(dds_pair$condition, dds_pair$dpi, sep = "_dpi")
  )

  design(dds_pair) <- ~ 0 + dpi_condition

  # IMPORTANT:
  # Do NOT estimate dispersions again.
  # Reuse the dispersions estimated from ~ dpi + condition.
  # --------------------------------------------------------------------------

  mcols(dds_pair)$dispGeneEst <- gene_wise_dispersion
  mcols(dds_pair)$dispersion <- final_dispersion
  dispersionFunction(dds_pair) <- dispersion_function

  dds_pair <- nbinomWaldTest(dds_pair)

  list(
    dds = dds_pair,
    meta = meta_pair,
    dispersion_function = dispersion_function
  )
}

# ============================================================================
# 6. Section 6: padj-gated DE
#
# Criteria:
#   |log2FC| >= 1
#   padj <= 0.05
# ============================================================================

lfc_threshold <- 1
padj_threshold <- 0.05

for (pair in condition_pairs) {

  ctrl <- pair[["control"]]
  trt <- pair[["treatment"]]
  comparison_dir <- paste0(trt, "_vs_", ctrl)

  dir.create(file.path(outdir, "SECTION_6", comparison_dir),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "SECTION_7", comparison_dir),
             recursive = TRUE, showWarnings = FALSE)

  # --------------------------------------------------------------------------
  # Fit pair
  # --------------------------------------------------------------------------

  mod <- fit_pair_deseq2(dds_raw = dds_raw, meta = meta, ctrl = ctrl, trt = trt)
  dds_pair <- mod$dds
  meta_pair <- mod$meta
  norm_counts_pair <- counts(dds_pair, normalized = TRUE)

  # --------------------------------------------------------------------------
  # Timepoints
  # --------------------------------------------------------------------------

  dpi_values <- as.character(sort(unique(as.numeric(as.character(meta_pair$dpi)))))

  for (day in dpi_values) {

    trt_lvl <- paste0(trt, "_dpi", day)
    ctrl_lvl <- paste0(ctrl, "_dpi", day)

    if (!all(c(trt_lvl, ctrl_lvl) %in% levels(dds_pair$dpi_condition))) {
      message("Skipping ", trt, "_vs_", ctrl, " at dpi ", day, " - level missing.")
      next
    }

    trt_idx <- as.character(dds_pair$dpi_condition) == trt_lvl
    ctrl_idx <- as.character(dds_pair$dpi_condition) == ctrl_lvl

    if (!any(trt_idx) || !any(ctrl_idx)) {
      message("Skipping ", trt, "_vs_", ctrl, " at dpi ", day,
              " - no samples in one level.")
      next
    }

    # ------------------------------------------------------------------------
    # Sample means
    # ------------------------------------------------------------------------

    mean_trt <- rowMeans(norm_counts_pair[, trt_idx, drop = FALSE])
    mean_ctrl <- rowMeans(norm_counts_pair[, ctrl_idx, drop = FALSE])

    # ------------------------------------------------------------------------
    # DESeq2 contrast
    # ------------------------------------------------------------------------

    res <- results(
      dds_pair,
      contrast = c("dpi_condition", trt_lvl, ctrl_lvl),
      alpha = padj_threshold,
      lfcThreshold = lfc_threshold,
      altHypothesis = "greaterAbs",
      independentFiltering = TRUE
    )

    # ------------------------------------------------------------------------
    # Checks
    # ------------------------------------------------------------------------

    if (nrow(res) != nrow(dds_pair))
      stop("Result row count mismatch for ", trt_lvl, " vs ", ctrl_lvl)

    if (length(mean_trt) != nrow(res) || length(mean_ctrl) != nrow(res)) {
      stop("Mean vector length mismatch for ", trt_lvl, " vs ", ctrl_lvl,
           "\nnrow(res) = ", nrow(res),
           "\nlength(mean_trt) = ", length(mean_trt),
           "\nlength(mean_ctrl) = ", length(mean_ctrl))
    }

    # ------------------------------------------------------------------------
    # Output table
    # ------------------------------------------------------------------------

    z <- data.frame(
      gene_id = rownames(res),
      mean_normalized_treatment = mean_trt,
      mean_normalized_control = mean_ctrl,
      as.data.frame(res),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    z$padj[!is.na(z$pvalue) & z$pvalue == 1] <- 1

    # ------------------------------------------------------------------------
    # Significant up/down genes
    # ------------------------------------------------------------------------

    z_up <- subset(z, !is.na(padj) & padj <= padj_threshold & log2FoldChange >= lfc_threshold)
    z_down <- subset(z, !is.na(padj) & padj <= padj_threshold & log2FoldChange <= -lfc_threshold)

    tag <- paste0(species, "_", trt, "_vs_", ctrl, "_dpi", day)
    section_dir <- file.path(outdir, "SECTION_6", comparison_dir)

    # ------------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------------

    write.table(z, file.path(section_dir, paste0(tag, "_DEGs.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(z_up, file.path(section_dir, paste0(tag, "_upDEGs.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(z_down, file.path(section_dir, paste0(tag, "_downDEGs.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    cat("Done:", tag, "-", nrow(z_up), "up,", nrow(z_down), "down\n")
  }

  # --------------------------------------------------------------------------
  # Save DESeq2 object
  # --------------------------------------------------------------------------

  saveRDS(
    dds_pair,
    file.path(outdir, "SECTION_6", comparison_dir,
              paste0("dds_", species, "_", trt, "_vs_", ctrl, ".rds"))
  )
}

cat("Section 6 (padj-gated) contrasts complete.\n")

# ============================================================================
# 7. Exploratory per-timepoint DE
#
# Candidate definition:
#   mean normalized count >= 10 in at least one condition
#   AND |log2FC| >= 1.5
#
# Technical-replicate consistency is QC only.
# ============================================================================

min_mean_normalized_count <- 10
lfc_candidate_threshold <- 1.5

# ============================================================================
# 7a. Technical-replicate QC
# ============================================================================

dds_raw_qc <- dds_raw
dds_raw_qc <- estimateSizeFactors(dds_raw_qc)
raw_norm_counts <- counts(dds_raw_qc, normalized = TRUE)

group_levels <- unique(meta$group)
tech_diff_list <- vector("list", length(group_levels))
names(tech_diff_list) <- group_levels

for (g in group_levels) {

  runs <- meta$sample[meta$group == g]

  if (length(runs) != 2) {
    warning("Group ", g, " has ", length(runs),
            " sequencing runs; expected exactly 2.")
    next
  }

  tech_diff_list[[g]] <- abs(
    log2(raw_norm_counts[, runs[1]] + 1) -
    log2(raw_norm_counts[, runs[2]] + 1)
  )
}

tech_diff_list <- tech_diff_list[!vapply(tech_diff_list, is.null, logical(1))]

if (length(tech_diff_list) > 0) {
  tech_diff_matrix <- do.call(cbind, tech_diff_list)
} else {
  tech_diff_matrix <- NULL
}

# ============================================================================
# 7b. Per-condition-pair DESeq2 analysis
# ============================================================================

for (pair in condition_pairs) {

  ctrl <- pair[["control"]]
  trt <- pair[["treatment"]]
  comparison_dir <- paste0(trt, "_vs_", ctrl)

  # --------------------------------------------------------------------------
  # Fit model
  # --------------------------------------------------------------------------

  mod <- fit_pair_deseq2(dds_raw = dds_raw, meta = meta, ctrl = ctrl, trt = trt)
  dds_pair <- mod$dds
  meta_pair <- mod$meta
  norm_counts_pair <- counts(dds_pair, normalized = TRUE)

  # --------------------------------------------------------------------------
  # Timepoints
  # --------------------------------------------------------------------------

  dpi_values <- as.character(sort(unique(as.numeric(as.character(meta_pair$dpi)))))

  # --------------------------------------------------------------------------
  # LFC matrix
  # --------------------------------------------------------------------------

  lfc_matrix <- matrix(
    NA_real_,
    nrow = nrow(dds_pair),
    ncol = length(dpi_values),
    dimnames = list(rownames(dds_pair), paste0("dpi", dpi_values))
  )

  all_results <- vector("list", length(dpi_values))
  names(all_results) <- paste0("dpi", dpi_values)

  # ==========================================================================
  # 7c. Per-dpi contrasts
  # ==========================================================================

  for (i in seq_along(dpi_values)) {

    day <- dpi_values[i]
    trt_lvl <- paste0(trt, "_dpi", day)
    ctrl_lvl <- paste0(ctrl, "_dpi", day)

    if (!all(c(trt_lvl, ctrl_lvl) %in% levels(dds_pair$dpi_condition))) {
      message("Skipping ", trt, "_vs_", ctrl, " at dpi ", day, " - level missing.")
      next
    }

    trt_idx <- as.character(dds_pair$dpi_condition) == trt_lvl
    ctrl_idx <- as.character(dds_pair$dpi_condition) == ctrl_lvl

    if (!any(trt_idx) || !any(ctrl_idx)) {
      message("Skipping ", trt, "_vs_", ctrl, " at dpi ", day,
              " - no samples in one level.")
      next
    }

    # ------------------------------------------------------------------------
    # Normalized expression means
    # ------------------------------------------------------------------------

    mean_trt <- rowMeans(norm_counts_pair[, trt_idx, drop = FALSE])
    mean_ctrl <- rowMeans(norm_counts_pair[, ctrl_idx, drop = FALSE])

    # ------------------------------------------------------------------------
    # DESeq2 exploratory contrast
    # ------------------------------------------------------------------------

    res <- results(
      dds_pair,
      contrast = c("dpi_condition", trt_lvl, ctrl_lvl),
      lfcThreshold = lfc_candidate_threshold,
      altHypothesis = "greaterAbs",
      independentFiltering = TRUE
    )

    # ------------------------------------------------------------------------
    # Technical replicate QC
    # ------------------------------------------------------------------------

    ctrl_group <- unique(meta_pair$group[
      meta_pair$condition == ctrl & meta_pair$dpi == day
    ])

    trt_group <- unique(meta_pair$group[
      meta_pair$condition == trt & meta_pair$dpi == day
    ])

    if (length(ctrl_group) == 1 && !is.null(tech_diff_matrix) &&
        ctrl_group %in% colnames(tech_diff_matrix)) {
      tech_ctrl <- as.numeric(tech_diff_matrix[, ctrl_group])
    } else {
      tech_ctrl <- rep(NA_real_, nrow(dds_pair))
    }

    if (length(trt_group) == 1 && !is.null(tech_diff_matrix) &&
        trt_group %in% colnames(tech_diff_matrix)) {
      tech_trt <- as.numeric(tech_diff_matrix[, trt_group])
    } else {
      tech_trt <- rep(NA_real_, nrow(dds_pair))
    }

    tech_max_diff <- pmax(tech_ctrl, tech_trt, na.rm = TRUE)
    tech_max_diff[is.infinite(tech_max_diff)] <- NA_real_

    # ------------------------------------------------------------------------
    # Build result table
    # ------------------------------------------------------------------------

    z <- data.frame(
      gene_id = rownames(res),
      mean_normalized_control = mean_ctrl,
      mean_normalized_infected = mean_trt,
      as.data.frame(res),
      technical_max_abs_log2_difference = tech_max_diff,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    z$padj[!is.na(z$pvalue) & z$pvalue == 1] <- 1

    # ------------------------------------------------------------------------
    # Candidate filtering
    # ------------------------------------------------------------------------

    expressed <- (
      z$mean_normalized_control >= min_mean_normalized_count |
      z$mean_normalized_infected >= min_mean_normalized_count
    )

    sufficient_effect <- (
      !is.na(z$log2FoldChange) &
      abs(z$log2FoldChange) >= lfc_candidate_threshold
    )

    z$candidate <- expressed & sufficient_effect

    # ------------------------------------------------------------------------
    # Technical replicate consistency
    # ------------------------------------------------------------------------

    z$technical_replicates_consistent <- (
      !is.na(z$technical_max_abs_log2_difference) &
      z$technical_max_abs_log2_difference <= 1
    )

    z$technical_replicate_status <- ifelse(
      is.na(z$technical_max_abs_log2_difference),
      "not_available",
      ifelse(z$technical_max_abs_log2_difference <= 1, "consistent", "inspect")
    )

    z$statistical_evidence <- "exploratory_DESeq2_pooled_dispersion"

    # ------------------------------------------------------------------------
    # Candidate up/down
    # ------------------------------------------------------------------------

    z_up <- subset(z, candidate & log2FoldChange >= lfc_candidate_threshold)
    z_down <- subset(z, candidate & log2FoldChange <= -lfc_candidate_threshold)

    tag <- paste0(species, "_", trt, "_vs_", ctrl, "_dpi", day)
    section_dir <- file.path(outdir, "SECTION_7", comparison_dir)

    # ------------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------------

    write.table(
      z,
      file.path(section_dir, paste0(tag, "_DESeq2_exploratory.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
      z_up,
      file.path(section_dir, paste0(tag, "_candidate_up.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
      z_down,
      file.path(section_dir, paste0(tag, "_candidate_down.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    # ------------------------------------------------------------------------
    # Save for combined table and LFC matrix
    # ------------------------------------------------------------------------

    all_results[[paste0("dpi", day)]] <- z
    lfc_matrix[, i] <- res$log2FoldChange

    cat("Done:", tag, "- candidates:", nrow(z_up), "up,", nrow(z_down), "down\n")
  }

  # ==========================================================================
  # 7d. Infection-effect matrix
  # ==========================================================================

  lfc_output <- data.frame(gene_id = rownames(lfc_matrix), lfc_matrix, check.names = FALSE)

  write.table(
    lfc_output,
    file.path(
      outdir, "SECTION_7", comparison_dir,
      paste0(species, "_", trt, "_vs_", ctrl, "_infection_LFC_matrix.tsv")
    ),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # ==========================================================================
  # 7e. Combined candidate table
  # ==========================================================================

  combined_candidates <- bind_rows(all_results, .id = "dpi") %>%
    filter(candidate) %>%
    arrange(gene_id, desc(abs(log2FoldChange)))

  write.table(
    combined_candidates,
    file.path(
      outdir, "SECTION_7", comparison_dir,
      paste0(species, "_", trt, "_vs_", ctrl, "_candidate_genes_all_DPI.tsv")
    ),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # ==========================================================================
  # 7f. Temporal summary
  # ==========================================================================

  max_abs_lfc <- apply(abs(lfc_matrix), 1, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })

  max_lfc_dpi <- apply(lfc_matrix, 1, function(x) {
    if (all(is.na(x))) NA_character_ else colnames(lfc_matrix)[which.max(abs(x))]
  })

  temporal_summary <- data.frame(
    gene_id = rownames(lfc_matrix),
    max_abs_log2FC = max_abs_lfc,
    DPI_of_max_abs_log2FC = max_lfc_dpi,
    stringsAsFactors = FALSE
  )

  candidate_genes <- unique(combined_candidates$gene_id)

  temporal_summary$passed_candidate_filter_any_DPI <-
    temporal_summary$gene_id %in% candidate_genes

  write.table(
    temporal_summary,
    file.path(
      outdir, "SECTION_7", comparison_dir,
      paste0(species, "_", trt, "_vs_", ctrl, "_temporal_summary.tsv")
    ),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # ==========================================================================
  # Save DESeq2 object
  # ==========================================================================

  saveRDS(
    dds_pair,
    file.path(
      outdir, "SECTION_7", comparison_dir,
      paste0("dds_", species, "_", trt, "_vs_", ctrl, ".rds")
    )
  )
}

cat("All section 7 contrasts complete. Outputs written to:", outdir, "\n")
