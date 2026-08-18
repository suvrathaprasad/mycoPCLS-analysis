#!/usr/bin/env Rscript
# ============================================================================
# run_deseq2_pca.R
#
# Final stage of the combined pipeline: import the featureCounts count
# matrix (STAR alignment + featureCounts quantification, upstream in
# main.nf), collapse technical replicates, PCA across all samples, then
# per-timepoint NK-vs-infection differential expression.
#
# There is 1 animal per condition x dpi cell (confirmed). No per-cell
# biological replication anywhere in this dataset. Per-timepoint DE
# therefore borrows dispersion across all dpi within a condition_pair
# (real residual df there, borrowing strength across timepoints), then
# switch to a saturated cell-means design (~0 + dpi_condition) and
# pull out each day's specific contrast by name, reusing the dispersion
# already estimated rather than re-estimating it on the saturated
# (0 residual df) design.
#
# Usage:
#   run_deseq2_pca.R <counts_tsv> <metadata_csv> <species> <outdir>
#
#   counts_tsv    - merged featureCounts count matrix from main.nf's
#                   FEATURECOUNTS step (counts.tsv: Geneid + one column per
#                   sample, column names are "<sample>.Aligned.sortedByCoord.out.bam")
#   metadata_csv  - metadata_rind.csv or metadata_schaf.csv. Must have
#                   'sample', 'group' (technical-replicate grouping),
#                   'condition', and 'dpi' columns.
#   species       - "rind" or "schaf" — used for output file naming
#   outdir        - output directory (created if it doesn't exist)
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: run_deseq2_pca.R <counts_tsv> <metadata_csv> <species> <outdir>")
}
counts_tsv   <- args[1]
metadata_csv <- args[2]
species      <- args[3]   # "rind" or "schaf" — output naming only
outdir       <- args[4]

suppressPackageStartupMessages({
  library(DESeq2)
  library(PCAtools)
  library(ggplot2)
  library(dplyr)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Create subdirectories for organised output
dir.create(file.path(outdir, "PCA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "SECTION_6"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "SECTION_7"), recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load metadata -----------------------------------------------------

meta <- read.csv(metadata_csv, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample

# ---- 2. Import counts from the featureCounts matrix ------------------------

counts_raw <- read.delim(counts_tsv, header = TRUE, comment.char = "#",
                          check.names = FALSE, row.names = 1)

# Column names are the BAM filenames STAR/featureCounts produced
# ("<sample>.Aligned.sortedByCoord.out.bam") so strip that back down to the
# plain sample ID so columns line up with metadata_csv's `sample` column.
colnames(counts_raw) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "",
                             colnames(counts_raw))

if (!all(meta$sample %in% colnames(counts_raw))) {
  stop("Missing count columns for: ",
       paste(setdiff(meta$sample, colnames(counts_raw)), collapse = ", "),
       " — check that sample IDs in the metadata match the samplesheet used for STAR/featureCounts.")
}

# Explicit reorder to metadata's sample order and featureCounts' column order
# follows whatever order Nextflow's .collect() handed it the BAMs in
counts_mat <- as.matrix(counts_raw[, meta$sample, drop = FALSE])
storage.mode(counts_mat) <- "integer"

# ---- 3. Collapse technical replicates ---------------------------------------
# Each biological sample was isolated/processed twice (technical replicate,
# same biological material, confirmed) so collapse before DE testing so the
# replicate pair isn't treated as two independent biological samples, which
# would silently overstate statistical power.

dds_raw <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta, design = ~1)

# 'group' column (condition + dpi) identifies which runs belong to the same
# biological sample and should be merged together.
dds <- collapseReplicates(dds_raw, groupby = dds_raw$group, run = colnames(dds_raw))
collapsed_meta <- as.data.frame(colData(dds))

# ---- 4. Species-specific condition pairs -------------------------------------
# Uncomment ONE block below matching the species you're running. Each entry
# is one NK-vs-infection contrast to run at every dpi found in the data.

## --- RIND (cattle): two independent challenge/control pairs -----------------
# condition_pairs <- list(
#   c(control = "NK_oP", treatment = "oP"),
#   c(control = "NK_mP", treatment = "mP"),
#   c(control = "NK_oP", treatment = "NK_mP"),
#   c(control = "oP",    treatment = "mP"   )
# )

## --- SCHAF (sheep): one challenge/control pair -------------------------------
condition_pairs <- list(
  c(control = "NK", treatment = "infection")
)

if (!exists("condition_pairs")) {
  stop("No species-specific condition_pairs block was uncommented in section 4. ",
       "Please edit the script and uncomment either the RIND or SCHAF block before running.")
}

# ---- 5. PCA plots (purely descriptive/QC) -------------------------------------

pca_meta <- collapsed_meta
pca_meta$dpi <- factor(pca_meta$dpi, levels = sort(unique(pca_meta$dpi)))

vsd <- vst(dds, blind = TRUE)
p <- PCAtools::pca(assay(vsd), metadata = pca_meta)

pca_plot <- PCAtools::biplot(
  p,
  colby = "condition",
  shape = "dpi",
  legendPosition = "right",
  title = paste("PCA", species)
)

ggsave(file.path(outdir, "PCA", paste0("pca_", species, ".pdf")), pca_plot, width = 8, height = 6)
saveRDS(dds, file.path(outdir, "PCA", paste0("dds_", species, "_collapsed.rds")))
write.csv(collapsed_meta, file.path(outdir, "PCA", paste0("collapsed_metadata_", species, ".csv")),
          row.names = FALSE)

# Same PCA, but on dds_raw (pre-collapseReplicates) every technical
# replicate shown as its own point, rather than summed into one biological
# sample.
pca_meta_raw <- as.data.frame(colData(dds_raw))
pca_meta_raw$dpi <- factor(pca_meta_raw$dpi, levels = sort(unique(pca_meta_raw$dpi)))

vsd_raw <- vst(dds_raw, blind = TRUE)
p_raw <- PCAtools::pca(assay(vsd_raw), metadata = pca_meta_raw)

pca_plot_raw <- PCAtools::biplot(
  p_raw,
  colby = "condition",
  shape = "dpi",
  # lab = paste0("rep", pca_meta_raw$replicate), # "rep" or "s" ?
  lab = pca_meta_raw$replicate,
  legendPosition = "right",
  title = paste("PCA with technical replicates of", species)
)

ggsave(file.path(outdir, "PCA", paste0("pca_", species, "_technical_replicates.pdf")), pca_plot_raw, width = 8, height = 6)

# ---- 6. Per-timepoint differential expression --------------------------------
# not used for pathway analysis: padj-gated lists here run
# ~20 genes up/down per contrast, below what GO/KEGG enrichment tools need
# to do anything meaningful with (most want 50+ for real power). Section 7
# below is the one actually driving candidate-gene selection now.
#
# confirmed: 1 animal per condition x dpi cell (the "Replicate" pair is the
# same animal sequenced twice, not two animals) - so there is zero
# biological replication in any single cell, and deseq2 cannot estimate
# per-gene dispersion within one dpi/pair subset. This gave an error.
#
# workaround (standard deseq2 practice for designs without per-cell
# replication): estimate dispersion once per condition_pair using all dpi
# for that pair together (real residual df here, borrowing strength across
# timepoints), then switch to a saturated cell-means design
# (~0 + dpi_condition) and pull out each day's specific contrast by name,
# reusing the dispersion already estimated rather than re-estimating it on
# the saturated (0 residual df) design.
#
# caveat: this assumes dispersion is roughly stable across dpi within a
# condition_pair. Reasonable given the alternative (no test at all), but
# every resulting p-value is an approximation, not a per-cell-replicated
# estimate, hence we should treat hits as leads, not settled results, until 
# validated another way.

lfc_threshold  <- 1
padj_threshold <- 0.05

for (pair in condition_pairs) {
  ctrl <- pair[["control"]]
  trt  <- pair[["treatment"]]

  # Create comparison directory name (e.g., "NK_oP_vs_oP")
  comparison_dir <- paste0(trt, "_vs_", ctrl)

  # Create directories for this comparison in SECTION_6 and SECTION_7
  dir.create(file.path(outdir, "SECTION_6", comparison_dir), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "SECTION_7", comparison_dir), recursive = TRUE, showWarnings = FALSE)

  pair_samples <- rownames(collapsed_meta)[collapsed_meta$condition %in% c(ctrl, trt)]
  dds_pair <- dds[, pair_samples]
  dds_pair$condition <- droplevels(factor(dds_pair$condition, levels = c(ctrl, trt)))
  dds_pair$dpi <- factor(dds_pair$dpi)

  design(dds_pair) <- ~ dpi + condition
  dds_pair <- estimateSizeFactors(dds_pair)
  dds_pair <- estimateDispersions(dds_pair)

  dds_pair$dpi_condition <- factor(paste(dds_pair$condition, dds_pair$dpi, sep = "_dpi"))
  design(dds_pair) <- ~ 0 + dpi_condition
  dds_pair <- nbinomWaldTest(dds_pair)

  for (day in sort(unique(collapsed_meta$dpi[collapsed_meta$condition %in% c(ctrl, trt)]))) {
    trt_lvl  <- paste0(trt, "_dpi", day)
    ctrl_lvl <- paste0(ctrl, "_dpi", day)

    if (!all(c(trt_lvl, ctrl_lvl) %in% levels(dds_pair$dpi_condition))) {
      message("Skipping ", trt, "_vs_", ctrl, " at dpi ", day, " - level missing.")
      next
    }

    bmab <- sapply(c(trt_lvl, ctrl_lvl), function(lvl)
      rowMeans(counts(dds_pair, normalized = TRUE)[, dds_pair$dpi_condition == lvl, drop = FALSE]))

    res <- results(dds_pair, contrast = c("dpi_condition", trt_lvl, ctrl_lvl),
                   alpha = padj_threshold, lfcThreshold = lfc_threshold)

    z <- cbind(gene_id = rownames(res), bmab, res)
    z$padj[!is.na(z$pvalue) & z$pvalue == 1] <- 1

    z_up   <- subset(z, padj <= padj_threshold & log2FoldChange >=  lfc_threshold)
    z_down <- subset(z, padj <= padj_threshold & log2FoldChange <= -lfc_threshold)

    tag <- paste0(species, "_", trt, "_vs_", ctrl, "_dpi", day)

    # Full differential expression results with base means and statistics
    write.table(z,      file.path(outdir, "SECTION_6", comparison_dir, paste0(tag, "_DEGs.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # Significantly up-regulated genes (padj < 0.05, log2FC >= 1)
    write.table(z_up,   file.path(outdir, "SECTION_6", comparison_dir, paste0(tag, "_upDEGs.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # Significantly down-regulated genes (padj < 0.05, log2FC <= -1)
    write.table(z_down, file.path(outdir, "SECTION_6", comparison_dir, paste0(tag, "_downDEGs.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    cat("Done:", tag, "-", nrow(z_up), "up,", nrow(z_down), "down\n")
  }

  # Save DESeq2 object for this condition pair across all timepoints
  saveRDS(dds_pair, file.path(outdir, "SECTION_6", comparison_dir, paste0("dds_", species, "_", trt, "_vs_", ctrl, ".rds")))
}

cat("Section 6 (padj-gated) contrasts complete.\n")

# ---- 7. Exploratory per-timepoint differential expression --------------------
# padj here is only exploratory (see section 6's caveat) here, the primary candidate
# filtering instead uses: expression in at least one condition, effect
# size (|log2FC|), technical-replicate consistency as a QC metric, and
# temporal consistency across dpi. padj/pvalue are retained as supporting
# stats only, not the primary definition of a candidate gene.

min_mean_normalized_count <- 10
lfc_candidate_threshold   <- 1.5

# ---- 7a. Technical-replicate QC ------------------------------------------
# Absolute log2 difference between the two technical sequencing runs of the
# same biological sample, on the uncollapsed data.

dds_raw_qc <- estimateSizeFactors(dds_raw)
raw_norm_counts <- counts(dds_raw_qc, normalized = TRUE)

group_levels <- unique(meta$group)
tech_diff_list <- vector("list", length(group_levels))
names(tech_diff_list) <- group_levels

for (g in group_levels) {
  runs <- meta$sample[meta$group == g]
  if (length(runs) != 2) {
    warning("Group ", g, " has ", length(runs), " technical replicates; expected exactly 2.")
    next
  }
  tech_diff_list[[g]] <- abs(
    log2(raw_norm_counts[, runs[1]] + 1) - log2(raw_norm_counts[, runs[2]] + 1)
  )
}
tech_diff_matrix <- do.call(cbind, tech_diff_list)

# ---- 7b. Per-condition-pair DESeq2 analysis -------------------------------

for (pair in condition_pairs) {
  ctrl <- pair[["control"]]
  trt  <- pair[["treatment"]]

  # Reuse comparison directory name (e.g., "NK_oP_vs_oP")
  comparison_dir <- paste0(trt, "_vs_", ctrl)

  pair_samples <- rownames(collapsed_meta)[collapsed_meta$condition %in% c(ctrl, trt)]
  dds_pair <- dds[, pair_samples]
  dds_pair$condition <- droplevels(factor(dds_pair$condition, levels = c(ctrl, trt)))
  dds_pair$dpi <- factor(dds_pair$dpi)

  design(dds_pair) <- ~ dpi + condition
  dds_pair <- estimateSizeFactors(dds_pair)
  dds_pair <- estimateDispersions(dds_pair)

  dds_pair$dpi_condition <- factor(paste(dds_pair$condition, dds_pair$dpi, sep = "_dpi"))
  design(dds_pair) <- ~ 0 + dpi_condition
  dds_pair <- nbinomWaldTest(dds_pair)

  norm_counts_pair <- counts(dds_pair, normalized = TRUE)

  dpi_values <- sort(unique(collapsed_meta$dpi[collapsed_meta$condition %in% c(ctrl, trt)]))

  lfc_matrix <- matrix(NA_real_, nrow = nrow(dds_pair), ncol = length(dpi_values))
  rownames(lfc_matrix) <- rownames(dds_pair)
  colnames(lfc_matrix) <- paste0("dpi", dpi_values)

  all_results <- vector("list", length(dpi_values))
  names(all_results) <- paste0("dpi", dpi_values)

  # ---- 7c. Per-dpi contrasts ----------------------------------------------

  for (i in seq_along(dpi_values)) {
    day <- dpi_values[i]
    trt_lvl  <- paste0(trt, "_dpi", day)
    ctrl_lvl <- paste0(ctrl, "_dpi", day)

    if (!all(c(trt_lvl, ctrl_lvl) %in% levels(dds_pair$dpi_condition))) {
      message("Skipping ", trt, "_vs_", ctrl, " at dpi ", day, " - level missing.")
      next
    }

    ctrl_idx <- dds_pair$dpi_condition == ctrl_lvl
    trt_idx  <- dds_pair$dpi_condition == trt_lvl
    mean_ctrl <- rowMeans(norm_counts_pair[, ctrl_idx, drop = FALSE])
    mean_trt  <- rowMeans(norm_counts_pair[, trt_idx, drop = FALSE])

    res <- results(dds_pair, contrast = c("dpi_condition", trt_lvl, ctrl_lvl),
                   lfcThreshold = lfc_candidate_threshold, independentFiltering = TRUE)

    ctrl_group <- collapsed_meta$group[collapsed_meta$condition == ctrl & collapsed_meta$dpi == day]
    trt_group  <- collapsed_meta$group[collapsed_meta$condition == trt  & collapsed_meta$dpi == day]

    tech_ctrl <- if (length(ctrl_group) == 1 && ctrl_group %in% colnames(tech_diff_matrix)) {
      tech_diff_matrix[, ctrl_group]
    } else {
      rep(NA_real_, nrow(dds_pair))
    }
    tech_trt <- if (length(trt_group) == 1 && trt_group %in% colnames(tech_diff_matrix)) {
      tech_diff_matrix[, trt_group]
    } else {
      rep(NA_real_, nrow(dds_pair))
    }
    tech_max_diff <- pmax(tech_ctrl, tech_trt, na.rm = TRUE)
    tech_max_diff[is.infinite(tech_max_diff)] <- NA_real_

    z <- data.frame(
      gene_id = rownames(res),
      mean_normalized_control = mean_ctrl,
      mean_normalized_infected = mean_trt,
      res,
      technical_max_abs_log2_difference = tech_max_diff,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    z$padj[!is.na(z$pvalue) & z$pvalue == 1] <- 1

    lfc_matrix[, i] <- res$log2FoldChange

    expressed <- (
      z$mean_normalized_control >= min_mean_normalized_count |
      z$mean_normalized_infected >= min_mean_normalized_count
    )
    sufficient_effect <- (
      !is.na(z$log2FoldChange) &
      abs(z$log2FoldChange) >= lfc_candidate_threshold
    )
    candidate <- expressed & sufficient_effect
    z$candidate <- candidate

    z$technical_replicates_consistent <- (
      !is.na(z$technical_max_abs_log2_difference) &
      z$technical_max_abs_log2_difference <= 1
    )
    z$technical_replicate_status <- ifelse(
      is.na(z$technical_max_abs_log2_difference), "not_available",
      ifelse(z$technical_max_abs_log2_difference <= 1, "consistent", "inspect")
    )
    z$statistical_evidence <- "exploratory_pooled_dispersion"

    z_up   <- subset(z, candidate & log2FoldChange >=  lfc_candidate_threshold)
    z_down <- subset(z, candidate & log2FoldChange <= -lfc_candidate_threshold)

    tag <- paste0(species, "_", trt, "_vs_", ctrl, "_dpi", day)

    # Full exploratory results with statistics, fold changes, and QC metrics
    write.table(z, file.path(outdir, "SECTION_7", comparison_dir, paste0(tag, "_DESeq2_exploratory.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # Significantly up/down-regulated genes passing all filters (expression + effect size + QC)
    write.table(z_up, file.path(outdir, "SECTION_7", comparison_dir, paste0(tag, "_candidate_up.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(z_down, file.path(outdir, "SECTION_7", comparison_dir, paste0(tag, "_candidate_down.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    all_results[[paste0("dpi", day)]] <- z

    cat("Done:", tag, "- candidates:", nrow(z_up), "up,", nrow(z_down), "down\n")
  }

  # ---- 7d. Infection-effect matrix across the complete time course --------

  # Log2 fold changes for each gene across all timepoints (columns = dpi0, dpi1, etc.)
  lfc_output <- data.frame(gene_id = rownames(lfc_matrix), lfc_matrix, check.names = FALSE)
  write.table(lfc_output, file.path(outdir, "SECTION_7", comparison_dir, paste0(species, "_", trt, "_vs_", ctrl, "_infection_LFC_matrix.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- 7e. Combined candidate table across all dpi ------------------------

  combined_candidates <- bind_rows(all_results, .id = "dpi") %>%
    filter(candidate) %>%
    arrange(gene_id, desc(abs(log2FoldChange)))

  # Combined significant candidates across all timepoints, sorted by effect size
  write.table(combined_candidates, file.path(outdir, "SECTION_7", comparison_dir, paste0(species, "_", trt, "_vs_", ctrl, "_candidate_genes_all_DPI.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- 7f. Temporal summary for each gene ----------------------------------

  max_abs_lfc <- apply(abs(lfc_matrix), 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
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
  temporal_summary$passed_candidate_filter_any_DPI <- temporal_summary$gene_id %in% candidate_genes

  # Maximum absolute log2FC and corresponding timepoint for each gene across timecourse
  write.table(temporal_summary, file.path(outdir, "SECTION_7", comparison_dir, paste0(species, "_", trt, "_vs_", ctrl, "_temporal_summary.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  saveRDS(dds_pair, file.path(outdir, "SECTION_7", comparison_dir, paste0("dds_", species, "_", trt, "_vs_", ctrl, ".rds")))
}

cat("All section 7 contrasts complete. Outputs written to:", outdir, "\n")
