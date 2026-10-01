#!/usr/bin/env Rscript
# ============================================================================
# run_cross_species_comparison.R
#
# Answers "what pathways are commonly affected across cattle and sheep"
#
# Usage:
#   run_cross_species_comparison.R <cattle_enrichment_dir> <sheep_enrichment_dir> <outdir>
# ============================================================================

# ---- EDITABLE SETTINGS -----------------------------------------------------
SHARED_DPI       <- c("0", "2", "6")
TOP_N_PLOT_TERMS <- 25   # most-significant shared terms shown per plot
PLOT_WIDTH  <- 11
PLOT_HEIGHT <- 9
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_cross_species_comparison.R <cattle_enrichment_dir> <sheep_enrichment_dir> <outdir>")
}
cattle_dir <- args[1]
sheep_dir  <- args[2]
outdir     <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DOSE)
})

# ---- Path parsing (same convention as run_tables.R) ------------------------
# ORA: enrichment_dir/<comparison>/dpi<N>/ORA/<direction>/<file>.rds
# GSEA: enrichment_dir/<comparison>/dpi<N>/GSEA/<file>.rds

parse_ora_path <- function(path) {
  direction  <- basename(dirname(path))
  ora_folder <- basename(dirname(dirname(path)))
  dpi_folder <- basename(dirname(dirname(dirname(path))))
  comparison <- basename(dirname(dirname(dirname(dirname(path)))))
  list(comparison = comparison, dpi = sub("^dpi", "", dpi_folder), direction = direction)
}
parse_gsea_path <- function(path) {
  dpi_folder <- basename(dirname(dirname(path)))
  comparison <- basename(dirname(dirname(dirname(path))))
  list(comparison = comparison, dpi = sub("^dpi", "", dpi_folder))
}

strip_kegg_prefix <- function(id) sub("^[a-z]+", "", id)  # "bta04610"/"hsa04610" -> "04610"

# ---- GO ORA: read + pool cattle's direct/ortholog buckets ------------------
# For cattle only, both buckets in one _GO_results.rds get rbind'ed into one
# per-contrast table — see header note on why this is the right call here.

read_ora_go <- function(enrichment_dir, species_label) {
  # Matches both _result.rds and _results.rds
  files <- list.files(enrichment_dir, pattern = "_GO_results?\\.rds$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("_GSEA_", files)]
  
  if (length(files) == 0) return(NULL)

  rows <- lapply(files, function(f) {
    info <- parse_ora_path(f)
    
    # Skip if this DPI isn't in our target analysis list
    if (!info$dpi %in% SHARED_DPI) return(NULL)
    
    go_results <- readRDS(f)
    if (is.null(go_results)) return(NULL)
    
    # Case A: Direct object structure
    if (inherits(go_results, "enrichResult") || is.data.frame(go_results) || !is.list(go_results)) {
      df <- as.data.frame(go_results)
      if (nrow(df) == 0) return(NULL)
      df$term_type <- if ("ONTOLOGY" %in% names(df)) df$ONTOLOGY else "GO"
      df$term_id <- df$ID
      df$annotation_source <- "Direct_Object"
      df$comparison <- info$comparison; df$dpi <- info$dpi; df$direction <- info$direction
      df$species <- species_label
      return(df)
    }
    
    # Case B: Nested list structure
    bind_rows(lapply(names(go_results), function(src) {
      df <- as.data.frame(go_results[[src]])
      if (nrow(df) == 0) return(NULL)
      df$term_type <- if ("ONTOLOGY" %in% names(df)) df$ONTOLOGY else "GO"
      df$term_id <- df$ID
      df$annotation_source <- src
      df$comparison <- info$comparison; df$dpi <- info$dpi; df$direction <- info$direction
      df$species <- species_label
      df
    }))
  })
  bind_rows(rows)
}

# ---- KEGG ORA: read + strip organism prefix for comparability -------------

read_ora_kegg <- function(enrichment_dir, species_label) {
  # Matches both _result.rds and _results.rds
  files <- list.files(enrichment_dir, pattern = "_KEGG_results?\\.rds$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("_GSEA_", files)]
  
  if (length(files) == 0) return(NULL)

  rows <- lapply(files, function(f) {
    info <- parse_ora_path(f)
    if (!info$dpi %in% SHARED_DPI) return(NULL)
    
    df <- as.data.frame(readRDS(f))
    if (nrow(df) == 0) return(NULL)
    df$term_type <- "KEGG"
    df$term_id <- strip_kegg_prefix(df$ID)
    df$annotation_source <- NA_character_
    df$comparison <- info$comparison; df$dpi <- info$dpi; df$direction <- info$direction
    df$species <- species_label
    df
  })
  bind_rows(rows)
}

# ---- Same two readers for GSEA ---

read_gsea_go <- function(enrichment_dir, species_label) {
  files <- list.files(enrichment_dir, pattern = "_GSEA_GO_results\\.rds$", recursive = TRUE, full.names = TRUE)
  rows <- lapply(files, function(f) {
    info <- parse_gsea_path(f)
    if (!info$dpi %in% SHARED_DPI) return(NULL)
    go_results <- readRDS(f)
    bind_rows(lapply(names(go_results), function(src) {
      df <- as.data.frame(go_results[[src]])
      if (nrow(df) == 0) return(NULL)
      df$term_type <- if ("ONTOLOGY" %in% names(df)) df$ONTOLOGY else "GO"
      df$term_id <- df$ID
      df$annotation_source <- src
      df$comparison <- info$comparison; df$dpi <- info$dpi
      df$species <- species_label
      df
    }))
  })
  bind_rows(rows)
}

read_gsea_kegg <- function(enrichment_dir, species_label) {
  files <- list.files(enrichment_dir, pattern = "_GSEA_KEGG_result\\.rds$", recursive = TRUE, full.names = TRUE)
  rows <- lapply(files, function(f) {
    info <- parse_gsea_path(f)
    if (!info$dpi %in% SHARED_DPI) return(NULL)
    df <- as.data.frame(readRDS(f))
    if (nrow(df) == 0) return(NULL)
    df$term_type <- "KEGG"
    df$term_id <- strip_kegg_prefix(df$ID)
    df$annotation_source <- NA_character_
    df$comparison <- info$comparison; df$dpi <- info$dpi
    df$species <- species_label
    df
  })
  bind_rows(rows)
}

cat("Reading cattle ORA...\n")
cattle_ora <- bind_rows(
  read_ora_go(cattle_dir, "cattle"),
  read_ora_kegg(cattle_dir, "cattle")
)
cat("Reading sheep ORA...\n")
sheep_ora <- bind_rows(
  read_ora_go(sheep_dir, "sheep"),
  read_ora_kegg(sheep_dir, "sheep")
)
ora_long <- bind_rows(cattle_ora, sheep_ora)

cat("Reading cattle GSEA...\n")
cattle_gsea <- bind_rows(
  read_gsea_go(cattle_dir, "cattle"),
  read_gsea_kegg(cattle_dir, "cattle")
)
cat("Reading sheep GSEA...\n")
sheep_gsea <- bind_rows(
  read_gsea_go(sheep_dir, "sheep"),
  read_gsea_kegg(sheep_dir, "sheep")
)
gsea_long <- bind_rows(cattle_gsea, sheep_gsea)

ensure_required_columns <- function(df) {
  req_cols <- c("term_type", "term_id", "species", "comparison", "dpi", "Description", "p.adjust")
  if (is.null(df) || nrow(df) == 0) {
    # Create an empty dataframe with the exact expected structure
    df <- as.data.frame(matrix(0, nr = 0, nc = length(req_cols)))
    names(df) <- req_cols
    return(df)
  }
  # Add missing columns as NA if they don't exist
  for (col in req_cols) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  return(df)
}

ora_long  <- ensure_required_columns(ora_long)
gsea_long <- ensure_required_columns(gsea_long)

# ---- Long (full-detail) tables ---------------------------------------------
# One row per (term, species, comparison, dpi, direction) hit, every
# contrast that found this term, not collapsed. Kept alongside the summary
# so nothing is hidden behind the found/not-found roll-up below.

write.table(ora_long %>% dplyr::arrange(term_type, term_id, species, comparison, dpi),
            file.path(outdir, "cross_species_ora_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_long %>% dplyr::arrange(term_type, term_id, species, comparison, dpi),
            file.path(outdir, "cross_species_gsea_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote long tables (", nrow(ora_long), "ORA rows,", nrow(gsea_long), "GSEA rows ).\n")

# ---- Summary tables: found in both? ----------------------------------------

build_summary <- function(long_df) {
  if (nrow(long_df) == 0) {
    return(data.frame(
      term_id = character(), term_type = character(), Description = character(),
      cattle_found = logical(), cattle_comparisons = character(), cattle_dpis = character(), cattle_min_padj = numeric(),
      sheep_found = logical(), sheep_dpis = character(), sheep_min_padj = numeric(), found_in_both = logical()
    ))
  }
  
  long_df %>%
    group_by(term_id, term_type, Description) %>%
    summarise(
      cattle_found = any(species == "cattle" & !is.na(species)),
      cattle_comparisons = paste(sort(unique(comparison[species == "cattle" & !is.na(species)])), collapse = ","),
      cattle_dpis = paste(sort(unique(dpi[species == "cattle" & !is.na(species)])), collapse = ","),
      cattle_min_padj = if(any(species == "cattle")) suppressWarnings(min(p.adjust[species == "cattle"], na.rm = TRUE)) else NA_real_,
      sheep_found = any(species == "sheep" & !is.na(species)),
      sheep_dpis = paste(sort(unique(dpi[species == "sheep" & !is.na(species)])), collapse = ","),
      sheep_min_padj = if(any(species == "sheep")) suppressWarnings(min(p.adjust[species == "sheep"], na.rm = TRUE)) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      found_in_both = cattle_found & sheep_found,
      cattle_min_padj = ifelse(is.infinite(cattle_min_padj), NA_real_, cattle_min_padj),
      sheep_min_padj = ifelse(is.infinite(sheep_min_padj), NA_real_, sheep_min_padj)
    ) %>%
    dplyr::arrange(desc(found_in_both), pmin(cattle_min_padj, sheep_min_padj, na.rm = TRUE))
}

ora_summary <- build_summary(ora_long)
gsea_summary <- build_summary(gsea_long)

write.table(ora_summary, file.path(outdir, "cross_species_ora_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_summary, file.path(outdir, "cross_species_gsea_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote summary tables. Shared (found_in_both) terms ORA:",
    sum(ora_summary$found_in_both), ", GSEA:", sum(gsea_summary$found_in_both), "\n")

# ---- Plots: shared terms only, one dot per (species, comparison) ----------
# x = group (species + comparison, so cattle's 4 comparisons stay visible
# rather than collapsed into one "cattle" column), y = term, size = gene
# count / setSize, color = significance.

plot_shared_terms <- function(long_df, summary_df, value_col, value_label, title, out_file) {
  shared_ids <- summary_df %>% filter(found_in_both) %>% dplyr::arrange(pmin(cattle_min_padj, sheep_min_padj, na.rm = TRUE)) %>%
    slice_head(n = TOP_N_PLOT_TERMS) %>% pull(term_id)
  if (length(shared_ids) == 0) {
    message("No shared terms to plot for: ", title)
    return(invisible(NULL))
  }

  plot_df <- long_df %>%
    filter(term_id %in% shared_ids) %>%
    mutate(group_label = ifelse(species == "cattle", paste0("Cattle: ", comparison), "Sheep"))

  # ORA rows have "Count", GSEA rows have "setSize" — never both in the same
  # call (ora_long/gsea_long are passed in separately) so pick whichever is
  # actually present rather than referencing a column that doesn't exist.
  size_col <- if ("Count" %in% names(plot_df)) "Count" else "setSize"
  plot_df$size_value <- plot_df[[size_col]]

  p <- ggplot(plot_df, aes(x = group_label, y = reorder(Description, .data[[value_col]]),
                            size = size_value, color = .data[[value_col]])) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(name = value_label, low = "steelblue", high = "firebrick") +
     labs(
      title = gsub("\\s*[—-]\\s*", " ", iconv(title, to = "ASCII//TRANSLIT")), 
      x = NULL, 
      y = NULL, 
      size = "Gene count"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(out_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  cat("Wrote", out_file, "\n")
}

plot_shared_terms(ora_long, ora_summary, "p.adjust", "p.adjust",
                   "Cross-species ORA shared significant terms", file.path(outdir, "cross_species_ora_dotplot.pdf"))
plot_shared_terms(gsea_long, gsea_summary, "NES", "NES",
                   "Cross-species GSEA shared significant terms", file.path(outdir, "cross_species_gsea_dotplot.pdf"))

# ----- Directional filtering ----------------------------------

cat("Filtering directional ORA pathways...\n")

# Extract pathways exclusively Up-regulated in both species
ora_both_up <- ora_long %>%
  filter(direction == "up") %>%
  group_by(term_id, term_type, Description) %>%
  summarise(
    cattle_up = any(species == "cattle"),
    sheep_up  = any(species == "sheep"),
    cattle_dpis = paste(sort(unique(dpi[species == "cattle"])), collapse = ","),
    sheep_dpis  = paste(sort(unique(dpi[species == "sheep"])), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(cattle_up & sheep_up) %>%
  dplyr::arrange(term_type, term_id)

# Extract pathways exclusively Down-regulated in both species
ora_both_down <- ora_long %>%
  filter(direction == "down") %>%
  group_by(term_id, term_type, Description) %>%
  summarise(
    cattle_down = any(species == "cattle"),
    sheep_down  = any(species == "sheep"),
    cattle_dpis = paste(sort(unique(dpi[species == "cattle"])), collapse = ","),
    sheep_dpis  = paste(sort(unique(dpi[species == "sheep"])), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(cattle_down & sheep_down) %>%
  dplyr::arrange(term_type, term_id)

write.table(ora_both_up, file.path(outdir, "cross_species_ora_both_up.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora_both_down, file.path(outdir, "cross_species_ora_both_down.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved directional tables. Both Up:", nrow(ora_both_up), ", Both Down:", nrow(ora_both_down), "\n")


# ---- Immune response extraction --------------------------------------

cat("Extracting immune-related pathways...\n")

# Define comprehensive immune keywords (case-insensitive regex)
immune_regex <- paste0(
  "(immune|cytokine|interferon|signaling|t[ -]cell|b[ -]cell|",
  "chemokine|toll|interleukin|defense|virus|viral|inflam|",
  "leukocyte|macrophage|lymphocyte|antigen|chemotaxis|bacter|complement)"
)

# Filter the ORA shared terms for immune terms
ora_immune_shared <- ora_summary %>%
  filter(found_in_both == TRUE) %>%
  filter(grepl(immune_regex, Description, ignore.case = TRUE))

# Filter the GSEA shared terms for immune terms
gsea_immune_shared <- gsea_summary %>%
  filter(found_in_both == TRUE) %>%
  filter(grepl(immune_regex, Description, ignore.case = TRUE))

write.table(ora_immune_shared, file.path(outdir, "cross_species_ora_summary_immune.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_immune_shared, file.path(outdir, "cross_species_gsea_summary_immune.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved immune pathways. Shared ORA Immune:", nrow(ora_immune_shared), 
    ", Shared GSEA Immune:", nrow(gsea_immune_shared), "\n")

# ---- Plots: Dedicated Immune Pathway Visualization -------------

generate_immune_dotplots <- function(long_df, immune_summary_df, method_label) {
  # Isolate the top significant shared immune terms
  shared_immune_ids <- immune_summary_df %>%
    arrange(pmin(cattle_min_padj, sheep_min_padj, na.rm = TRUE)) %>%
    head(TOP_N_PLOT_TERMS) %>%
    pull(term_id)
  
  if (length(shared_immune_ids) == 0) {
    cat("No shared immune terms found to plot for", method_label, "\n")
    return(NULL)
  }
  
  # Filter the long dataframe and safely extract metric tracking
  is_gsea <- (method_label == "GSEA")
  
  plot_data <- long_df %>%
    filter(term_id %in% shared_immune_ids) %>%
    group_by(term_id, Description, species, dpi) %>%
    summarise(
      min_padj = min(p.adjust, na.rm = TRUE),
      # Calculate ORA Count or take absolute GSEA NES score depending on the runtime context
      size_metric = if(is_gsea) {
        if("NES" %in% names(.)) max(abs(NES), na.rm = TRUE) else 1.0
      } else {
        if("Count" %in% names(.)) sum(Count, na.rm = TRUE) else 1.0
      },
      .groups = "drop"
    )
  
  # Maintain strict factor levels based on significance ranking
  term_order <- immune_summary_df %>%
    filter(term_id %in% shared_immune_ids) %>%
    arrange(pmin(cattle_min_padj, sheep_min_padj, na.rm = TRUE)) %>%
    pull(Description)
  
  plot_data <- plot_data %>%
    mutate(Description = factor(Description, levels = rev(term_order)))
  
  # Define labels dynamically based on analytical method
  size_title <- if(is_gsea) "Absolute NES" else "Gene Count"
  raw_title  <- paste0("Top Shared Immune Pathways (", method_label, ")")
  clean_title <- gsub("\\s*[—-]\\s*", " ", iconv(raw_title, to = "ASCII//TRANSLIT"))

  p <- ggplot(plot_data, aes(x = factor(dpi), y = Description)) +
    geom_point(aes(size = size_metric, color = -log10(min_padj)), alpha = 0.85) +
    facet_wrap(~species, scales = "fixed") +
    scale_color_gradient(low = "steelblue", high = "firebrick") +
    labs(
      title = clean_title,
      x = "Days Post Infection (DPI)",
      y = "Immune Pathway / Term",
      color = "-log10(p.adjust)",
      size = size_title
    ) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 0, hjust = 0.5),
      strip.text   = element_text(face = "bold", size = 12),
      panel.grid.major = element_line(color = "gray92"),
      plot.title   = element_text(hjust = 0.5, face = "bold")
    )
  
  outfile <- file.path(outdir, paste0("cross_species_", tolower(method_label), "_immune_dotplot.pdf"))
  ggsave(outfile, plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  cat("Saved immune dotplot to:", outfile, "\n")
}

cat("Generating cross-species immune plots...\n")
generate_immune_dotplots(ora_long, ora_immune_shared, "ORA")
generate_immune_dotplots(gsea_long, gsea_immune_shared, "GSEA")

cat("All cross-species comparison outputs written to:", outdir, "\n")


# ---- Plots: Clean Targeted Comparison Visualization -----------------------

generate_targeted_comparison_plots <- function(long_df, summary_df, method_label, immune_only = FALSE) {
  # Determine which summary source to use for the top N ranking
  target_summary <- if (immune_only) ora_immune_shared else summary_df
  if (method_label == "GSEA" && immune_only) target_summary <- gsea_immune_shared
  
  # Isolate top shared terms
  shared_terms <- target_summary %>%
    filter(found_in_both == TRUE) %>%
    arrange(pmin(cattle_min_padj, sheep_min_padj, na.rm = TRUE)) %>%
    head(TOP_N_PLOT_TERMS) %>%
    pull(term_id)
  
  if (length(shared_terms) == 0) return(NULL)
  
  # Keep all sheep, but restrict cattle only to 'oP_vs_NK_oP'
  filtered_long <- long_df %>%
    filter(
      term_id %in% shared_terms,
      (species == "sheep") | (species == "cattle" & comparison == "oP_vs_NK_oP")
    )
  
  if (nrow(filtered_long) == 0) {
    cat("No data found for Cattle 'oP_vs_NK_oP' in top", method_label, "terms.\n")
    return(NULL)
  }
  
  # Calculate plotting metrics
  is_gsea <- (method_label == "GSEA")
  plot_data <- filtered_long %>%
    group_by(term_id, Description, species, dpi) %>%
    summarise(
      min_padj = min(p.adjust, na.rm = TRUE),
      size_metric = if(is_gsea) {
        if("NES" %in% names(.)) max(abs(NES), na.rm = TRUE) else 1.0
      } else {
        if("Count" %in% names(.)) sum(Count, na.rm = TRUE) else 1.0
      },
      .groups = "drop"
    )
  
  term_order <- target_summary %>%
    filter(term_id %in% shared_terms) %>%
    arrange(pmin(cattle_min_padj, sheep_min_padj, na.rm = TRUE)) %>%
    pull(Description)
  
  plot_data <- plot_data %>%
    mutate(Description = factor(Description, levels = rev(term_order)))

  size_title <- if(is_gsea) "Absolute NES" else "Gene Count"
  suffix     <- if(immune_only) " Immune Pathways" else " Pathways"
  
  raw_title  <- paste0("Top Shared", suffix, " across Cattle and Sheep (", method_label, ")")
  clean_title <- gsub("\\s*[—-]\\s*", " ", iconv(raw_title, to = "ASCII//TRANSLIT"))
  
  p <- ggplot(plot_data, aes(x = factor(dpi), y = Description)) +
    geom_point(aes(size = size_metric, color = -log10(min_padj)), alpha = 0.85) +
    facet_wrap(~species, scales = "fixed") +
    scale_color_gradient(low = "steelblue", high = "firebrick") +
    labs(
      title = clean_title,
      x = "Days Post Infection (DPI)",
      y = "Biological Pathway / Term",
      color = "-log10(p.adjust)",
      size = size_title
    ) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 0, hjust = 0.5),
      strip.text   = element_text(face = "bold", size = 12),
      panel.grid.major = element_line(color = "gray92"),
      plot.title   = element_text(hjust = 0.5, face = "bold")
    )
  
  file_suffix <- if(immune_only) "_targeted_immune_dotplot.pdf" else "_targeted_global_dotplot.pdf"
  outfile <- file.path(outdir, paste0("cross_species_", tolower(method_label), file_suffix))
  ggsave(outfile, plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  cat("Saved targeted plot to:", outfile, "\n")
}

cat("Generating targeted comparison plots (Cattle filtered vs All Sheep)...\n")
# Global Shared Plots
generate_targeted_comparison_plots(ora_long, ora_summary, "ORA", immune_only = FALSE)
generate_targeted_comparison_plots(gsea_long, gsea_summary, "GSEA", immune_only = FALSE)

# Immune Sub-plots
generate_targeted_comparison_plots(ora_long, ora_summary, "ORA", immune_only = TRUE)
generate_targeted_comparison_plots(gsea_long, gsea_summary, "GSEA", immune_only = TRUE)
