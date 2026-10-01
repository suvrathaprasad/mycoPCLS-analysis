#!/usr/bin/env Rscript
# ============================================================================
# run_plots_custom.R : plots and annotation for a user-supplied gene list
#
# Usage:
#   run_plots_custom.R <section7_dir> <species> <outdir> <cache_dir> <gene_list.tsv>
#
# The gene-list TSV must contain a column named "gene_id".
# It may optionally contain a column named "gene_symbol".
#
# Example:
#   gene_id
#   IFNG
#   IL1A
#   IL1B
#   ...
#
# IMPORTANT:
# - gene_id may contain Ensembl IDs or gene symbols.
# - For bovine data, symbols are resolved using the cached annotation map
#   and org.Bt.eg.db.
# - GO annotation uses the same direct/ortholog mapping logic as the main
#   enrichment pipeline.
# - No enrichment is rerun here.
#
# Outputs are written to:
#   <outdir>/custom_gene_plots/
# ============================================================================


# ============================================================================
# 1. Arguments
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: run_plots_custom.R <section7_dir> <species> <outdir> <cache_dir> <gene_list.tsv>")
}

section7_dir <- args[1]
species <- tolower(args[2])
outdir <- args[3]
cache_dir <- args[4]
gene_list_file <- args[5]

if (!species %in% c("cattle", "sheep")) {
  stop("species must be 'cattle' or 'sheep', got: ", species)
}

if (!dir.exists(section7_dir)) {
  stop("SECTION_7 directory not found: ", section7_dir)
}

if (!file.exists(gene_list_file)) {
  stop("Gene-list file not found: ", gene_list_file)
}


# ============================================================================
# 2. Display settings
# ============================================================================

SPECIES_DISPLAY <- list(
  cattle = "Bovine PCLS transcriptomes",
  sheep = "Sheep PCLS transcriptomes"
)

CONDITION_DISPLAY <- list(
  cattle = c(
    NK_oP = "NC w/o PANTA",
    NK_mP = "NC w/ PANTA",
    oP = "23MA0307 w/o PANTA",
    mP = "23MA0307 w/ PANTA"
  ),
  sheep = c(
    NK = "NC",
    infection = "23MA0307"
  )
)

COMPARISON_DISPLAY <- list(
  cattle = c(
    NK_mP_vs_NK_oP = "NC w/ PANTA vs. NC w/o PANTA",
    oP_vs_NK_oP = "M. bovis Infection w/o PANTA vs. NC w/o PANTA",
    mP_vs_NK_mP = "M. bovis Infection w/ PANTA vs. NC w/ PANTA",
    mP_vs_oP = "M. bovis Infection w/ PANTA vs. M. bovis Infection w/o PANTA"
  ),
  sheep = c(
    infection_vs_NK = "M. bovis Infection vs. NC"
  )
)

comparison_label <- function(x, species) {
  labels <- COMPARISON_DISPLAY[[species]]
  out <- unname(labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

PLOT_WIDTH <- 10
PLOT_HEIGHT <- 8

THEME <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )


# ============================================================================
# 3. Packages
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(GO.db)
  library(org.Hs.eg.db)

  if (species == "cattle") library(org.Bt.eg.db)
})


species_label <- SPECIES_DISPLAY[[species]]

custom_dir <- outdir
dir.create(custom_dir, recursive = TRUE, showWarnings = FALSE)

plots_data_dir <- file.path(custom_dir, "custom_plots_data")
dir.create(plots_data_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# 4. Read user gene list
# ============================================================================

if (!file.exists(gene_list_file)) {
  stop("Gene list file not found: ", gene_list_file)
}

gene_input <- read.delim(
  gene_list_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!nrow(gene_input)) {
  stop("Gene list file is empty: ", gene_list_file)
}

if (!"gene_id" %in% names(gene_input)) {
  stop("Gene list must contain a column named 'gene_id'.")
}

user_symbols <- unique(trimws(as.character(gene_input$gene_id)))
user_symbols <- user_symbols[!is.na(user_symbols) & nzchar(user_symbols)]

if (!length(user_symbols)) {
  stop("No valid gene symbols found in ", gene_list_file)
}

cat("Read", length(user_symbols), "user-provided gene symbols.\n")

# ============================================================================
# 5. Annotation cache
# ============================================================================

annot_map_file <- file.path(
  cache_dir,
  paste0(species, "_annotation_map.rds")
)

if (!file.exists(annot_map_file)) {
  stop("Annotation cache not found: ", annot_map_file)
}

annot_map <- readRDS(annot_map_file)

required_annotation_columns <- c(
  "gene_id",
  "entrez_go",
  "source_go",
  "entrez_kegg",
  "source_kegg"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(annot_map)
)

if (length(missing_annotation_columns)) {
  stop(
    "Annotation map is missing required columns: ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}


# ============================================================================
# 6. Gene symbol / ID resolution
# ============================================================================

direct_ids <- unique(
  na.omit(
    as.character(
      annot_map$entrez_go[annot_map$source_go == "direct"]
    )
  )
)

ortholog_ids <- unique(
  na.omit(
    as.character(
      annot_map$entrez_go[annot_map$source_go == "human_ortholog"]
    )
  )
)

sym_direct <- if (species == "cattle" && length(direct_ids) > 0) {
  AnnotationDbi::select(
    org.Bt.eg.db,
    keys = direct_ids,
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
} else {
  data.frame(
    ENTREZID = character(0),
    SYMBOL = character(0)
  )
}

sym_ortholog <- if (length(ortholog_ids) > 0) {
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ortholog_ids,
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
} else {
  data.frame(
    ENTREZID = character(0),
    SYMBOL = character(0)
  )
}

sym_lookup <- bind_rows(sym_direct, sym_ortholog) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENTREZID, .keep_all = TRUE)

symbol_for_entrez <- setNames(
  sym_lookup$SYMBOL,
  as.character(sym_lookup$ENTREZID)
)

resolved_symbol <- symbol_for_entrez[
  as.character(annot_map$entrez_go)
]

gene_labels <- setNames(
  ifelse(
    !is.na(resolved_symbol) & resolved_symbol != "",
    resolved_symbol,
    annot_map$gene_id
  ),
  annot_map$gene_id
)

label_genes <- function(ids) {
  out <- unname(gene_labels[ids])
  out[is.na(out)] <- ids[is.na(out)]
  out
}

# ============================================================================
# 7. Resolve user genes against annotation map
# ============================================================================

user_genes <- unique(
  trimws(as.character(gene_input$gene_id))
)

user_genes <- user_genes[
  !is.na(user_genes) & nzchar(user_genes)
]

# Directly supplied annotation IDs are accepted if present.

direct_matches <- annot_map %>%
  filter(gene_id %in% user_genes) %>%
  distinct(gene_id, .keep_all = TRUE)

# Resolve supplied gene symbols through the labels generated in Section 6.

symbol_map <- data.frame(
  gene_id = annot_map$gene_id,
  gene_symbol = unname(gene_labels[annot_map$gene_id]),
  stringsAsFactors = FALSE
) %>%
  filter(
    !is.na(gene_symbol),
    nzchar(gene_symbol)
  ) %>%
  distinct(gene_symbol, .keep_all = TRUE)

symbol_matches <- symbol_map %>%
  filter(gene_symbol %in% user_genes)

annotation_match <- bind_rows(
  direct_matches,
  annot_map %>%
    filter(gene_id %in% symbol_matches$gene_id) %>%
    distinct(gene_id, .keep_all = TRUE)
) %>%
  distinct(gene_id, .keep_all = TRUE)

# Build a requested-symbol -> annotation-ID lookup.

resolved_lookup <- symbol_matches %>%
  select(
    requested_gene = gene_symbol,
    gene_id
  )

resolved_input <- data.frame(
  requested_gene = user_genes,
  stringsAsFactors = FALSE
) %>%
  left_join(
    resolved_lookup,
    by = "requested_gene"
  )

# Preserve directly supplied annotation IDs.

resolved_input$gene_id[
  is.na(resolved_input$gene_id) &
    resolved_input$requested_gene %in% direct_matches$gene_id
] <- resolved_input$requested_gene[
  is.na(resolved_input$gene_id) &
    resolved_input$requested_gene %in% direct_matches$gene_id
]

resolved_input <- resolved_input %>%
  left_join(
    annot_map %>%
      select(
        gene_id,
        entrez_go,
        source_go,
        entrez_kegg,
        source_kegg
      ) %>%
      distinct(gene_id, .keep_all = TRUE),
    by = "gene_id"
  )

resolved_input$gene_symbol <- ifelse(
  !is.na(resolved_input$gene_id),
  label_genes(resolved_input$gene_id),
  resolved_input$requested_gene
)

write.table(
  resolved_input,
  file.path(
    custom_dir,
    paste0("custom_gene_annotation_", species, ".tsv")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

n_resolved <- sum(
  !is.na(resolved_input$gene_id) &
    nzchar(resolved_input$gene_id)
)

cat(
  "Resolved",
  n_resolved,
  "of",
  length(user_genes),
  "user-provided genes.\n"
)

if (n_resolved < length(user_genes)) {
  unresolved <- resolved_input$requested_gene[
    is.na(resolved_input$gene_id)
  ]

  cat(
    "Unresolved genes:",
    paste(unresolved, collapse = ", "),
    "\n"
  )
}

# ============================================================================
# 8. Discover DESeq2 exploratory files
# ============================================================================

expl_files <- list.files(
  section7_dir,
  pattern = "_DESeq2_exploratory\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (!length(expl_files)) {
  stop("No DESeq2 exploratory TSV files found under ", section7_dir)
}

parse_exploratory <- function(path) {
  comparison <- basename(dirname(path))
  fname <- basename(path)

  m <- regmatches(
    fname,
    regexec(
      "_dpi([^_]+)_DESeq2_exploratory\\.tsv$",
      fname
    )
  )[[1]]

  if (length(m) != 2) {
    return(NULL)
  }

  data.frame(
    comparison = comparison,
    dpi = m[2],
    file = path,
    stringsAsFactors = FALSE
  )
}

expl_index <- bind_rows(
  lapply(expl_files, parse_exploratory)
)

expl_index$dpi_num <- suppressWarnings(
  as.numeric(expl_index$dpi)
)

expl_index <- expl_index %>%
  arrange(comparison, dpi_num)

cat(
  "Found",
  nrow(expl_index),
  "comparison/dpi exploratory files.\n"
)


# ============================================================================
# 9. Collect custom-gene statistics
# ============================================================================

statistics_list <- list()

for (i in seq_len(nrow(expl_index))) {

  info <- expl_index[i, ]
  df <- read.delim(
    info$file,
    stringsAsFactors = FALSE
  )

  if (!"gene_id" %in% names(df)) next

  hits <- df %>%
    filter(gene_id %in% resolved_input$gene_id)

  if (!nrow(hits)) next

  hits <- hits %>%
    mutate(
      comparison = info$comparison,
      dpi = info$dpi,
      dpi_num = info$dpi_num,
      gene_symbol = label_genes(gene_id)
    )

  statistics_list[[length(statistics_list) + 1]] <- hits
}

if (!length(statistics_list)) {
  stop(
    "None of the supplied genes were found in the DESeq2 exploratory files."
  )
}

custom_stats <- bind_rows(statistics_list)

# Put useful columns first when they exist.

preferred_columns <- c(
  "gene_id",
  "gene_symbol",
  "comparison",
  "dpi",
  "dpi_num",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "candidate"
)

preferred_columns <- intersect(
  preferred_columns,
  names(custom_stats)
)

custom_stats <- custom_stats %>%
  select(
    all_of(preferred_columns),
    everything()
  ) %>%
  arrange(
    gene_symbol,
    comparison,
    dpi_num
  )

write.table(
  custom_stats,
  file.path(
    custom_dir,
    paste0("custom_gene_statistics_", species, ".tsv")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  custom_stats,
  file.path(
    plots_data_dir,
    paste0("custom_gene_statistics_", species, ".rds")
  )
)

cat(
  "Collected statistics for",
  n_distinct(custom_stats$gene_id),
  "of the supplied genes.\n"
)


# ============================================================================
# 10. Summary table
# ============================================================================

summary_cols <- intersect(
  c(
    "gene_id",
    "gene_symbol",
    "comparison",
    "dpi",
    "dpi_num",
    "log2FoldChange",
    "pvalue",
    "padj",
    "candidate"
  ),
  names(custom_stats)
)

custom_summary <- custom_stats %>%
  select(all_of(summary_cols)) %>%
  arrange(gene_symbol, comparison, dpi_num)

write.table(
  custom_summary,
  file.path(
    custom_dir,
    paste0("custom_gene_summary_", species, ".tsv")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================================
# 11. Expression heatmap
# ============================================================================

dds_files <- list.files(
  section7_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

dds_files <- dds_files[
  basename(dirname(dds_files)) %in%
    names(COMPARISON_DISPLAY[[species]])
]

if (length(dds_files)) {

  for (dds_file in dds_files) {

    comparison <- basename(dirname(dds_file))

    cat(
      "Processing expression heatmap for",
      comparison,
      "\n"
    )

    dds <- readRDS(dds_file)

    if (!inherits(dds, "DESeqDataSet")) {
      message(
        "Skipping expression heatmap for ",
        comparison,
        " — RDS is not a DESeqDataSet."
      )
      next
    }

    # Ensure size factors are available
    if (is.null(sizeFactors(dds))) {
      dds <- DESeq2::estimateSizeFactors(dds)
    }

    norm_counts <- counts(
      dds,
      normalized = TRUE
    )

    genes_present <- intersect(
      resolved_input$gene_id,
      rownames(norm_counts)
    )

    if (length(genes_present) < 2) {
      message(
        "Skipping expression heatmap for ",
        comparison,
        " fewer than 2 supplied genes were found in normalized counts."
      )
      next
    }

    mat <- log2(
      norm_counts[
        genes_present,
        ,
        drop = FALSE
      ] + 1
    )

    mat_z <- t(
      scale(t(mat))
    )

    keep <- apply(
      mat_z,
      1,
      function(x) all(is.finite(x))
    )

    mat_z <- mat_z[
      keep,
      ,
      drop = FALSE
    ]

    if (nrow(mat_z) < 2) {
      message(
        "Skipping expression heatmap for ",
        comparison,
        " fewer than 2 genes have finite z-scores."
      )
      next
    }

    rownames(mat_z) <- make.unique(
      label_genes(rownames(mat_z))
    )

    group_values <- as.character(
      colData(dds)$condition
    )

    group_labels <- CONDITION_DISPLAY[[species]][
      group_values
    ]

    group_labels[
      is.na(group_labels)
    ] <- group_values[
      is.na(group_labels)
    ]

    dpi_values <- as.character(
      colData(dds)$dpi
    )

    # Initial temporary dataframe to compute the sorting layout
    ann_col <- data.frame(
      group = group_labels,
      dpi = dpi_values,
      dpi_num = as.numeric(dpi_values),
      row.names = colnames(dds),
      stringsAsFactors = FALSE
    )

    # Order sample identifiers by experimental group first, then by ascending numeric DPI
    sample_order <- rownames(ann_col[order(ann_col$group, ann_col$dpi_num), ])
    
    # Re-order the matrix data and metadata layout synchronously
    mat_z   <- mat_z[, sample_order, drop = FALSE]
    ann_col <- ann_col[sample_order, , drop = FALSE]
    
    # Drop the temporary numerical helper column
    ann_col$dpi_num <- NULL

    # Isolate uniquely present timepoints and enforce sorted numeric factor order
    present_dpis_numeric <- sort(unique(as.numeric(ann_col$dpi)))
    
    ann_col$group <- factor(ann_col$group)
    ann_col$dpi   <- factor(ann_col$dpi, levels = as.character(present_dpis_numeric))

    group_levels <- levels(ann_col$group)
    dpi_levels   <- levels(ann_col$dpi)

    annotation_colors <- list(
      group = setNames(
        grDevices::hcl.colors(
          length(group_levels),
          "Dark 3"
        ),
        group_levels
      ),
      dpi = setNames(
        grDevices::hcl.colors(
          length(dpi_levels),
          "Viridis"
        ),
        dpi_levels
      )
    )

    outfile <- file.path(
      custom_dir,
      paste0(
        "custom_gene_expression_heatmap_",
        species,
        "_",
        comparison,
        ".pdf"
      )
    )

    datafile <- file.path(
      plots_data_dir,
      paste0(
        "custom_gene_expression_heatmap_",
        species,
        "_",
        comparison,
        ".rds"
      )
    )

    pdf(
      outfile,
      width = 12,
      height = max(
        7,
        0.25 * nrow(mat_z) + 5
      )
    )

    pheatmap(
      mat_z,
      cluster_cols = FALSE,               
      annotation_col = ann_col,
      annotation_colors = annotation_colors,
      main = paste(
        species_label,
        comparison_label(comparison, species),
        "user-selected genes"
      ),
      show_rownames = TRUE,
      angle_col = 45,
      fontsize_row = 8,
      fontsize_col = 6,
      fontsize = 8,
      annotation_names_col = TRUE,
      annotation_legend = TRUE
    )

    dev.off()

    saveRDS(
      list(
        mat = mat,
        mat_z = mat_z,
        annotation_col = ann_col,
        comparison = comparison
      ),
      datafile
    )

    cat(
      "Wrote custom-gene expression heatmap for",
      comparison,
      "\n"
    )
  }

} else {

  message(
    "Skipping expression heatmaps no comparison-specific DESeq2 RDS files found."
  )
}

# ============================================================================
# 12. Per-comparison trajectory plots
# ============================================================================

for (cd in unique(custom_stats$comparison)) {

  cd_label <- comparison_label(cd, species)

  traj <- custom_stats %>%
    filter(comparison == cd) %>%
    filter(
      !is.na(log2FoldChange),
      is.finite(log2FoldChange)
    ) %>%
    mutate(
      dpi_num = as.numeric(dpi_num),
      gene_label = gene_symbol
    )

  if (!nrow(traj)) next

  gene_levels <- unique(traj$gene_label)
  n_genes <- length(gene_levels)

 base_colors <- c(
    "#0072B2",  # blue
    "#E41A1C",  # red
    "#009E73",  # green
    "#984EA3",  # purple
    "#FF7F00",  # orange
    "#00BFC4",  # cyan
    "#F781BF",  # pink
    "#A65628",  # brown
    "#B3B300",  # olive/yellow-green
    "#6A3D9A",  # violet
    "#1B9E77",  # teal
    "#E7298A",  # magenta
    "#66A61E",  # yellow-green
    "#E6AB02",  # mustard
    "#7570B3",  # blue-violet
    "#D95F02",  # burnt orange
    "#1F78B4",  # medium blue
    "#B15928",  # sienna
    "#8DD3C7",  # turquoise
    "#FB8072"   # salmon
  )

  if (n_genes <= length(base_colors)) {
    gene_colors <- base_colors[seq_len(n_genes)]
  } else {
    # If there are more genes than the predefined palette,
    # generate additional well-separated colours.
    gene_colors <- grDevices::hcl.colors(
      n_genes,
      palette = "Dynamic"
    )
  }

  names(gene_colors) <- gene_levels

  p <- ggplot(
    traj,
    aes(
      x = dpi_num,
      y = log2FoldChange,
      color = gene_label,
      group = gene_id
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_line(
      linewidth = 0.9,
      alpha = 0.9
    ) +
    geom_point(
      size = 2.5,
      alpha = 1
    ) +
    scale_color_manual(
      values = gene_colors,
      breaks = gene_levels
    ) +
    scale_x_continuous(
      breaks = sort(unique(traj$dpi_num))
    ) +
    labs(
      title = paste(
        species_label,
        cd_label,
        "user-selected gene trajectories"
      ),
      x = "dpi",
      y = "log2 fold change",
      color = "Gene"
    ) +
    THEME +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(0.45, "cm")
    )

  outfile <- file.path(
    custom_dir,
    paste0(
      "custom_gene_trajectories_",
      species,
      "_",
      cd,
      ".pdf"
    )
  )

  ggsave(
    outfile,
    p,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT
  )

  saveRDS(
    traj,
    file.path(
      plots_data_dir,
      paste0(
        "custom_gene_trajectories_",
        species,
        "_",
        cd,
        ".rds"
      )
    )
  )

  cat(
    "Wrote custom-gene trajectory plot for",
    cd,
    "\n"
  )
}

# ============================================================================
# 13. Significance overview
# ============================================================================

if (all(c("log2FoldChange", "padj") %in% names(custom_stats))) {

  dpi_levels <- sort(unique(custom_stats$dpi_num))

  sig_plot_data <- custom_stats %>%
    left_join(
      resolved_input %>%
        select(
          gene_id,
          source_go
        ) %>%
        distinct(gene_id, .keep_all = TRUE),
      by = "gene_id"
    ) %>%
    mutate(
      comparison_label = comparison_label(comparison, species),
      gene_label = label_genes(gene_id),

      source_go = case_when(
        source_go == "direct" ~ "Direct",
        source_go == "human_ortholog" ~ "Human ortholog",
        TRUE ~ "No GO source"
      ),

      significance = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Significant up",
        !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Significant down",
        TRUE ~ "Not significant"
      ),

      dpi_factor = factor(
        dpi_num,
        levels = dpi_levels
      )
    )

  p_sig <- ggplot(
    sig_plot_data,
    aes(
      x = dpi_factor,
      y = gene_label,
      color = source_go,
      shape = significance,
      size = -log10(pmax(padj, .Machine$double.xmin))
    )
  ) +
    geom_point(alpha = 0.85) +

    facet_wrap(
      ~ comparison_label,
      scales = "free_x"
    ) +

    scale_color_manual(
      values = c(
        "Direct" = "#377EB8",
        "Human ortholog" = "#E41A1C",
        "No GO source" = "grey70"
      )
    ) +

    scale_shape_manual(
      values = c(
        "Significant up" = 16,
        "Significant down" = 17,
        "Not significant" = 1
      )
    ) +

    labs(
      title = paste(
        species_label,
        "user-selected gene differential expression"
      ),
      x = "dpi",
      y = "Gene",
      color = "GO annotation source",
      shape = "Significance",
      size = "-log10(padj)"
    ) +

    THEME

  ggsave(
    file.path(
      custom_dir,
      paste0(
        "custom_gene_significance_",
        species,
        ".pdf"
      )
    ),
    p_sig,
    width = 11,
    height = max(
      7,
      0.3 * n_distinct(sig_plot_data$gene_id) + 5
    )
  )

  saveRDS(
    sig_plot_data,
    file.path(
      plots_data_dir,
      paste0(
        "custom_gene_significance_",
        species,
        ".rds"
      )
    )
  )

  cat("Wrote custom-gene significance overview.\n")
}

# ============================================================================
# 14. GO annotations for supplied genes
# ============================================================================

go_annotation <- resolved_input %>%
  filter(
    !is.na(entrez_go),
    !is.na(source_go)
  ) %>%
  mutate(
    entrez_go = as.character(entrez_go)
  )

if (nrow(go_annotation)) {

  go_rows <- list()

  for (src in unique(go_annotation$source_go)) {

    ids <- unique(
      go_annotation$entrez_go[
        go_annotation$source_go == src
      ]
    )

    if (!length(ids)) next

    orgdb <- if (
      species == "cattle" && src == "direct"
    ) {
      org.Bt.eg.db
    } else {
      org.Hs.eg.db
    }

    go_map <- tryCatch(
      AnnotationDbi::select(
        orgdb,
        keys = ids,
        keytype = "ENTREZID",
        columns = c(
          "SYMBOL",
          "GENENAME",
          "GO",
          "ONTOLOGY"
        )
      ),
      error = function(e) {
        message(
          "GO annotation failed for source ",
          src,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(go_map) || !nrow(go_map)) next

    go_map$source <- src

    go_rows[[length(go_rows) + 1]] <- go_map
  }

  if (length(go_rows)) {

    go_annotation_table <- bind_rows(go_rows) %>%
      left_join(
        go_annotation %>%
          select(
            gene_id,
            entrez_go,
            source_go
          ) %>%
          distinct(),
        by = c(
          "ENTREZID" = "entrez_go",
          "source" = "source_go"
        )
      ) %>%
      mutate(
        gene_symbol = label_genes(gene_id)
      ) %>%
      select(
        gene_id,
        gene_symbol,
        everything()
      )

    write.table(
      go_annotation_table,
      file.path(
        custom_dir,
        paste0(
          "custom_gene_GO_terms_",
          species,
          ".tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    saveRDS(
      go_annotation_table,
      file.path(
        plots_data_dir,
        paste0(
          "custom_gene_GO_terms_",
          species,
          ".rds"
        )
      )
    )

    cat(
      "Wrote GO annotations for supplied genes.\n"
    )

  } else {
    message(
      "No GO annotations could be resolved for the supplied genes."
    )
  }

} else {
  message(
    "No supplied genes have GO annotation IDs."
  )
}


# ============================================================================
# 15. GO term summary
# ============================================================================

go_file <- file.path(
  custom_dir,
  paste0(
    "custom_gene_GO_terms_",
    species,
    ".tsv"
  )
)

if (file.exists(go_file)) {

  go_terms <- read.delim(
    go_file,
    stringsAsFactors = FALSE
  )

  if ("GO" %in% names(go_terms)) {

    go_summary <- go_terms %>%
      filter(
        !is.na(GO),
        GO != ""
      ) %>%
      count(
        GO,
        ONTOLOGY,
        name = "n_genes"
      ) %>%
      arrange(
        desc(n_genes),
        GO
      )

    write.table(
      go_summary,
      file.path(
        custom_dir,
        paste0(
          "custom_gene_GO_term_summary_",
          species,
          ".tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}


# ============================================================================
# 16. Final report
# ============================================================================

#cat("\n")
#cat("============================================================\n")
#cat("Custom gene analysis complete\n")
#cat("Species: ", species_label, "\n", sep = "")
#cat("Input genes: ", nrow(gene_input), "\n", sep = "")
#cat(
#  "Genes found in DESeq2 data: ",
#  n_distinct(custom_stats$gene_id),
#  "\n",
#  sep = ""
#)
#cat("Output directory: ", custom_dir, "\n", sep = "")
#cat("============================================================\n")

# ============================================================================
# 17. Individual sample expression plots
# ============================================================================

dds_files <- list.files(
  section7_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

# Keep only RDS files located directly inside a comparison directory
dds_files <- dds_files[
  basename(dirname(dds_files)) %in%
    names(COMPARISON_DISPLAY[[species]])
]

if (length(dds_files)) {

  for (dds_file in dds_files) {

    comparison <- basename(dirname(dds_file))

    cat(
      "Processing individual sample expression plots for",
      comparison,
      "\n"
    )

    dds_custom <- readRDS(dds_file)

    if (!inherits(dds_custom, "DESeqDataSet")) {
      message(
        "Skipping individual sample expression plots for ",
        comparison,
        " — RDS is not a DESeqDataSet."
      )
      next
    }

    if (is.null(sizeFactors(dds_custom))) {
      dds_custom <- DESeq2::estimateSizeFactors(
        dds_custom
      )
    }

    norm_counts <- counts(
      dds_custom,
      normalized = TRUE
    )

    genes_present <- intersect(
      resolved_input$gene_id,
      rownames(norm_counts)
    )

    if (!length(genes_present)) {
      message(
        "Skipping individual sample expression plots for ",
        comparison,
        " — no supplied genes found."
      )
      next
    }

    expression_long <- as.data.frame(
      norm_counts[
        genes_present,
        ,
        drop = FALSE
      ]
    )

    expression_long$gene_id <- rownames(
      expression_long
    )

    expression_long <- expression_long %>%
      pivot_longer(
        cols = -gene_id,
        names_to = "sample",
        values_to = "normalized_expression"
      )

    sample_meta <- as.data.frame(
      colData(dds_custom)
    )

    sample_meta$sample <- rownames(
      sample_meta
    )

    expression_long <- expression_long %>%
      left_join(
        sample_meta,
        by = "sample"
      ) %>%
      mutate(
        gene_symbol = label_genes(gene_id),

        group_label = CONDITION_DISPLAY[[species]][
          as.character(condition)
        ],

        group_label = ifelse(
          is.na(group_label),
          as.character(condition),
          group_label
        ),

        dpi_num = suppressWarnings(
          as.numeric(as.character(dpi))
        )
      )

    saveRDS(
      expression_long,
      file.path(
        plots_data_dir,
        paste0(
          "17_individual_sample_expression_",
          species,
          "_",
          comparison,
          ".rds"
        )
      )
    )

    p_expression <- ggplot(
      expression_long,
      aes(
        x = factor(dpi_num),
        y = normalized_expression,
        color = group_label
      )
    ) +

      geom_point(
        position = position_jitter(
          width = 0.12,
          height = 0
        ),
        size = 2.2,
        alpha = 0.8
      ) +

      stat_summary(
        aes(
          group = group_label
        ),
        fun = mean,
        geom = "line",
        linewidth = 0.8,
        position = position_dodge(
          width = 0.35
        )
      ) +

      facet_wrap(
        ~ gene_symbol,
        scales = "free_y"
      ) +

      scale_y_continuous(
        trans = "log1p"
      ) +

      labs(
        title = paste(
          species_label,
          comparison_label(comparison, species),
          "expression of user-selected genes"
        ),
        x = "dpi",
        y = "Normalized expression",
        color = "Condition"
      ) +

      THEME

    outfile <- file.path(
      custom_dir,
      paste0(
        "custom_gene_individual_expression_",
        species,
        "_",
        comparison,
        ".pdf"
      )
    )

    ggsave(
      outfile,
      p_expression,
      width = 13,
      height = max(
        8,
        2.8 * ceiling(
          n_distinct(
            expression_long$gene_id
          ) / 3
        )
      )
    )

    cat(
      "Wrote individual sample expression plots for",
      comparison,
      "\n"
    )
  }

} else {

  message(
    "Skipping individual sample expression plots no comparison-specific DESeq2 RDS files found."
  )
}

# ============================================================================
# 18. Custom-gene log2FC response heatmap
# ============================================================================

if (nrow(custom_stats)) {

  custom_heatmap_data <- custom_stats %>%
    filter(
      !is.na(log2FoldChange),
      is.finite(log2FoldChange)
    ) %>%
    mutate(
      comparison_label = comparison_label(
        comparison,
        species
      ),
      column_label = paste(
        comparison_label,
        "dpi",
        dpi,
        sep = " | "
      )
    ) %>%
    select(
      gene_id,
      gene_symbol,
      comparison,
      dpi_num,
      column_label,
      log2FoldChange
    ) %>%
    distinct(
      gene_id,
      comparison,
      dpi_num,
      .keep_all = TRUE
    )

  if (nrow(custom_heatmap_data)) {

    custom_matrix <- custom_heatmap_data %>%
      select(
        gene_id,
        column_label,
        log2FoldChange
      ) %>%
      pivot_wider(
        names_from = column_label,
        values_from = log2FoldChange
      ) %>%
      as.data.frame()

    gene_names <- custom_matrix$gene_id

    custom_matrix$gene_id <- NULL

    custom_matrix <- as.matrix(
      custom_matrix
    )

    rownames(custom_matrix) <- make.unique(
      label_genes(gene_names)
    )

    # Order genes by maximum absolute response.
    row_order <- order(
      apply(
        abs(custom_matrix),
        1,
        max,
        na.rm = TRUE
      ),
      decreasing = TRUE
    )

    custom_matrix <- custom_matrix[
      row_order,
      ,
      drop = FALSE
    ]

    max_abs <- max(
      abs(custom_matrix),
      na.rm = TRUE
    )

    pdf(
      file.path(
        custom_dir,
        paste0(
          "custom_gene_log2FC_heatmap_",
          species,
          ".pdf"
        )
      ),
      width = 13,
      height = max(
        7,
        0.35 * nrow(custom_matrix) + 5
      )
    )

    pheatmap(
      custom_matrix,
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
        -max_abs,
        max_abs,
        length.out = 101
      ),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = NA,
      main = paste(
        species_label,
        "user-selected gene log2 fold changes"
      ),
      fontsize_row = 8,
      fontsize_col = 8,
    )

    dev.off()

    saveRDS(
      custom_matrix,
      file.path(
        plots_data_dir,
        paste0(
          "custom_gene_log2FC_heatmap_",
          species,
          ".rds"
        )
      )
    )

    cat(
      "Wrote custom log2FC response heatmap.\n"
    )
  }
}

# ============================================================================
# 19. GO annotation barplot for supplied genes
# this shows only the biological functions for the genes in the TSV
# not if these GO processes are significantly enriched in the genes in the TSV
# raw annotation frequency among the genes in the TSV
# ============================================================================

if (file.exists(go_file)) {

  go_terms <- read.delim(
    go_file,
    stringsAsFactors = FALSE
  )

  if (all(
    c("GO", "ONTOLOGY") %in% names(go_terms)
  )) {

    go_ids <- unique(
      go_terms$GO[
        !is.na(go_terms$GO) &
        go_terms$GO != ""
      ]
    )

    go_names <- tryCatch(
      AnnotationDbi::select(
        GO.db,
        keys = go_ids,
        keytype = "GOID",
        columns = c("TERM", "ONTOLOGY")
      ),
      error = function(e) {
        message(
          "Could not retrieve GO term names: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(go_names) && nrow(go_names)) {

      go_names <- go_names %>%
        select(
          GOID,
          TERM
        ) %>%
        distinct(
          GOID,
          .keep_all = TRUE
        )

      go_terms <- go_terms %>%
        left_join(
          go_names,
          by = c("GO" = "GOID")
        )

    } else {

      go_terms$TERM <- NA_character_
    }

    go_summary_plot <- go_terms %>%
      filter(
        !is.na(GO),
        GO != ""
      ) %>%
      group_by(
        GO,
        TERM,
        ONTOLOGY
      ) %>%
      summarise(
        n_genes = n_distinct(gene_symbol),
        gene_symbols = paste(
          sort(
            unique(
              gene_symbol[
                !is.na(gene_symbol) &
                gene_symbol != ""
              ]
            )
          ),
          collapse = ", "
        ),
        .groups = "drop"
      ) %>%
      mutate(
        # Fall back to GO ID if a term name could not be retrieved
        term_label = ifelse(
          !is.na(TERM) & TERM != "",
          TERM,
          GO
        )
      ) %>%
      arrange(
        desc(n_genes),
        term_label
      ) %>%
      slice_head(n = 30)

    write.table(
      go_summary_plot %>%
        select(
          GO,
          TERM,
          ONTOLOGY,
          n_genes,
          gene_symbols
        ),
      file.path(
        custom_dir,
        paste0(
          "custom_gene_GO_annotation_summary_",
          species,
          ".tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    saveRDS(
      go_summary_plot,
      file.path(
        plots_data_dir,
        paste0(
          "custom_gene_GO_annotation_summary_",
          species,
          ".rds"
        )
      )
    )

    if (nrow(go_summary_plot)) {

      p_go_custom <- ggplot(
        go_summary_plot,
        aes(
          x = n_genes,
          y = reorder(term_label, n_genes),
          fill = ONTOLOGY
        )
      ) +
        geom_col(
          width = 0.7
        ) +
        labs(
          title = paste(
            species_label,
            "GO annotations among user-selected genes"
          ),
          x = "Number of selected genes",
          y = "GO term",
          fill = "Ontology"
        ) +
        THEME +
        theme(
          axis.text.y = element_text(
            size = 9
          )
        )

      ggsave(
        file.path(
          custom_dir,
          paste0(
            "custom_gene_GO_annotation_barplot_",
            species,
            ".pdf"
          )
        ),
        p_go_custom,
        width = 11,
        height = max(
          7,
          0.30 * nrow(go_summary_plot) + 4
        )
      )

      cat(
        "Wrote custom-gene GO annotation barplot.\n"
      )

    } else {

      message(
        "No GO annotations available for plotting."
      )
    }
  }
}

# ============================================================================
# 20. Direct vs human-ortholog annotation provenance
# ============================================================================

provenance_data <- resolved_input %>%
  filter(
    !is.na(gene_id),
    !is.na(source_go)
  ) %>%
  mutate(
    source_label = case_when(
      source_go == "direct" ~ "Bovine direct",
      source_go == "human_ortholog" ~ "Human ortholog",
      TRUE ~ source_go
    ),
    gene_symbol = label_genes(gene_id)
  )

if (nrow(provenance_data)) {

  provenance_summary <- provenance_data %>%
    count(
      source_label,
      name = "n_genes"
    )

  p_provenance <- ggplot(
    provenance_summary,
    aes(
      x = source_label,
      y = n_genes,
      fill = source_label
    )
  ) +
    geom_col(
      width = 0.65,
      show.legend = FALSE
    ) +
    geom_text(
      aes(label = n_genes),
      vjust = -0.3,
      size = 4
    ) +
    scale_fill_manual(
      values = c(
        "Bovine direct" = "#0072B2",
        "Human ortholog" = "#E69F00"
      )
    ) +
    labs(
      title = paste(
        species_label,
        "GO annotation provenance"
      ),
      x = NULL,
      y = "Number of selected genes"
    ) +
    THEME

  ggsave(
    file.path(
      custom_dir,
      paste0(
        "custom_gene_annotation_provenance_",
        species,
        ".pdf"
      )
    ),
    p_provenance,
    width = 7,
    height = 6
  )

  # Also save the actual gene-level provenance.
  write.table(
    provenance_data %>%
      select(
        requested_gene,
        gene_id,
        gene_symbol,
        entrez_go,
        source_go
      ),
    file.path(
      custom_dir,
      paste0(
        "custom_gene_annotation_provenance_",
        species,
        ".tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(
    "Wrote annotation provenance plot.\n"
  )
}
