#!/usr/bin/env Rscript
# ============================================================================
# fetch_annotation_cache.R
#
# RUN THIS ON THE LOGIN NODE, MANUALLY, BEFORE SUBMITTING THE PIPELINE TO
# SLURM. NOT a Nextflow process — Slurm compute nodes on this HPC have no
# internet access, only the login node does, and this script is the only
# part of Phase 4 that needs a live network connection.
#
# It pre-fetches every network-dependent thing run_enrichment.R needs:
#   1. WHOLE-GENOME BioMart homology (source species -> human) + org.Hs.eg.db
#      Entrez mapping: every gene in the species' Ensembl gene set (release 116
#      archive), not just genes that happened to appear in a DESeq2 result. Built
#      once per species, so no further BioMart calls are needed for any later
#      analysis on that species, and it can be built before any DE data exists.
#      All three species (cattle, sheep, goat) use the SAME unconditional
#      human-ortholog path. org.Bt.eg.db is no longer used here — the legacy
#      cattle hybrid (direct + ortholog fallback) code is kept below as
#      commented-out blocks, not deleted.
#   2. KEGG pathway <-> gene mappings for organism "hsa" (all three species go
#      through human orthologs, so all three need the human pathway table).
#      clusterProfiler's enrichKEGG()/gseKEGG() fetch this live on every
#      call; clusterProfiler:::download_KEGG() is the same fetch, done once
#      here and cached, so the compute-node side can use enricher()/GSEA()
#      with a manually-supplied TERM2GENE/TERM2NAME table instead — same
#      statistics, no live call. Fetched once and reused: if kegg_hsa.rds
#      already exists in cache_dir it is NOT re-downloaded, so all species
#      share one identical KEGG snapshot.
#   GO enrichment (org.Hs.eg.db) is NOT fetched here — that is a local
#   Bioconductor SQLite package already installed in the conda env, no network
#   call happens at enrichment time regardless of node.
#
# Usage:
#   fetch_annotation_cache.R <species> <cache_dir>
#
#   species       - "cattle", "sheep" or "goat" (matches run_enrichment.R's
#                   argument convention)
#   cache_dir     - where the cache files get written; run this once per
#                   species. Filenames are species-qualified so species don't
#                   collide, BUT a cattle_annotation_map.rds from the legacy
#                   hybrid design would be overwritten by the new ortholog-only
#                   one — the script refuses to do that (see guard below), so
#                   use a NEW cache_dir if a pilot cache already exists.
#
# Run once per species, ahead of the actual pipeline run:
#   fetch_annotation_cache.R cattle annotation_cache_main/
#   fetch_annotation_cache.R sheep  annotation_cache_main/
#   fetch_annotation_cache.R goat   annotation_cache_main/
#
# Output (per species):
#   <cache_dir>/<species>_annotation_map.rds
#     data.frame: gene_id, entrez_go, source_go, entrez_kegg, source_kegg
#     — one row per gene in the genome; genes with no human ortholog carry NA.
#   <cache_dir>/kegg_hsa.rds
#     list(TERM2GENE = data.frame(pathway_id, entrez_id),
#          TERM2NAME = data.frame(pathway_id, pathway_name))
#     — ready to hand straight to clusterProfiler::enricher()/GSEA() in
#     place of enrichKEGG()/gseKEGG().
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
# --- legacy 3-arg form (gene universe read from a SECTION_7 output dir) ---
#if (length(args) != 3) {
#  stop("Usage: fetch_annotation_cache.R <species> <section7_dir> <cache_dir>")
#}
#species       <- tolower(args[1])
#section7_dir  <- args[2]
#cache_dir     <- args[3]
if (!length(args) %in% c(2, 3)) {
  stop("Usage: fetch_annotation_cache.R <species> <cache_dir> [gene_ids_file]")
}
species   <- tolower(args[1])
cache_dir <- args[2]
args[3]   <- if (length(args) == 3) args[3] else NA_character_  # optional gene_ids_file

if (!species %in% c("cattle", "sheep", "goat")) {
  stop("species must be 'cattle', 'sheep' or 'goat', got: ", species)
}

# ---- BioMart batch resilience settings -------------------------------------
# Ported from fac.R (Suvi's hardened sheep-specific version). NOTE: fac.R's
# own comment says "keep this at 50" but its actual BATCH_SIZE constant was
# 100 — kept at 100 here (the value actually in effect there), flagging the
# mismatch rather than silently picking one.
BATCH_SIZE       <- 100
MAX_ATTEMPTS     <- 5
RETRY_SLEEP_BASE <- 5
MIN_SPLIT_SIZE   <- 5   # below this, give up and record as failed rather than keep halving

suppressPackageStartupMessages({
  library(curl)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  #if (species == "cattle") library(org.Bt.eg.db)   # legacy hybrid path only
  library(clusterProfiler)
  library(dplyr)
})

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Guard: never silently overwrite a legacy hybrid (direct + ortholog) annotation
# map. The pilot's cattle cache is a validated artifact; the new ortholog-only
# map shares its filename, so refuse and make the caller pick a fresh cache_dir.
annot_out <- file.path(cache_dir, paste0(species, "_annotation_map.rds"))
if (file.exists(annot_out)) {
  existing_map <- readRDS(annot_out)
  if (any(existing_map$source_go == "direct", na.rm = TRUE)) {
    stop("Refusing to overwrite a legacy hybrid annotation map (contains 'direct' rows): ",
         annot_out, "\n  Use a new cache_dir for the ortholog-only cache.")
  }
}

# Pinned Ensembl archive host — must match run_enrichment.R exactly.
ARCHIVE_HOST <- "https://jun2026.archive.ensembl.org"
#ENSEMBL_MART <- "ENSEMBL_MART_ENSEMBL"
ARCHIVE_BIOMART_URL <- paste0(ARCHIVE_HOST, "/biomart/martservice")

biomart_dataset <- switch(species,
  cattle = "btaurus_gene_ensembl",
  sheep  = "oaries_gene_ensembl",
  goat   = "chircus_gene_ensembl"
)
# All three species go through human orthologs, so all three use the human
# ("hsa") KEGG pathway table. (Legacy: cattle = "bta" in the hybrid design.)
kegg_organism <- "hsa"

# ---- BioMart connection with retry + logged live-Ensembl fallback ---------
# Identical to run_enrichment.R's version — duplicated deliberately rather
# than shared/sourced, since it's two small self-contained helpers and this
# script needs to keep working standalone regardless of what else changes.
#connect_biomart <- function(dataset) {
#  tryCatch({
#    useEnsembl(biomart = ENSEMBL_MART, dataset = dataset, host = ARCHIVE_HOST)
#  }, error = function(e) {
#    for (m in c("www", "us#east", "asia")) {
#     message("  *** PINNED ARCHIVE HOST UNREACHABLE — falling back to LIVE Ensembl ('", m, "' mirror). ",
#              "This run's annotation will NOT match the pinned release 116 archive — flag this in any output. ***")
#      mart <- tryCatch(
#        useEnsembl(biomart = ENSEMBL_MART, dataset = dataset, mirror = m),
#        error = function(e2) NULL
#      )
#      if (!is.null(mart)) return(mart)
#    }
#    stop("Both the pinned archive host and every live Ensembl mirror failed for dataset: ", dataset)
#  })
#}

#call_with_retry <- function(fn, max_attempts = 5) {
#  last_error <- NULL
#  for (attempt in seq_len(max_attempts)) {
#    result <- tryCatch(fn(), error = function(e) e)
#    if (!inherits(result, "error")) return(result)
#    last_error <- result
#    message("  attempt ", attempt, " failed (", conditionMessage(result), "), retrying...")
#    Sys.sleep(5 * attempt)
#  }
#  stop("Failed after ", max_attempts, " attempts — last error: ", conditionMessage(last_error))
#}

#get_human_ortholog_entrez <- function(gene_ids, source_dataset) {
#  mart_source <- connect_biomart(source_dataset)
#  mart_human  <- connect_biomart("hsapiens_gene_ensembl")
#
#  homology <- call_with_retry(function() {
#    getLDS(
#      attributes  = "ensembl_gene_id", filters = "ensembl_gene_id",
#      values      = gene_ids, mart = mart_source,
#      attributesL = "ensembl_gene_id", martL = mart_human
#    )
#  })
#  colnames(homology) <- c("gene_id", "human_ensembl_gene_id")
#  homology <- homology[!duplicated(homology$gene_id), , drop = FALSE]

#  human_entrez <- AnnotationDbi::select(
#    org.Hs.eg.db, keys = unique(homology$human_ensembl_gene_id),
#    keytype = "ENSEMBL", columns = "ENTREZID"
#  )
#  human_entrez <- human_entrez[!duplicated(human_entrez$ENSEMBL) & !is.na(human_entrez$ENTREZID), ]

#  out <- merge(homology, human_entrez, by.x = "human_ensembl_gene_id", by.y = "ENSEMBL",
#               all.x = TRUE)
#  setNames(out$ENTREZID, out$gene_id)
#}

call_with_retry <- function(fn, max_attempts = MAX_ATTEMPTS, sleep_base = RETRY_SLEEP_BASE) {
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(fn(), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    last_error <- result
    message("  attempt ", attempt, " failed (", conditionMessage(result), "), retrying...")
    if (attempt < max_attempts) Sys.sleep(sleep_base * attempt)
  }
  stop("Failed after ", max_attempts, " attempts — last error: ", conditionMessage(last_error))
}

# One raw request for one batch of IDs. Ported from fac.R's
# fetch_biomart_batch(), generalized to any dataset/filter/attribute list
# (fac.R hardcoded the two homology columns; biomart_query() below reuses
# this for both the homology fetch AND the whole-genome gene-list fetch).
# Two failure modes fac.R found in practice, both checked explicitly:
#   - Ensembl sometimes returns an HTML page with HTTP 200 (not an error
#     status, so this needs its own check, not just response$status_code).
#   - a malformed/unexpected response won't have the expected column count.
fetch_biomart_batch <- function(dataset, filter_name, ids, attributes) {
  ids_string <- paste(ids, collapse = ",")
  attribute_xml <- paste0('<Attribute name="', attributes, '"/>', collapse = "")
  filter_xml <- if (!is.null(filter_name)) {
    paste0('<Filter name="', filter_name, '" value="', ids_string, '"/>')
  } else ""
  query <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" ',
    'formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">',
    '<Dataset name="', dataset, '" interface="default">', filter_xml, attribute_xml,
    '</Dataset></Query>'
  )
  url <- paste0(ARCHIVE_BIOMART_URL, "?query=", curl::curl_escape(query))
  response <- curl::curl_fetch_memory(url)
  if (response$status_code != 200) stop("BioMart HTTP status ", response$status_code)

  txt <- rawToChar(response$content)
  if (!nzchar(trimws(txt))) stop("BioMart returned an empty response.")

  if (grepl("^\\s*<html|^\\s*<!DOCTYPE html|^\\s*<HTML", txt, ignore.case = TRUE)) {
    stop("BioMart returned HTML instead of TSV. First response line: ",
         substr(trimws(strsplit(txt, "\n", fixed = TRUE)[[1]][1]), 1, 200))
  }

  result <- read.delim(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
  expected_cols <- length(attributes)
  if (ncol(result) != expected_cols) {
    stop("Unexpected BioMart column count: got ", ncol(result), ", expected ", expected_cols,
         ". First response line: ", substr(trimws(strsplit(txt, "\n", fixed = TRUE)[[1]][1]), 1, 200))
  }
  result
}

# Batched, self-healing query: a failed batch does NOT terminate the run.
# It's retried; if it still fails, it's split in half and each half tried
# independently (which retries and can split again); below MIN_SPLIT_SIZE,
# those IDs are recorded as failed and the run continues. Ported from
# fac.R's biomart_query()/process_ids(), generalized beyond the homology-
# specific 2-column case, and with filter_name/attributes made explicit
# arguments so this same function covers both the homology fetch and (with
# filter_name = NULL, an unfiltered whole-attribute dump) the whole-genome
# gene-list fetch.
biomart_query <- function(dataset, filter_name, filter_values, attributes,
                           batch_size = BATCH_SIZE, failed_ids_file = NULL) {
  empty_result <- function() {
    as.data.frame(setNames(replicate(length(attributes), character(0), simplify = FALSE), attributes))
  }

  if (is.null(filter_name)) {
    # Unfiltered dump (whole-genome gene list) — nothing to batch/split on,
    # just retry the single request.
    return(list(results = call_with_retry(function() {
      fetch_biomart_batch(dataset, NULL, NULL, attributes)
    }), failed_ids = character(0)))
  }

  filter_values <- unique(filter_values[!is.na(filter_values) & nzchar(filter_values)])
  if (!length(filter_values)) return(list(results = empty_result(), failed_ids = character(0)))

  batches <- split(filter_values, ceiling(seq_along(filter_values) / batch_size))
  results <- list()
  failed_ids <- character()

  process_ids <- function(ids, label = "") {
    if (!length(ids)) return(empty_result())

    result <- tryCatch(
      call_with_retry(function() fetch_biomart_batch(dataset, filter_name, ids, attributes)),
      error = function(e) e
    )
    if (!inherits(result, "error")) return(result)

    message("  BioMart request failed", if (nzchar(label)) paste0(" (", label, ")"),
            ": ", conditionMessage(result))

    if (length(ids) <= MIN_SPLIT_SIZE) {
      message("  Giving up on these ", length(ids), " IDs for now; recording as failed.")
      failed_ids <<- c(failed_ids, ids)
      return(empty_result())
    }

    mid <- floor(length(ids) / 2)
    left_ids <- ids[seq_len(mid)]
    right_ids <- ids[(mid + 1):length(ids)]
    message("  Splitting ", length(ids), " IDs -> ", length(left_ids), " + ", length(right_ids))

    dplyr::bind_rows(
      process_ids(left_ids,  label = paste0("left ",  length(left_ids),  " IDs")),
      process_ids(right_ids, label = paste0("right ", length(right_ids), " IDs"))
    )
  }

  for (i in seq_along(batches)) {
    message("BioMart batch ", i, "/", length(batches), " (", length(batches[[i]]), " IDs)")
    batch_result <- process_ids(batches[[i]], label = paste0("batch ", i))
    if (nrow(batch_result)) results[[length(results) + 1L]] <- batch_result

    # Written after every batch, not just at the end, so an interrupted job
    # still leaves a useful partial failure list.
    if (!is.null(failed_ids_file) && length(failed_ids)) {
      writeLines(unique(failed_ids), failed_ids_file)
    }
  }

  combined <- if (length(results)) dplyr::bind_rows(results) else empty_result()
  list(results = combined, failed_ids = unique(failed_ids))
}

get_human_ortholog_entrez <- function(gene_ids, source_dataset) {
  gene_ids <- unique(na.omit(as.character(gene_ids)))
  if (!length(gene_ids)) return(list(entrez = setNames(character(0), character(0)), failed_ids = character(0)))

  failed_ids_file <- file.path(cache_dir, paste0("failed_biomart_ids_", species, ".txt"))

  bm <- biomart_query(
    dataset=source_dataset,
    filter_name="ensembl_gene_id",
    filter_values=gene_ids,
    attributes=c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
    batch_size=BATCH_SIZE,
    failed_ids_file=failed_ids_file
  )
  homology <- bm$results
  failed_ids <- bm$failed_ids

  cat("BioMart homology fetch: ", nrow(homology), " rows returned, ",
      length(failed_ids), " gene(s) failed after retries/splitting.\n", sep = "")
  if (length(failed_ids)) {
    cat("Failed gene IDs written to: ", failed_ids_file, "\n", sep = "")
  }

  if (!nrow(homology)) {
    warning("BioMart returned no human orthologs for ", length(gene_ids), " genes.")
    return(list(entrez = setNames(rep(NA_character_, length(gene_ids)), gene_ids), failed_ids = failed_ids))
  }

  colnames(homology)[1:2] <- c("gene_id", "human_ensembl_gene_id")
  homology <- homology[!is.na(homology$human_ensembl_gene_id) & homology$human_ensembl_gene_id != "" & !duplicated(homology$gene_id), , drop=FALSE]

  human_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys=unique(homology$human_ensembl_gene_id),
    keytype="ENSEMBL",
    columns="ENTREZID"
  )
  human_entrez <- human_entrez[!is.na(human_entrez$ENTREZID) & !duplicated(human_entrez$ENSEMBL), , drop=FALSE]

  out <- merge(homology, human_entrez, by.x="human_ensembl_gene_id", by.y="ENSEMBL", all.x=TRUE)
  out <- out[!duplicated(out$gene_id), , drop=FALSE]

  result <- setNames(out$ENTREZID, out$gene_id)
  missing <- setdiff(gene_ids, names(result))
  if (length(missing)) result <- c(result, setNames(rep(NA_character_, length(missing)), missing))
  list(entrez = result[gene_ids], failed_ids = failed_ids)
}

# ---- 1. Gene universe -------------------------------------------------------
# --- legacy: gene universe read back from an existing DESeq2 result dir ---
# Same source run_enrichment.R itself uses: one _DESeq2_exploratory.tsv file
# is enough, since every file for a species shares the same gene_id universe.
#expl_files <- list.files(section7_dir, pattern = "_DESeq2_exploratory\\.tsv$",
#                         recursive = TRUE, full.names = TRUE)
#if (!length(expl_files)) stop("No _DESeq2_exploratory.tsv files found under ", section7_dir)
#gene_lists <- lapply(expl_files, function(f) {
#  x <- read.delim(f, stringsAsFactors = FALSE)
#  if (!"gene_id" %in% names(x)) stop("Missing gene_id column in: ", f)
#  x$gene_id
#})
#all_genes <- unique(unlist(gene_lists, use.names = FALSE))
#all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
#cat("Found", length(expl_files), "DESeq2 result files containing",
#    length(all_genes), "unique genes across all DPIs/comparisons.\n")

# NEW: whole-genome gene universe. TWO sources supported:
#   - optional 3rd CLI arg <gene_ids_file>: a plain-text list of gene IDs,
#     one per line (Suvi's own workflow: extracted from the reference GTF
#     already on disk, so this step needs no network call at all).
#   - otherwise: fetched directly from BioMart for this species' dataset
#     (needs network, but means the script is fully self-contained with no
#     other input required).
# Either way, this no longer depends on a completed DESeq2 run at all — it
# can (and should) be run once, ahead of any analysis, and the resulting
# cache covers every gene in the genome rather than only genes that
# happened to appear in a specific DE result.

if (!is.na(args[3]) && nzchar(args[3])) {
  gene_ids_file <- args[3]
  cat("Reading gene universe from file:", gene_ids_file, "\n")
  all_genes <- readLines(gene_ids_file)
  all_genes <- trimws(all_genes)
  all_genes <- unique(all_genes[!is.na(all_genes) & nzchar(all_genes)])
} else {
  cat("Fetching whole-genome gene list for", biomart_dataset, "from BioMart...\n")
  genome_fetch <- biomart_query(
    dataset = biomart_dataset, filter_name = NULL, filter_values = NULL,
    attributes = "ensembl_gene_id"
  )
  all_genes <- unique(as.character(genome_fetch$results[[1]]))
  all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
}
if (!length(all_genes)) stop("No genes found — check the gene_ids_file or the BioMart dataset: ", biomart_dataset)

cat("Gene universe:", length(all_genes), "genes for", species, ".\n")

# ---- 2. Annotation map (BioMart homology + Entrez) -------------------------
# ALL THREE SPECIES now go through the SAME unconditional human-ortholog path
# — the cattle hybrid (direct org.Bt.eg.db + ortholog fallback for gaps only)
# is retired per decision, kept below as a commented-out block, not deleted.

# --- legacy: cattle hybrid (org.Bt.eg.db direct + ortholog fallback) ---
#if (species == "cattle") {
#
#  direct_map <- AnnotationDbi::select(
#    org.Bt.eg.db, keys = all_genes, keytype = "ENSEMBL", columns = "ENTREZID"
#  )
#  direct_map <- direct_map[!duplicated(direct_map$ENSEMBL), ]
#
#  mapped_entrez <- na.omit(unique(direct_map$ENTREZID))
#  go_check <- AnnotationDbi::select(
#    org.Bt.eg.db, keys = mapped_entrez, keytype = "ENTREZID", columns = "GOALL"
#  )
#  entrez_with_go <- unique(go_check$ENTREZID[!is.na(go_check$GOALL)])
#  direct_map$has_direct_go <- !is.na(direct_map$ENTREZID) & direct_map$ENTREZID %in% entrez_with_go
#
#  gap_genes <- direct_map$ENSEMBL[!direct_map$has_direct_go]
#  cat(length(gap_genes), "of", length(all_genes),
#      "genes lack direct GO annotation in org.Bt.eg.db; fetching human ortholog for these...\n")
#
#  ortholog_entrez <- if (length(gap_genes) > 0) {
#    get_human_ortholog_entrez(gap_genes, biomart_dataset)
#  } else {
#    setNames(character(0), character(0))
#  }
#
#  annot_map <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE) %>%
#    left_join(direct_map, by = c("gene_id" = "ENSEMBL")) %>%
#    mutate(
#      entrez_go   = ifelse(has_direct_go, ENTREZID, ortholog_entrez[gene_id]),
#      source_go   = ifelse(has_direct_go, "direct",
#                            ifelse(!is.na(ortholog_entrez[gene_id]), "human_ortholog", NA)),
#      entrez_kegg = ENTREZID,
#      source_kegg = ifelse(!is.na(ENTREZID), "direct", NA)
#    ) %>%
#    select(gene_id, entrez_go, source_go, entrez_kegg, source_kegg)
#
#} else {
#
#  cat("Fetching human ortholog for all", length(all_genes), "sheep genes (unconditional)...\n")
#  ortholog_entrez <- get_human_ortholog_entrez(all_genes, biomart_dataset)
#
#  annot_map <- data.frame(
#    gene_id   = all_genes,
#    entrez_go = ortholog_entrez[all_genes],
#    source_go = ifelse(!is.na(ortholog_entrez[all_genes]), "human_ortholog", NA),
#    stringsAsFactors = FALSE
#  ) %>%
#    mutate(entrez_kegg = entrez_go, source_kegg = source_go)
#}

# NEW: unconditional human ortholog for cattle, sheep, AND goat — no
# species branch at all. This is the exact same code sheep already used;
# it now applies uniformly.
cat("Fetching human ortholog for all", length(all_genes), species, "genes (unconditional)...\n")
ortholog_result <- get_human_ortholog_entrez(all_genes, biomart_dataset)
ortholog_entrez <- ortholog_result$entrez

annot_map <- data.frame(
  gene_id   = all_genes,
  entrez_go = ortholog_entrez[all_genes],
  source_go = ifelse(!is.na(ortholog_entrez[all_genes]), "human_ortholog", NA),
  stringsAsFactors = FALSE
) %>%
  mutate(entrez_kegg = entrez_go, source_kegg = source_go)

if (length(ortholog_result$failed_ids)) {
  cat(length(ortholog_result$failed_ids),
      "gene(s) never resolved after retries/splitting — see failed_biomart_ids_",
      species, ".txt in ", cache_dir, ". They're kept in annot_map with NA, not silently dropped.\n", sep = "")
}

n_annotated <- sum(!is.na(annot_map$entrez_go))
cat("Annotation map:", n_annotated, "of", length(all_genes), "genes have a usable GO entrez ID.\n")

# annot_out already set near the top of the script (used there for the
# legacy-cache overwrite guard) — reused here rather than redefined.
saveRDS(annot_map, annot_out)
cat("Wrote", annot_out, "\n")

# ---- 3. KEGG term<->gene / term<->name tables -------------------------------
# download_KEGG() is the same internal fetch enrichKEGG()/gseKEGG() do on
# every live call — done once here and cached instead. Its return list's
# KEGGPATHID2EXTID / KEGGPATHID2NAME data frames are already exactly the
# TERM2GENE / TERM2NAME shape enricher()/GSEA() expect.
#
# All three species now use organism "hsa" (unconditional human ortholog),
# so this table is IDENTICAL across cattle/sheep/goat — fetched once and
# reused rather than re-downloaded per species. Whichever species runs this
# script first populates it; the other two just find it already there.

kegg_out <- file.path(cache_dir, paste0("kegg_", kegg_organism, ".rds"))

if (file.exists(kegg_out)) {

  cat("KEGG cache for organism '", kegg_organism, "' already exists at ", kegg_out,
      " — reusing it, not re-fetching.\n", sep = "")

} else {

  cat("Fetching KEGG pathway data for organism '", kegg_organism, "'...\n", sep = "")

  kegg_data <- call_with_retry(function() {
    clusterProfiler:::download_KEGG(kegg_organism, keggType = "KEGG", keyType = "kegg")
  })

  term2gene <- kegg_data$KEGGPATHID2EXTID
  term2name <- kegg_data$KEGGPATHID2NAME
  colnames(term2gene) <- c("pathway_id", "entrez_id")
  colnames(term2name) <- c("pathway_id", "pathway_name")

  saveRDS(list(TERM2GENE = term2gene, TERM2NAME = term2name), kegg_out)
  cat("Wrote", kegg_out, "(", nrow(term2gene), "pathway-gene pairs,",
      length(unique(term2gene$pathway_id)), "pathways )\n")
}

cat("Done. Cache ready in:", cache_dir, "\n")
