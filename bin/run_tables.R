#!/usr/bin/env Rscript
# ===================================================================================
# run_tables.R : flat, browsable tables TSVs
#
# Reads run_pca_deseq2.R's SECTION_7 TSVs and Phase 4's enrichResult/gseaResult .rds
# objects; nothing here recomputes anything, it only flattens/combines what
# already exists into single tables someone can open in Excel, filter, and
# sort instead of hunting through ~150+ small per-contrast files.
#
# Every gene-level table includes BOTH Ensembl and Entrez IDs, plus gene
# symbol and annotation_source, via Phase 4's cached annotation map — same
# cache run_enrichment.R and run_plots.R already use, nothing new fetched.
#
# Six tables:
#   1. master_candidate_genes  — every candidate gene, every comparison/dpi/
#      direction, one flat table
#   2. comparison_summary      — one row per (comparison, dpi): counts,
#      totals, QC — the "landing page" table
#   3. master_enrichment_ora   — every GO/KEGG ORA hit, every comparison/
#      dpi/direction, flattened
#   4. temporal_trajectories   — one row per candidate gene per comparison,
#      log2FC across dpi (wide format)
#   5. gsea_summary            — every GO/KEGG GSEA result, flattened
#      (kept separate from #3 — GSEA scores the whole ranked list, ORA
#      tests a candidate subset; different question, shouldn't be merged)
#   6. gene_pathway_crosswalk  — one row per (gene, GO/KEGG term) pair,
#      built by exploding table 3's gene lists rather than re-reading any
#      RDS — "what pathways is gene X in". ORA-scoped only (every gene in
#      here was already a candidate gene by construction, since ORA's
#      input gene sets ARE the candidate lists); GSEA's core_enrichment
#      genes are deliberately not mixed in here, same separation-of-
#      concerns reasoning as keeping table 5 apart from table 3.
#
# Usage:
#   run_tables.R <section7_outdir> <enrichment_dir> <species> <outdir> <cache_dir>
# ===================================================================================

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5) stop("Usage: run_tables.R <section7_dir> <enrichment_dir> <species> <outdir> <cache_dir>")

section7_dir <- args[1]
enrichment_dir <- args[2]
species <- tolower(args[3])
outdir <- args[4]
cache_dir <- args[5]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(AnnotationDbi)
  library(DOSE)
  library(org.Hs.eg.db)
  if(species=="cattle") library(org.Bt.eg.db)
})

tables_dir <- outdir
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Annotation map (Entrez ID, symbol, source) ----------------------------

annot_map_file <- file.path(cache_dir,paste0(species,"_annotation_map.rds"))
if(!file.exists(annot_map_file)) stop("Annotation cache not found: ",annot_map_file)
annot_map <- readRDS(annot_map_file)

direct_ids <- unique(na.omit(annot_map$entrez_go[annot_map$source_go=="direct"]))
ortholog_ids <- unique(na.omit(annot_map$entrez_go[annot_map$source_go=="human_ortholog"]))

sym_direct <- if(species=="cattle" && length(direct_ids)>0) {
  AnnotationDbi::select(org.Bt.eg.db,keys=direct_ids,keytype="ENTREZID",columns="SYMBOL")
} else {
  data.frame(ENTREZID=character(0),SYMBOL=character(0))
}

sym_ortholog <- if(length(ortholog_ids)>0) {
  AnnotationDbi::select(org.Hs.eg.db,keys=ortholog_ids,keytype="ENTREZID",columns="SYMBOL")
} else {
  data.frame(ENTREZID=character(0),SYMBOL=character(0))
}

sym_lookup <- rbind(sym_direct,sym_ortholog)
sym_lookup <- sym_lookup[!duplicated(sym_lookup$ENTREZID) & !is.na(sym_lookup$SYMBOL),]
symbol_for_entrez <- setNames(sym_lookup$SYMBOL,sym_lookup$ENTREZID)

gene_lookup <- annot_map %>%
  dplyr::transmute(
    gene_id,
    entrez_id=entrez_go,
    gene_symbol=unname(symbol_for_entrez[as.character(entrez_go)]),
    annotation_source=source_go
  )

entrez_to_ensembl <- gene_lookup %>%
  dplyr::filter(!is.na(entrez_id)) %>%
  dplyr::group_by(entrez_id) %>%
  dplyr::summarise(gene_id=paste(unique(gene_id),collapse=","),.groups="drop")

entrez_to_ensembl_map <- setNames(entrez_to_ensembl$gene_id,entrez_to_ensembl$entrez_id)

map_entrez_list <- function(entrez_str) {
  vapply(entrez_str, function(x) {
    if (is.na(x) || !nzchar(x)) return("")

    ids <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
    mapped <- unname(entrez_to_ensembl_map[ids])
    mapped <- mapped[!is.na(mapped)]

    paste(mapped, collapse = "/")
  }, character(1))
}
cat("Annotation map loaded:",sum(!is.na(gene_lookup$entrez_id)),"of",nrow(gene_lookup),"genes have an Entrez ID.\n")

# ============================================================================
# 1. Master candidate gene table
# ============================================================================

candidate_files <- list.files(section7_dir,pattern="_candidate_(up|down)\\.tsv$",recursive=TRUE,full.names=TRUE)

parse_candidate_path <- function(path) {
  comparison_dir <- basename(dirname(path))
  fname <- basename(path)
  m <- regmatches(fname,regexec("_dpi([^_]+)_candidate_(up|down)\\.tsv$",fname))[[1]]
  if(length(m)!=3) stop("Could not parse candidate filename: ",fname)
  list(comparison=comparison_dir,dpi=m[2],direction=m[3])
}

candidate_rows <- lapply(candidate_files,function(f) {
  info <- parse_candidate_path(f)
  df <- read.delim(f,stringsAsFactors=FALSE)
  df$comparison <- info$comparison
  df$dpi <- info$dpi
  df$direction <- info$direction
  df
})

master_candidates <- dplyr::bind_rows(candidate_rows) %>%
  dplyr::left_join(gene_lookup,by="gene_id") %>%
  dplyr::select(comparison,dpi,direction,gene_id,gene_symbol,entrez_id,annotation_source,dplyr::everything())

write.table(master_candidates,file.path(tables_dir,paste0("1_master_candidate_genes_",species,".tsv")),
            sep="\t",quote=FALSE,row.names=FALSE)
cat("Wrote master_candidate_genes (",nrow(master_candidates),"rows).\n")

# ============================================================================
# 2. Comparison summary table
# ============================================================================

expl_files <- list.files(section7_dir,pattern="_DESeq2_exploratory\\.tsv$",recursive=TRUE,full.names=TRUE)

parse_expl_path <- function(path) {
  comparison_dir <- basename(dirname(path))
  fname <- basename(path)
  m <- regmatches(fname,regexec("_dpi([^_]+)_DESeq2_exploratory\\.tsv$",fname))[[1]]
  if(length(m)!=2) stop("Could not parse exploratory filename: ",fname)
  list(comparison=comparison_dir,dpi=m[2])
}

summary_rows <- lapply(expl_files,function(f) {
  info <- parse_expl_path(f)
  df <- read.delim(f,stringsAsFactors=FALSE)
  data.frame(
    comparison=info$comparison,
    dpi=info$dpi,
    genes_tested=sum(!is.na(df$padj)),
    candidates_up=sum(df$candidate & df$log2FoldChange>0,na.rm=TRUE),
    candidates_down=sum(df$candidate & df$log2FoldChange<0,na.rm=TRUE),
    pct_technical_replicates_consistent=round(100*mean(df$technical_replicates_consistent,na.rm=TRUE),1),
    stringsAsFactors=FALSE
  )
})

comparison_summary <- dplyr::bind_rows(summary_rows) %>%
  dplyr::arrange(comparison,as.numeric(dpi))

write.table(comparison_summary,file.path(tables_dir,paste0("2_comparison_summary_",species,".tsv")),
              sep="\t",quote=FALSE,row.names=FALSE)
cat("Wrote comparison_summary (",nrow(comparison_summary),"rows).\n")

# ============================================================================
# 3. Master enrichment results table (ORA — GO + KEGG)
# ============================================================================

go_ora_files <- list.files(enrichment_dir,pattern="_GO_results\\.rds$",recursive=TRUE,full.names=TRUE)
go_ora_files <- go_ora_files[!grepl("_GSEA_",go_ora_files)]

kegg_ora_files <- list.files(enrichment_dir,pattern="_KEGG_result\\.rds$",recursive=TRUE,full.names=TRUE)
kegg_ora_files <- kegg_ora_files[!grepl("_GSEA_",kegg_ora_files)]

# ORA structure:
# enrichment/<comparison>/dpi<N>/ORA/<direction>/<file>.rds

parse_enrichment_path <- function(path) {
  direction <- basename(dirname(path))
  ora_folder <- basename(dirname(dirname(path)))
  dpi_folder <- basename(dirname(dirname(dirname(path))))
  comparison <- basename(dirname(dirname(dirname(dirname(path)))))
  if(ora_folder!="ORA") stop("Expected ORA directory in enrichment path: ",path)
  list(comparison=comparison,dpi=sub("^dpi","",dpi_folder),direction=direction)
}

go_ora_rows <- lapply(go_ora_files,function(f) {
  info <- parse_enrichment_path(f)
  go_results <- readRDS(f)

  dplyr::bind_rows(lapply(names(go_results),function(src) {
    df <- as.data.frame(go_results[[src]])
    if(nrow(df)==0) return(NULL)

    df$comparison <- info$comparison
    df$dpi <- info$dpi
    df$direction <- info$direction
    df$annotation_source <- src
    df$term_type <- if("ONTOLOGY" %in% names(df)) df$ONTOLOGY else "GO"
    df$gene_id <- map_entrez_list(df$geneID)
    df
  }))
})

kegg_ora_rows <- lapply(kegg_ora_files,function(f) {
  info <- parse_enrichment_path(f)
  df <- as.data.frame(readRDS(f))
  if(nrow(df)==0) return(NULL)

  df$comparison <- info$comparison
  df$dpi <- info$dpi
  df$direction <- info$direction
  df$annotation_source <- if(species=="cattle") "direct" else "human_ortholog"
  df$term_type <- "KEGG"
  df$gene_id <- map_entrez_list(df$geneID)
  df
})

master_ora <- dplyr::bind_rows(c(go_ora_rows,kegg_ora_rows))

if(nrow(master_ora)>0) {
  master_ora <- master_ora %>%
    dplyr::select(
      comparison,dpi,direction,term_type,annotation_source,
      ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue,
      Count,geneID,gene_id,dplyr::everything()
    )
}

write.table(master_ora,file.path(tables_dir,paste0("3_master_enrichment_ora_",species,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
cat("Wrote master_enrichment_ora (",nrow(master_ora),"rows).\n")

# ============================================================================
# 4. Gene-to-pathway crosswalk
# ============================================================================

if(nrow(master_ora)>0 && "gene_id" %in% names(master_ora)) {

  gene_pathway_crosswalk <- master_ora %>%
    dplyr::filter(!is.na(gene_id) & gene_id!="") %>%
    tidyr::separate_rows(gene_id,sep="/") %>%
    tidyr::separate_rows(gene_id,sep=",") %>%
    dplyr::filter(gene_id!="") %>%
    dplyr::left_join(gene_lookup %>% dplyr::select(gene_id,gene_symbol,entrez_id),by="gene_id") %>%
    dplyr::select(gene_id,gene_symbol,entrez_id,annotation_source,comparison,dpi,direction,term_type,ID,Description,p.adjust) %>%
    dplyr::distinct() %>%
    dplyr::arrange(gene_id,comparison,dpi,direction)

} else {

  gene_pathway_crosswalk <- data.frame(
    gene_id=character(0),
    gene_symbol=character(0),
    entrez_id=character(0),
    annotation_source=character(0),
    comparison=character(0),
    dpi=character(0),
    direction=character(0),
    term_type=character(0),
    ID=character(0),
    Description=character(0),
    p.adjust=numeric(0),
    stringsAsFactors=FALSE
  )

  cat("No ORA enrichment results found; writing an empty gene_pathway_crosswalk table.\n")
}

write.table(gene_pathway_crosswalk,file.path(tables_dir,paste0("4_gene_pathway_crosswalk_",species,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
cat("Wrote gene_pathway_crosswalk (",nrow(gene_pathway_crosswalk),"rows).\n")

# ============================================================================
# 5. Temporal trajectory table
# ============================================================================

comparison_dirs <- list.dirs(section7_dir,recursive=FALSE,full.names=FALSE)

trajectory_rows <- lapply(comparison_dirs,function(cd) {
  temporal_file <- list.files(file.path(section7_dir,cd),pattern="_temporal_summary\\.tsv$",full.names=TRUE)
  lfc_file <- list.files(file.path(section7_dir,cd),pattern="_infection_LFC_matrix\\.tsv$",full.names=TRUE)

  if(length(temporal_file)==0 || length(lfc_file)==0) return(NULL)

  temporal <- read.delim(temporal_file[1],stringsAsFactors=FALSE)
  lfc_mat <- read.delim(lfc_file[1],stringsAsFactors=FALSE)

  candidates <- temporal %>% dplyr::filter(passed_candidate_filter_any_DPI)

  lfc_mat %>%
    dplyr::filter(gene_id %in% candidates$gene_id) %>%
    dplyr::left_join(candidates %>% dplyr::select(gene_id,max_abs_log2FC,DPI_of_max_abs_log2FC),by="gene_id") %>%
    dplyr::mutate(comparison=cd)
})

temporal_trajectories <- dplyr::bind_rows(trajectory_rows) %>%
  dplyr::left_join(gene_lookup,by="gene_id") %>%
  dplyr::select(comparison,gene_id,gene_symbol,entrez_id,annotation_source,max_abs_log2FC,DPI_of_max_abs_log2FC,dplyr::everything())

write.table(temporal_trajectories,file.path(tables_dir,paste0("5_temporal_trajectories_",species,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
cat("Wrote temporal_trajectories (",nrow(temporal_trajectories),"rows).\n")

# ============================================================================
# 6. GSEA summary table — GO + KEGG
# ============================================================================

go_gsea_files <- list.files(enrichment_dir,pattern="_GSEA_GO_results\\.rds$",recursive=TRUE,full.names=TRUE)
kegg_gsea_files <- list.files(enrichment_dir,pattern="_GSEA_KEGG_result\\.rds$",recursive=TRUE,full.names=TRUE)

# GSEA structure:
# enrichment/<comparison>/dpi<N>/GSEA/<file>.rds

parse_gsea_path <- function(path) {
  gsea_folder <- basename(dirname(path))
  dpi_folder <- basename(dirname(dirname(path)))
  comparison <- basename(dirname(dirname(dirname(path))))
  if(gsea_folder!="GSEA") stop("Expected GSEA directory in enrichment path: ",path)
  list(comparison=comparison,dpi=sub("^dpi","",dpi_folder))
}

go_gsea_rows <- lapply(go_gsea_files,function(f) {
  info <- parse_gsea_path(f)
  go_results <- readRDS(f)

  dplyr::bind_rows(lapply(names(go_results),function(src) {
    df <- as.data.frame(go_results[[src]])
    if(nrow(df)==0) return(NULL)

    df$comparison <- info$comparison
    df$dpi <- info$dpi
    df$annotation_source <- src
    df$term_type <- if("ONTOLOGY" %in% names(df)) df$ONTOLOGY else "GO"

    gene_col <- if("core_enrichment" %in% names(df)) "core_enrichment" else "geneID"
    df$gene_id <- map_entrez_list(df[[gene_col]])
    df
  }))
})

kegg_gsea_rows <- lapply(kegg_gsea_files,function(f) {
  info <- parse_gsea_path(f)
  df <- as.data.frame(readRDS(f))
  if(nrow(df)==0) return(NULL)

  df$comparison <- info$comparison
  df$dpi <- info$dpi
  df$annotation_source <- if(species=="cattle") "direct" else "human_ortholog"
  df$term_type <- "KEGG"

  gene_col <- if("core_enrichment" %in% names(df)) "core_enrichment" else "geneID"
  df$gene_id <- map_entrez_list(df[[gene_col]])
  df
})

gsea_summary <- dplyr::bind_rows(c(go_gsea_rows,kegg_gsea_rows))

if(nrow(gsea_summary)>0) {
  if(!"qvalues" %in% names(gsea_summary)) {
    if("qvalue" %in% names(gsea_summary)) gsea_summary$qvalues <- gsea_summary$qvalue
    else gsea_summary$qvalues <- NA_real_
  }

  gsea_summary <- gsea_summary %>%
    dplyr::select(
      comparison,dpi,term_type,annotation_source,
      ID,Description,setSize,enrichmentScore,NES,
      pvalue,p.adjust,qvalues,gene_id,dplyr::everything()
    )
}

write.table(gsea_summary,file.path(tables_dir,paste0("6_gsea_summary_",species,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
cat("Wrote gsea_summary (",nrow(gsea_summary),"rows).\n")

cat("All tables written to: ",tables_dir,"\n")
