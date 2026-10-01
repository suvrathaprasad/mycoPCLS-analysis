#!/usr/bin/env Rscript
# -------------------------------------------------------------------------------
# run_enrichment.R : functional enrichment (GO/KEGG, ORA + GSEA)
#
# plots produced here are individual enrichment result which are standard
# clusterProfiler visualisations. plots here are ONLY for QC. the 
# RDS objects produced here are more important as they contain enrichment results
#
# Usage: run_enrichment.R <deseq2_outdir> <species> <outdir> <cache_dir>
# -------------------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4) stop("Usage: run_enrichment.R <deseq2_outdir> <species> <outdir> <cache_dir>")

deseq2_outdir <- args[1]
species <- tolower(args[2])
outdir <- args[3]
cache_dir <- args[4]

if(!species %in% c("cattle","sheep")) stop("species must be 'cattle' or 'sheep', got: ",species)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  if(species=="cattle") library(org.Bt.eg.db)
  library(dplyr)
  library(ggplot2)
})

section7_dir <- file.path(deseq2_outdir,"SECTION_7")
if(!dir.exists(section7_dir)) stop("SECTION_7 directory not found: ",section7_dir)

enrichment_dir <- file.path(outdir,"enrichment")
dir.create(enrichment_dir,recursive=TRUE,showWarnings=FALSE)

# ---- 1. Annotation caches --------------------------------------------------

annotation_file <- file.path(cache_dir,paste0(species,"_annotation_map.rds"))
if(!file.exists(annotation_file)) stop("Annotation cache not found: ",annotation_file)
annot_map <- readRDS(annotation_file)

kegg_organism <- if(species=="cattle") "bta" else "hsa"
kegg_file <- file.path(cache_dir,paste0("kegg_",kegg_organism,".rds"))
if(!file.exists(kegg_file)) stop("KEGG cache not found: ",kegg_file)
kegg_cache <- readRDS(kegg_file)

TERM2GENE <- kegg_cache$TERM2GENE
TERM2NAME <- kegg_cache$TERM2NAME

TERM2GENE <- data.frame(pathway_id=as.character(TERM2GENE[[1]]),entrez_id=as.character(TERM2GENE[[2]]),stringsAsFactors=FALSE)
TERM2NAME <- data.frame(pathway_id=as.character(TERM2NAME[[1]]),pathway_name=as.character(TERM2NAME[[2]]),stringsAsFactors=FALSE)

cat("Loaded annotation cache:",annotation_file,"\n")
cat("Loaded KEGG cache:",kegg_file,"\n")
cat("Annotation map:",nrow(annot_map),"genes\n")
cat("KEGG:",nrow(TERM2GENE),"pathway-gene pairs\n")

# ---- 2. Discover DESeq2 outputs -------------------------------------------

candidate_files <- list.files(section7_dir,pattern="_candidate_(up|down)\\.tsv$",recursive=TRUE,full.names=TRUE)
if(length(candidate_files)==0) stop("No _candidate_up.tsv/_candidate_down.tsv files found under ",section7_dir)

exploratory_files <- list.files(section7_dir,pattern="_DESeq2_exploratory\\.tsv$",recursive=TRUE,full.names=TRUE)
if(length(exploratory_files)==0) stop("No _DESeq2_exploratory.tsv files found under ",section7_dir)

parse_candidate_path <- function(path) {
  comparison_dir <- basename(dirname(path)); fname <- basename(path)
  m <- regmatches(fname,regexec("_dpi([^_]+)_candidate_(up|down)\\.tsv$",fname))[[1]]
  if(length(m)!=3) stop("Could not parse dpi/direction from filename: ",fname)
  list(comparison_dir=comparison_dir,dpi=m[2],direction=m[3])
}

parse_exploratory_path <- function(path) {
  comparison_dir <- basename(dirname(path)); fname <- basename(path)
  m <- regmatches(fname,regexec("_dpi([^_]+)_DESeq2_exploratory\\.tsv$",fname))[[1]]
  if(length(m)!=2) stop("Could not parse dpi from filename: ",fname)
  list(comparison_dir=comparison_dir,dpi=m[2])
}

MIN_CANDIDATES <- 5
EXPRESSION_FLOOR <- 10

# ---- 3. Annotation helpers -------------------------------------------------

to_entrez_subsets <- function(gene_ids,id_col,source_col) {
  sub <- annot_map[annot_map$gene_id %in% gene_ids,c("gene_id",id_col,source_col)]
  names(sub) <- c("gene_id","entrez","source")
  sub <- sub[!is.na(sub$entrez),]
  list(
    direct=unique(as.character(sub$entrez[sub$source=="direct"])),
    ortholog=unique(as.character(sub$entrez[sub$source=="human_ortholog"]))
  )
}

# ---- 4. GO ORA -------------------------------------------------------------

run_go_ora <- function(candidate_genes,background_genes,tag) {
  cand <- to_entrez_subsets(candidate_genes,"entrez_go","source_go")
  bg <- to_entrez_subsets(background_genes,"entrez_go","source_go")
  results <- list()

  if(length(cand$direct)>0 && length(bg$direct)>0) {
    res <- tryCatch(
      enrichGO(gene=cand$direct,universe=bg$direct,OrgDb=org.Bt.eg.db,keyType="ENTREZID",ont="ALL",pAdjustMethod="BH"),
      error=function(e){message(tag," GO ORA (direct) failed: ",conditionMessage(e));NULL}
    )
    if(!is.null(res) && nrow(as.data.frame(res))>0) results$direct <- res
  }

  if(length(cand$ortholog)>0 && length(bg$ortholog)>0) {
    res <- tryCatch(
      enrichGO(gene=cand$ortholog,universe=bg$ortholog,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="ALL",pAdjustMethod="BH"),
      error=function(e){message(tag," GO ORA (ortholog) failed: ",conditionMessage(e));NULL}
    )
    if(!is.null(res) && nrow(as.data.frame(res))>0) results$ortholog <- res
  }

  results
}

# ---- 5. KEGG ORA -----------------------------------------------------------

run_kegg_ora <- function(candidate_genes,background_genes,tag) {
  cand <- unique(na.omit(as.character(annot_map$entrez_kegg[annot_map$gene_id %in% candidate_genes])))
  bg <- unique(na.omit(as.character(annot_map$entrez_kegg[annot_map$gene_id %in% background_genes])))
  if(length(cand)==0 || length(bg)==0) return(NULL)

  src_label <- if(species=="cattle") "direct" else "human_ortholog"

  res <- tryCatch(
    enricher(gene=cand,universe=bg,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,pAdjustMethod="BH"),
    error=function(e){message(tag," KEGG ORA failed: ",conditionMessage(e));NULL}
  )

  if(is.null(res) || nrow(as.data.frame(res))==0) return(NULL)
  list(result=res,source_label=src_label)
}

# ---- 6. Result writers -----------------------------------------------------

write_go_results <- function(go_results,out_dir,tag) {
  if(length(go_results)==0) return(invisible(NULL))
  tables <- list()

  for(src in names(go_results)) {
    df <- as.data.frame(go_results[[src]])
    if(nrow(df)==0) next
    df$annotation_source <- src
    tables[[src]] <- df

    p <- tryCatch(
      enrichplot::dotplot(go_results[[src]],showCategory=20)+ggtitle(paste0(tag," GO (",src,")")),
      error=function(e) NULL
    )

    if(!is.null(p)) ggsave(file.path(out_dir,paste0(tag,"_GO_",src,"_dotplot.pdf")),p,width=9,height=7)
  }

  if(length(tables)==0) return(invisible(NULL))
  write.table(bind_rows(tables),file.path(out_dir,paste0(tag,"_GO.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  saveRDS(go_results,file.path(out_dir,paste0(tag,"_GO_results.rds")))
}

write_kegg_results <- function(kegg_out,out_dir,tag) {
  if(is.null(kegg_out)) return(invisible(NULL))

  df <- as.data.frame(kegg_out$result)
  if(nrow(df)==0) return(invisible(NULL))
  df$annotation_source <- kegg_out$source_label

  write.table(df,file.path(out_dir,paste0(tag,"_KEGG.tsv")),sep="\t",quote=FALSE,row.names=FALSE)

  p <- tryCatch(
    enrichplot::dotplot(kegg_out$result,showCategory=20)+ggtitle(paste0(tag," KEGG (",kegg_out$source_label,")")),
    error=function(e) NULL
  )

  if(!is.null(p)) ggsave(file.path(out_dir,paste0(tag,"_KEGG_dotplot.pdf")),p,width=9,height=7)
  saveRDS(kegg_out$result,file.path(out_dir,paste0(tag,"_KEGG_result.rds")))
}

# ---- 7. ORA ----------------------------------------------------------------

expl_index <- new.env()

for(ef in exploratory_files) {
  einfo <- parse_exploratory_path(ef)
  assign(paste(einfo$comparison_dir,einfo$dpi,sep="||"),ef,envir=expl_index)
}

for(f in candidate_files) {
  info <- parse_candidate_path(f)
  candidates <- read.delim(f,stringsAsFactors=FALSE)$gene_id

  if(length(candidates)<MIN_CANDIDATES) {
    message("Skipping ORA for ",info$comparison_dir," dpi",info$dpi," ",info$direction," — only ",length(candidates)," candidates (< ",MIN_CANDIDATES,").")
    next
  }

  expl_key <- paste(info$comparison_dir,info$dpi,sep="||")
  if(!exists(expl_key,envir=expl_index,inherits=FALSE)) {
    message("No matching exploratory file for ",f," — skipping (background unavailable).")
    next
  }

  expl_file <- get(expl_key,envir=expl_index)
  expl <- read.delim(expl_file,stringsAsFactors=FALSE)

  background <- expl$gene_id[
    (!is.na(expl$mean_normalized_control)&expl$mean_normalized_control>=EXPRESSION_FLOOR) |
    (!is.na(expl$mean_normalized_infected)&expl$mean_normalized_infected>=EXPRESSION_FLOOR)
  ]

  out_dir <- file.path(
   enrichment_dir,
   info$comparison_dir,
   paste0("dpi",info$dpi),
   "ORA",
   info$direction
  )
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

  tag <- paste0(species,"_",info$comparison_dir,"_dpi",info$dpi,"_",info$direction)

  go_res <- run_go_ora(candidates,background,tag)
  kegg_res <- run_kegg_ora(candidates,background,tag)

  write_go_results(go_res,out_dir,tag)
  write_kegg_results(kegg_res,out_dir,tag)

  cat("Done ORA:",tag,"\n")
}

cat("ORA complete.\n")

# ---- 8. Build ranked lists -------------------------------------------------

build_ranked_list <- function(expl,id_col,source_col,want_source) {
  sub <- expl[is.finite(expl$stat),c("gene_id","stat")]
  sub <- merge(sub,annot_map[,c("gene_id",id_col,source_col)],by="gene_id")
  names(sub)[names(sub)==id_col] <- "entrez"
  names(sub)[names(sub)==source_col] <- "source"
  sub <- sub[!is.na(sub$entrez)&sub$source==want_source,]
  sub$entrez <- as.character(sub$entrez)
  sub <- sub[!duplicated(sub$entrez),]
  rl <- as.numeric(sub$stat)
  names(rl) <- sub$entrez
  rl <- rl[is.finite(rl)&!is.na(names(rl))&names(rl)!=""]
  sort(rl,decreasing=TRUE)
}

# ---- 9. GO GSEA ------------------------------------------------------------

run_go_gsea <- function(ranked_by_source,tag,out_dir) {
  results <- list()

  for(src in names(ranked_by_source)) {
    rl <- ranked_by_source[[src]]
    if(length(rl)<10) next
    
    rl <- rl[is.finite(rl) & rl != 0]
    rl <- sort(rl, decreasing = TRUE)
    if (length(rl) < 10) next

    orgdb <- if(src=="direct") org.Bt.eg.db else org.Hs.eg.db
#    debug_file <- file.path(out_dir,paste0(tag,"_ranked_go_debug_",src,".rds"))
#    saveRDS(rl,debug_file)

    cat("GSEA GO:",tag,"source =",src,"genes =",length(rl),"range =",paste(range(rl),collapse=" "),"\n")

    res <- tryCatch(
      gseGO(
        geneList=rl,
        OrgDb=orgdb,
        keyType="ENTREZID",
        ont="ALL",
        minGSSize=10,
        maxGSSize=500,
        pAdjustMethod="BH",
        pvalueCutoff=1,
        verbose=FALSE
      ),
      error=function(e){
        message(tag," GSEA GO (",src,") failed: ",conditionMessage(e))
        NULL
      }
    )

    if(!is.null(res) && nrow(as.data.frame(res))>0) results[[src]] <- res
  }

  results
}

# ---- 10. KEGG GSEA ---------------------------------------------------------

run_kegg_gsea <- function(ranked,src_label,tag,out_dir) {
  if(length(ranked)<10) return(NULL)

  ranked <- as.numeric(ranked) |> setNames(names(ranked))
  ranked <- ranked[is.finite(ranked) & ranked != 0 & !is.na(names(ranked)) & names(ranked) != ""]
  ranked <- ranked[names(ranked)%in%TERM2GENE$entrez_id]
  ranked <- ranked[!duplicated(names(ranked))]
  ranked <- sort(ranked,decreasing=TRUE)

#  saveRDS(ranked,file.path(out_dir,paste0(tag,"_ranked_kegg_mapped_debug.rds")))
  cat("KEGG mapped:",length(ranked),"\n")

  if(length(ranked)<10) return(NULL)

  res <- tryCatch(
    GSEA(
      geneList=ranked,
      TERM2GENE=TERM2GENE,
      TERM2NAME=TERM2NAME,
      minGSSize=10,
      maxGSSize=500,
      pAdjustMethod="BH",
      pvalueCutoff=1,
      verbose=FALSE
    ),
    error=function(e){
      message(tag," GSEA KEGG failed: ",conditionMessage(e))
      NULL
    }
  )

  if(is.null(res) || nrow(as.data.frame(res))==0) return(NULL)
  list(result=res,source_label=src_label)
}

# ---- 11. GSEA per comparison -----------------------------------------------

for(f in exploratory_files) {
  info <- parse_exploratory_path(f)
  expl <- read.delim(f,stringsAsFactors=FALSE)

  ranked_go <- list(
    direct=build_ranked_list(expl,"entrez_go","source_go","direct"),
    ortholog=build_ranked_list(expl,"entrez_go","source_go","human_ortholog")
  )
  ranked_go <- ranked_go[lengths(ranked_go)>0]

  tag <- paste0(species,"_",info$comparison_dir,"_dpi",info$dpi)
  out_dir <- file.path(enrichment_dir,info$comparison_dir,paste0("dpi",info$dpi),"GSEA")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

  go_gsea <- run_go_gsea(ranked_go,tag,out_dir)
  write_go_results(go_gsea,out_dir,paste0(tag,"_GSEA"))

  src <- if(species=="cattle") "direct" else "human_ortholog"
  ranked_kegg <- build_ranked_list(expl,"entrez_kegg","source_kegg",src)

#  saveRDS(ranked_kegg,file.path(out_dir,paste0(tag,"_ranked_kegg_debug.rds")))

  kegg_gsea <- run_kegg_gsea(ranked_kegg,src,tag,out_dir)
  write_kegg_results(kegg_gsea,out_dir,paste0(tag,"_GSEA"))

  cat("Done GSEA:",tag,"\n")
}

cat("All enrichment complete. Outputs written to: ",enrichment_dir,"\n")
