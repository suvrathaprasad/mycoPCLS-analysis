#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * End-to-end RNAseq pipeline with single conda env and no containers.
 *
 * Raw FASTQ -> FastQC -> Trim Galore -> STAR -> featureCounts -> MultiQC
 *                                                      \-> DESeq2 + PCA
 *
 *
 * Usage:
 *   nextflow run rnaseq_pipeline.nf \
 *     -c nextflow.config \
 *     -params-file params_<species>.yaml \
 *     -resume 
 */

params.samplesheet  = null
params.fasta        = null
params.gtf          = null
params.star_index   = null   // optional — built from fasta/gtf if not given
params.metadata     = null
params.species       = null  // "rind" or "schaf" output naming; must match
                              // the uncommented design block in
                              // bin/run_deseq2_pca.R
params.strandedness = 0      // featureCounts -s — confirm strandedness of the experiment before trusting default
params.outdir        = "results"

// -----------------------------------------------------------------------
// 1. Raw read QC
// -----------------------------------------------------------------------
process FASTQC_RAW {
    tag "$sample"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    cache true

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc ${fastq_1} ${fastq_2}
    """
}

// -----------------------------------------------------------------------
// 2. Adapter trimming
// -----------------------------------------------------------------------
process TRIM_GALORE {
    tag "$sample"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'
    cache true

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*trimming_report.txt", emit: reports

    script:
    """
    trim_galore --paired --cores ${task.cpus} ${fastq_1} ${fastq_2}
    """
}

// -----------------------------------------------------------------------
// 3. STAR genome index : built once, shared across every sample this run
// -----------------------------------------------------------------------
process STAR_INDEX {
    tag "star_index"
    publishDir "${params.outdir}/star_index", mode: 'copy'
    label 'process_high'
    cache true

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    """
    mkdir star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100 \\
        --runThreadN ${task.cpus}
    """
}

// -----------------------------------------------------------------------
// 4. Alignment
// -----------------------------------------------------------------------
process STAR_ALIGN {
    tag "$sample"
    publishDir "${params.outdir}/star", mode: 'copy'
    label 'process_high'
    cache true

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)
    path star_index

    output:
    tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample}.Log.final.out", emit: reports

    script:
    """
    STAR \\
        --runMode alignReads \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fastq_1} ${fastq_2} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${sample}. \\
        --runThreadN ${task.cpus}
    """
}

// ----------------------------------------------------------------------------
// 5. Quantification : all samples counted together into one matrix for DESeq2
// ----------------------------------------------------------------------------
process FEATURECOUNTS {
    tag "all_samples"
    publishDir "${params.outdir}/featurecounts", mode: 'copy'
    cache true

    input:
    path bams
    path gtf

    output:
    path "counts.tsv", emit: counts
    path "counts_raw.tsv.summary", emit: summary

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -p --countReadPairs \\
        -s ${params.strandedness} \\
        -a ${gtf} \\
        -o counts_raw.tsv \\
        ${bams}

    # columns of read counts needed for DESeq2
    cut -f1,7- counts_raw.tsv > counts.tsv
    """
}

// ------------------------------------------------------------------------
// 6. MultiQC : aggregate FastQC, Trim Galore, STAR and featureCounts logs
// ------------------------------------------------------------------------
process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cache true

    input:
    path '*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

// ----------------------------------------------------------------------
// 7. DESeq2 + PCA : species-specific design still lives in the R script
// ----------------------------------------------------------------------
process DESEQ2_PCA {
    tag "${params.species}"
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    cache true

    input:
    path counts
    path metadata

    output:
    path "PCA/pca_${params.species}.pdf"
    path "PCA/pca_${params.species}_technical_replicates.pdf"
    path "PCA/dds_${params.species}_collapsed.rds"
    path "PCA/collapsed_metadata_${params.species}.csv"
    path "SECTION_6/*/dds_${params.species}_*_vs_*.rds"
    path "SECTION_6/*/${params.species}_*_DEGs.tsv"
    path "SECTION_6/*/${params.species}_*_upDEGs.tsv"
    path "SECTION_6/*/${params.species}_*_downDEGs.tsv"
    path "SECTION_7/*/${params.species}_*_DESeq2_exploratory.tsv"
    path "SECTION_7/*/${params.species}_*_candidate_up.tsv"
    path "SECTION_7/*/${params.species}_*_candidate_down.tsv"
    path "SECTION_7/*/${params.species}_*_infection_LFC_matrix.tsv"
    path "SECTION_7/*/${params.species}_*_candidate_genes_all_DPI.tsv"
    path "SECTION_7/*/${params.species}_*_temporal_summary.tsv"
    path "SECTION_7/*/dds_${params.species}_*_vs_*.rds"

    script:
    """
    run_deseq2_pca.R ${counts} ${metadata} ${params.species} .
    """
}

workflow {

    reads_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }

    fasta_ch    = Channel.fromPath(params.fasta, checkIfExists: true)
    gtf_ch      = Channel.fromPath(params.gtf, checkIfExists: true)
    metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)

    // ---- 1. raw QC --------------------------------------------------------
    FASTQC_RAW(reads_ch)

    // ---- 2. trimming --------------------------------------------------------
    TRIM_GALORE(reads_ch)

    // ---- 3. genome index : build unless one was already supplied ----------
    if (params.star_index) {
        index_ch = Channel.fromPath(params.star_index, type: 'dir', checkIfExists: true)
    } else {
        STAR_INDEX(fasta_ch, gtf_ch)
        index_ch = STAR_INDEX.out.index
    }

    // ---- 4. alignment : same index reused across every sample -------------
    STAR_ALIGN(TRIM_GALORE.out.trimmed_reads, index_ch.first())

    // ---- 5. quantification : one matrix across all samples ----------------
    // Sorted by sample name before collecting: STAR_ALIGN tasks finish in
    // whatever order they happen to (varies run to run, especially with
    // partial cache hits) and an unsorted .collect() would give
    // FEATURECOUNTS a differently-ordered input list each time which
    // changes its task hash and breaks -resume caching even when every
    // underlying BAM is identical.
    all_bams_ch = STAR_ALIGN.out.bam
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { list -> list.collect { sample, bam -> bam } }
    FEATURECOUNTS(all_bams_ch, gtf_ch)

    // ---- 6. DESeq2 + PCA ------------------------------------------------
    DESEQ2_PCA(FEATURECOUNTS.out.counts, metadata_ch)

    // ---- 7. MultiQC : everything with a log/report to aggregate -----------
    // Same reasoning as above, sort by filename so the input list order
    // (and therefore this task's cache hash) is stable across resumes.
    multiqc_inputs = FASTQC_RAW.out.reports
        .mix(TRIM_GALORE.out.reports)
        .mix(STAR_ALIGN.out.reports)
        .mix(FEATURECOUNTS.out.summary)
        .flatten()
        .toSortedList { a, b -> a.name <=> b.name }

    MULTIQC(multiqc_inputs)
}
