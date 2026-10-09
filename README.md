# mycoPCLS RNAseq Functional Annotation Pipeline

An end-to-end Nextflow pipeline for RNAseq analysis of precision-cut lung slices (PCLS) from cattle (Rind) and sheep (Schaf) infected with Mycobacterium bovis. The pipeline performs quality control, alignment, quantification, differential expression analysis, functional enrichment, and generates publication-ready plots and tables.

## Overview

This pipeline processes raw FASTQ files through the following stages:
1. **FASTQC_RAW** - Raw read quality control
2. **TRIM_GALORE** - Adapter trimming
3. **STAR_INDEX** - Genome index generation (optional)
4. **STAR_ALIGN** - Read alignment to reference genome
5. **FEATURECOUNTS** - Gene-level quantification
6. **MULTIQC** - Aggregate QC reports
7. **PCA_DESEQ2** - Principal component analysis and differential expression (DESeq2)
8. **ENRICHMENT** - Functional enrichment analysis (GO/KEGG/GSEA)
9. **TABLES** - Flat, browsable TSV summaries
10. **PLOTS** - Publication-ready visualizations
11. **CUSTOM_PLOTS** - User-supplied gene list plots (optional)

## Key Features

- **Single conda environment** - No containers required
- **Nextflow orchestration** - Robust, reproducible workflow with `-resume` capability
- **R-only analysis scripts** - All statistical analysis performed in R
- **Species-specific handling** - Annotation strategy per species: currently cattle hybrid (direct + ortholog fallback) and sheep human ortholog-only; moving to human ortholog-only for all species (planned)
- **Technical replicate-aware** - Properly handles the experimental design (1 biological unit per condition×dpi cell)
- **Comprehensive outputs** - From raw QC to publication-ready figures and tables

## Installation

### Prerequisites
- [Nextflow](https://www.nextflow.io/) (included in conda environment)
- Java JDK 17+ (included in conda environment)
- [Mamba](https://mamba.readthedocs.io/) or [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Setup
1. Clone this repository
2. Create the reference folders and place the genome FASTA and GTF there (`ref/` is not tracked in git), e.g. `ref/rind/`, `ref/schaf/`
3. Create the conda environment:
   ```bash
   mamba env create -f rnaseq-pipeline-env.yml
   ```
4. Activate the environment:
   ```bash
   conda activate rnaseq-pipeline
   ```

## Usage

### Basic Execution
Run the pipeline for either cattle or sheep:

```bash
# For cattle (Rind)
nextflow run rnaseq_pipeline.nf \
  -c nextflow.config \
  -params-file params_rind.yaml \
  -resume

# For sheep (Schaf)
nextflow run rnaseq_pipeline.nf \
  -c nextflow.config \
  -params-file params_schaf.yaml \
  -resume
```

### Parameters
Key parameters can be overridden via command line or in the YAML parameter files:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--samplesheet` | CSV with sample information (sample, fastq_1, fastq_2) | `samplesheet_rind.csv` / `samplesheet_schaf.csv` |
| `--fasta` | Reference genome FASTA | `ref/rind/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa` / `ref/schaf/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa` |
| `--gtf` | Reference genome annotation GTF | `ref/rind/Bos_taurus.ARS-UCD2.0.116.gtf` / `ref/schaf/Ovis_aries.ARS-UI_Ramb_v3.0.116.gtf` |
| `--star_index` | Pre-built STAR index (optional) | Built automatically if not provided |
| `--metadata` | Experimental design metadata | `metadata_rind.csv` / `metadata_schaf.csv` |
| `--species` | Species identifier (must match R scripts) | `"cattle"` / `"sheep"` |
| `--strandedness` | RNA-seq strandedness (0=unstranded, 1=stranded, 2=reverse) | `0` |
| `--outdir` | Base output directory | `results` (creates `results/rind/` or `results/schaf/`) |
| `--annotation_cache` | Not yet configurable: the pipeline currently reads the cache from `ref/annotation_cache` (fixed in `rnaseq_pipeline.nf`) | `${projectDir}/ref/annotation_cache` |
| `--gene_list` | Optional TSV with gene IDs for custom plots | `null` |

### Example with Overrides
```bash
nextflow run rnaseq_pipeline.nf \
  -c nextflow.config \
  -params-file params_rind.yaml \
  --star_index /path/to/prebuilt/star_index \
  --outdir /path/to/output \
  -resume
```

## Pipeline Stages

### 1. Quality Control
- **FASTQC_RAW**: Per-sample fastQC on raw reads
- **TRIM_GALORE**: Adapter and quality trimming using Trim Galore!

### 2. Alignment & Quantification
- **STAR_INDEX**: Genome index generation (if needed)
- **STAR_ALIGN**: Splice-aware alignment with STAR
- **FEATURECOUNTS**: Gene-level quantification using Subread featureCounts
- **MULTIQC**: Aggregated QC report

### 3. Differential Expression
- **PCA_DESEQ2**: 
  - PCA visualization (both regular and technical replicate-focused)
  - DESeq2 analysis with two-stage design:
    1. Estimate dispersion on `~ dpi + condition` (across time course)
    2. Refit to saturated `~0 + dpi_condition` design
    3. Extract Wald test statistics for day-specific contrasts
  - Outputs: PCA plots, DESeq2 objects, exploratory results, and candidate gene lists

### 4. Functional Enrichment
- **ENRICHMENT**: 
  - ORA (Over-Representation Analysis) and GSEA (Gene Set Enrichment Analysis)
  - GO (Gene Ontology) and KEGG pathways
  - Species-specific annotation strategy:
    - Cattle: Hybrid approach (direct `org.Bt.eg.db` + human ortholog fallback); being replaced by human ortholog-only (planned)
    - Sheep: Human ortholog-only approach (`org.Hs.eg.db`)
  - Uses precomputed annotation cache to avoid internet calls on compute nodes

### 5. Summarization & Visualization
- **TABLES**: Flat TSV files for easy browsing:
  1. `1_master_candidate_genes` - All candidates across comparisons
  2. `2_comparison_summary` - Counts per comparison/dpi
  3. `3_master_enrichment_ora` - All ORA hits
  4. `4_gene_pathway_crosswalk` - Gene-to-pathway mappings
  5. `5_temporal_trajectories` - log2FC across timepoints
  6. `6_gsea_summary` - All GSEA results
- **PLOTS**: Publication-ready visualizations (12 plot types; underlying data saved as RDS in `plots_data/` for later editing with `helpeR_scripts/modify_plots.R`). Main groups:
  1. PCA plots
  2. Heatmaps (top candidates)
  3. Volcano plots
  4. Temporal trajectories
  5. Candidate overlap Venn (cattle only)
  6. enrichment plots (emapplot, etc.)
  7. Cross-comparison overview plots (candidates over time, trajectory overlay, etc.)
- **CUSTOM_PLOTS** (optional): User-defined gene list analysis

## Output Structure
```
results/
└── [rind|schaf]/
    ├── fastqc_raw/           # FastQC reports (raw reads)
    ├── trim_galore/          # Trim Galore! outputs
    ├── star_index/           # STAR genome index (if built)
    ├── star/                 # STAR alignments
    ├── featurecounts/        # Gene count matrices
    ├── multiqc/              # MultiQC report
    ├── deseq2/               # DESeq2 and enrichment results
    │   ├── PCA/              # PCA plots and objects
    │   ├── enrichment/       # Functional enrichment results
    │   │   ├── [comparison]/dpi[N]/ORA/{up,down}/  # ORA results
    │   │   └── [comparison]/dpi[N]/GSEA/       # GSEA results
    │   └── ...               # DESeq2 exploratory outputs
    ├── tables/               # Summary TSV files
    ├── plots/                # Publication-ready figures
    └── custom_gene_plots/    # User gene list plots (if --gene_list provided)
```

## Annotation Cache
The pipeline uses an annotation cache to avoid live BioMart/KEGG calls on HPC compute nodes:
1. Run `fetch_annotation_cache.R` once per species (cattle, sheep, goat) into the **same** cache directory. Usage: `Rscript helpeR_scripts/fetch_annotation_cache.R <species> <cache_dir> [gene_ids_file]`. The optional third argument restricts the query to a list of gene IDs; without it, the whole genome is fetched:
   ```bash
   Rscript helpeR_scripts/fetch_annotation_cache.R cattle /path/to/cache
   Rscript helpeR_scripts/fetch_annotation_cache.R sheep /path/to/cache
   Rscript helpeR_scripts/fetch_annotation_cache.R goat /path/to/cache
   ```
2. The pipeline uses the cache directory specified in parameters
3. Cache includes:
   - Species-specific annotation maps (Ensembl ↔ Entrez ↔ gene symbols)
   - Shared KEGG pathway↔gene tables (`kegg_hsa.rds`)
   - Failed BioMart ID logs for troubleshooting

## Experimental Design Notes
These notes describe the PCLS pilot experiment. The main multi-strain experiment and goat are planned.

- **Technical replicates only**: 1 biological unit per condition×dpi cell
- **Two sequencing runs per animal/dpi/condition**: Treated as technical replicates (not same library)
- **No `collapseReplicates()`**: Runs kept as separate columns in DESeq2 analysis
- **Dispersion borrowing**: Estimated across timepoints within condition pairs
- **Species argument**: English ("cattle"/"sheep") in R scripts for publication-ready labels

## Dependencies
See `rnaseq-pipeline-env.yml` for the complete conda environment specification, including:
- **Pipeline orchestration**: Nextflow, OpenJDK
- **Preprocessing**: FastQC, Trim Galore!, STAR, Subread, MultiQC
- **Differential expression**: Bioconductor DESeq2, PCAtools
- **Enrichment**: Bioconductor clusterProfiler, enrichplot, OrgDb packages
- **Visualization**: ggplot2, ggrepel, pheatmap, ComplexHeatmap, RColorBrewer, enhancedvolcano, ggVennDiagram

## Related Scripts
Located in `bin/` and `helpeR_scripts/`:
- `bin/run_pca_deseq2.R` - DESeq2 and PCA analysis
- `bin/run_enrichment.R` - Functional enrichment (ORA/GSEA)
- `bin/run_tables.R` - Summary table generation
- `bin/run_plots.R` - Publication-ready plotting
- `bin/run_plots_custom.R` - Custom gene list plots
- `helpeR_scripts/fetch_annotation_cache.R` - Annotation cache generation
- `helpeR_scripts/annotation_coverage_check.R` - Annotation coverage check
- `helpeR_scripts/modify_plots.R` - Re-style plots from saved `plots_data/` RDS files
- `helpeR_scripts/run_cross_species_comparison.R` - Cross-species pathway analysis

## Author
Developed by Suvi (Suvratha Jayaprasad), Postdoc at Friedrich-Loeffler-Institut (FLI) Jena.
