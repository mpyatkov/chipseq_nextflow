#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Summary_normalization_factors.xlsx: normalization factors of all three
// methods (diffReps incl. RIPPM, DESeq2, CSAW) for every comparison.
// Waits for all of them, so the workbook appears near the end of a run.
process NORM_FACTORS_REPORT {

    executor 'local'
    container params.rlang_sif
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "*.xlsx", overwrite: true

    input:
    path(diffreps_reports, stageAs: "diffreps/*")  // raw diffReps reports (headers hold the factors)
    path(sample_stats)                             // Norm_Factors.tsv (RIPPM)
    path(fq_num_reads)                             // numreads.csv
    path(deseq2_factors, stageAs: "deseq2/*")      // *_norm_factors.tsv from calc_deseq2
    path(csaw_factors, stageAs: "csaw/*")          // *_norm_factors.tsv from calc_csaw

    output:
    path("*.xlsx")

    script:
    """
    diffreps_output_parser.R --path diffreps \
        --rippm_report ${sample_stats} \
        --fastq_num_reads ${fq_num_reads} \
        --deseq2_factors deseq2 \
        --csaw_factors csaw
    """
}
