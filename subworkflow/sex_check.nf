workflow SEX_CHECK {
    take:
    all_bams
    labels
    sex_marker_genes
    norm_factors
    
    main:
    only_bam_files = all_bams
        .collect{it -> it[2]}

    COUNT_SEX_MARKERS(only_bam_files, sex_marker_genes)

    GENERATE_REPORT(
        labels,
        norm_factors,
        COUNT_SEX_MARKERS.out
    )
    emit:
    report = GENERATE_REPORT.out
}

process COUNT_SEX_MARKERS {

    executor 'sge'
    cpus 16
    memory '32 GB'
    time '4h'

    beforeScript 'source $HOME/.bashrc'
    // echo true
    
    input:
    path(bam_files)
    path(sex_marker_genes)

    output:
    path("sex_markers.counts")
    
    script:
    """
    module load miniconda
    conda activate /projectnb/wax-es/routines/condaenv/featureCounts2
    
    featureCounts -F SAF -a $sex_marker_genes -p --countReadPairs -Q 30 -T $task.cpus \
        -o sex_markers.counts *.bam
    """
}


process GENERATE_REPORT {
    tag "${group_name}"
    cpus 1
    beforeScript 'source $HOME/.bashrc'
    executor 'local'
    // publishDir path: "${params.output_dir}/mumerge_output/", mode: "copy", pattern: "${group_name}_mumerge.bed", overwrite: true
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "*.xlsx", overwrite: true
    
    input:
    path(labels)
    path(norm_factors)
    path(counts)

    output:
    path("*.xlsx")
    
    script:
    """
    module load R/${params.rversion}
    sex_check_atac.R \
        --counts $counts \
        --norm_factors $norm_factors \
        --samples $labels \
        --output_file sex_confirmation.xlsx \
        --warnings_file sex_check_warnings.txt
    """
}


