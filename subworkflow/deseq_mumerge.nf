#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def parse_config(row) {
    def meta = [:]
    meta.treatment_name    = row[1].trim()
    meta.control_name      = row[2].trim()
    meta.treatment_samples = row[3].trim()
    meta.control_samples   = row[4].trim()
    meta.group_name = "${meta.treatment_name}_vs_${meta.control_name}"

    return meta
}

// Counting of reads/fragments inside MUMERGE peaks
process calc_coverage {
    tag "${meta.group_name}"
    cpus 4
    time '2h'
    beforeScript 'source $HOME/.bashrc'    
        
    input:
    tuple val(meta), path("*"), path(MUMERGE)
    
    output:
    tuple val(meta), path("${meta.group_name}_coverage.tsv")

    script:
    featurecounts_options = null
    if (meta.library == "paired") {
        featurecounts_options = "-p --countReadPairs -B -C"
    } else {
        featurecounts_options = "-s 0"
    }
    
    """
    module load subread
    TRBAM=\$(find . -name "*.bam" | grep -E "${meta.treatment_samples}" | xargs -n1 basename | sort | paste -s -d " ")
    CTRLBAM=\$(find . -name "*.bam" | grep -E "${meta.control_samples}" | xargs -n1 basename | sort | paste -s -d " ")
    awk 'OFS="\t" {print \$1"_"\$2"_"\$3, \$1, \$2+1, \$3, "."}' ${MUMERGE} > mumerge_peaks.saf
    featureCounts -T $task.cpus -F SAF ${featurecounts_options} -a mumerge_peaks.saf -o ${meta.group_name}_coverage.tsv \${TRBAM} \${CTRLBAM}
    """
}

process calc_deseq2 {
    tag "${meta.group_name}"
    // cpus 1
    // time '1h'
    executor "local"
    
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/deseq2_output/", mode: "copy", pattern: "*.xlsx", overwrite: true

        
    input:
    tuple val(meta), path(coverage)

    output:
    tuple val(meta), path("*.xlsx")

    script:
    output_prefix="${meta.group_name}_MuMerge_DEseq2"

    """
    module load R/${params.rversion}

    deseq2_diffpeak_analysis.R --coverage_file ${coverage} \
        --control_samples "${meta.control_samples}" \
        --treatment_samples "${meta.treatment_samples}" \
        --output_prefix ${output_prefix}
    """
}

workflow DESEQ_MUMERGE {
    take:
    diffreps_config
    bams
    mumerge_peaks
    
    main:

    bambai_only = bams
        .map { id, p, bam, bai -> [p, [bam, bai]] }
        .groupTuple(by: 0)
        .map { p, files -> [p, files.flatten()] }
    
    input_params = diffreps_config
        .splitCsv() 
        .map{it->parse_config(it)}
        .unique { it.group_name }
        .combine(bambai_only)
        .map{meta,lib,bambai ->
            def new_meta = meta
            new_meta.library = lib
            return [new_meta, bambai]} 
        .map{meta,bambai ->
            def new_bambai = bambai.findAll{it =~ "${meta.treatment_samples}|${meta.control_samples}"}
            return [meta.group_name, meta, new_bambai]}
        .combine(mumerge_peaks,by: 0)
        .map{group_name,meta, bambai, mumerge_peaks -> [meta,bambai,mumerge_peaks]} 

    calc_coverage(input_params) | calc_deseq2
    
    emit:
    deseq_results = calc_deseq2.out
}
