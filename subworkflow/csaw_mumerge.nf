#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process calc_csaw {
    tag "${meta.group_name}"
    cpus 4
    time '4h'
    memory '16 GB'
    
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/csaw_output/", mode: "copy", pattern: "*.xlsx", overwrite: true
        
    input:
    tuple val(meta), path("*"), path(MUMERGE)

    output:
    tuple val(meta), path("*.xlsx")

    script:
    output_prefix="${meta.group_name}_MuMerge_CSAW"

    """
    module load R/${params.rversion}
 
    csaw_diffpeak_analysis.R --control_samples "${meta.control_samples}" \
        --treatment_samples "${meta.treatment_samples}" \
        --library "${meta.library}" \
        --mumerge "${MUMERGE}" \
        --output_prefix ${output_prefix}
    """
}

workflow CSAW_MUMERGE {
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

    
    calc_csaw(input_params)
    
    emit:
    full_report = calc_csaw.out
}
