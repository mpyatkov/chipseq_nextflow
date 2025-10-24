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

process mumerge {
    tag "${group_name}"
    cpus 8
    beforeScript 'source $HOME/.bashrc'    
    publishDir path: "${params.output_dir}/mumerge_output/", mode: "copy", pattern: "${group_name}_mumerge.bed", overwrite: true

    input:
    tuple val(group_name), path("*")

    output:
    tuple val(group_name), path("${group_name}_mumerge.bed")

    script:
    """
    module load R/${params.rversion}
    run_mumerge.R --files_path ./ --pattern "bed" --ncores $task.cpus --output_file "${group_name}_mumerge.bed"
    """
}

process total_union {

    executor "local"
    // echo true
    
    beforeScript 'source $HOME/.bashrc'
    
    input:
    path("*")
    
    output:
    path("peak_union.bed"), emit: union

    script:
    """
    module load bedtools
    cat * > tmp.bed
    sort -k1,1 -k2,2n tmp.bed > tmp1.bed
    bedtools merge -i tmp1.bed > peak_union.bed
    """
    
    stub:
    """
    cat * > peak_union.bed
    """
}

workflow MUMERGE {
    take:
    diffreps_config
    peakcaller_peaks
    
    main:
    
    input_params = diffreps_config
        .splitCsv() 
        .map{it->parse_config(it)}
        .unique { it.group_name }
        .combine(peakcaller_peaks.map{it -> it[1]}.collect().toSortedList())
    .map{meta,xls ->
            def new_xls = xls.findAll{it =~ "${meta.treatment_samples}|${meta.control_samples}"} 
            return [meta.group_name, new_xls]} 

    input_params | mumerge

    // mumerge.out.count().view{n -> println "Number of mumerge group comparisons: ${n}"}

    total_overlap = null
    if (mumerge.out.count() == 1) {
        total_overlap = mumerge.out
    } else {
        total_union(mumerge.out
                    .map{group_name, mumerge_file  -> [mumerge_file]}
                    .collect())
        
        total_overlap = total_union.out
    }
    
    emit:
    mumerge_peaks = mumerge.out
    mumerge_overlap = total_overlap
}
