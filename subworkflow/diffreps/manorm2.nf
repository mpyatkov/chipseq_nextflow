#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def create_diffreps_channel(row) {
    def meta = [:]
    // 1, Hnf6_Male, Hnf6_Female, G73M03|G73M04|G76M12|G76M13, G73M01|G73M02|G76M10|G76M11, RIPPM, 1000
    meta.num               = row[0].trim()
    meta.treatment_name    = row[1].trim()
    meta.control_name      = row[2].trim()
    meta.treatment_samples = row[3].trim()
    meta.control_samples   = row[4].trim()
    meta.normalization     = row[5].trim()
    meta.window_size       = row[6].trim()
    meta.report_name = "diffReps_${meta.treatment_name}.vs.${meta.control_name}_${meta.normalization}_${meta.window_size}"
    meta.group_name = "${meta.treatment_name}_vs_${meta.control_name}"
    return meta
}

workflow MANORM2 {
    // take:
    // diffreps_config      //parse_configuration_xls.out.diffreps_config
    // sample_labels_config //parse_configuration_xls.out.sample_labels_config
    // fragments_bed6       //bam_count.out.fragments_bed6
    // norm_factors         //calc_norm_factors.out
    // peakcaller_xls       //macs2_callpeak.out.xls
    // peaks_for_manorm2     //macs2_callpeak.out.narrow_bed
    // mm9_chrom_sizes      //mm9_chrom_sizes
    
    take:
    diffreps_config   // diffreps configuration file (MANORM normalization)
    fragments         // bed3 fragments obtained from bam files
    peaks_for_manorm2 //macs2_callpeak.out.narrow_bed
    
    main:

    input_params = diffreps_config
        .splitCsv() 
        .map{it->create_diffreps_channel(it)}
        .branch {
            manorm2: it.normalization =~"MANORM"
            other: true
        }
        
    // only path to fragments
    manorm2_fragments_ch = fragments
        .map{it->it[1]}
        .collect().toList()

    // only path to MACS2/EPIC2 peaks
    manorm2_peaks_ch = peaks_for_manorm2
        .map{it->it[1]}
        .collect().toList()
    
    manorm2_params = input_params.manorm2
        .combine(manorm2_fragments_ch)
        .combine(manorm2_peaks_ch)
        .map{meta,fr,pks ->
            def new_fr = fr.findAll{it =~ "${meta.treatment_samples}|${meta.control_samples}"}
            def new_pks = pks.findAll{it =~ "${meta.treatment_samples}|${meta.control_samples}"}
            return [meta, new_fr, new_pks]
        } | manorm2_create_profile
        
    manorm2_diffexp(manorm2_create_profile.out.profile)

    manorm2_diffexp.out.manorm2_historagms.collect() | aggregate_manorm_pdf
    
    emit:
    profile = manorm2_create_profile.out.profile
    diff_table = manorm2_diffexp.out.diff_table

}

process manorm2_create_profile {
    tag "${meta.group_name}"
    executor 'sge'
    memory '16 GB'
    cpus 1
    time '1h'
    beforeScript 'source $HOME/.bashrc'
    
    input:
    tuple val(meta), path(fragments), path(peaks)
    
    output:
    tuple val(meta), path("${meta.group_name}_profile_bins.xls"), emit: profile
    
    script:
    """
    module load miniconda
    conda activate /projectnb/wax-es/routines/condaenv/manorm2_utils

    peaks=`find . -name "*narrow*" | sort |xargs -n1 basename | paste -s -d ","`
    reads=`find . -name "*fragments.bed" | sort | xargs -n1 basename| paste -s -d ","`
    names=`find . -name "*fragments.bed" | sort | xargs -n1 basename | replace '_fragments.bed' '' | paste -s -d ","`

    profile_bins --peaks=\${peaks} \
        --reads=\${reads} \
        --labs=\${names} \
        -n ${meta.group_name}
    """
}

process manorm2_diffexp {
    tag "${meta.group_name}"

    executor 'local'
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/manorm2_output/${meta.num}_MANORM2_${meta.group_name}/", mode: "copy", pattern: "*.{xlsx,pdf,bed}", overwrite: true
    
    input:
    tuple val(meta), path(profile)

    output:
    tuple val(meta), path("*.xlsx"), emit: diff_table
    tuple val(meta), path("*.bed"), emit: bed_track
    path("*.pdf"), emit: manorm2_historagms
    
    
    script:
    output_prefix = "${meta.group_name}_MANORM2"
    
    """
    module load R
    
    manorm2_diffexpr.R \
        --manorm2_profile ${profile} \
        --control_samples '${meta.control_samples}' \
        --treatment_samples '${meta.treatment_samples}' \
        --control_name ${meta.control_name} \
        --treatment_name ${meta.treatment_name} \
        --output_prefix ${output_prefix} 
    """
}

process aggregate_manorm_pdf {
    tag("${group_name}")
    executor 'local'
    publishDir path: "${params.output_dir}/manorm2_output/", mode: "copy", pattern: "*.pdf", overwrite: true
    
    input:
    path(hist)

    output:
    path("Aggregated_Manorm2_Histograms.pdf")

    script:
    """
    module load poppler
    pdfunite $hist "Aggregated_Manorm2_Histograms.pdf"
    """
}
