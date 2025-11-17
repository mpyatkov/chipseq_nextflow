#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow QUALITY_PCA {

    take:
    diffreps_config      //parse_configuration_xls.out.diffreps_config
    peakcaller_xls       //macs2_callpeak.out.xls
    peakcaller_name

    main:

    pairs_to_compare = diffreps_config
        .map{meta -> ["${meta.treatment_name}_vs_${meta.control_name}",meta.treatment_name, meta.control_name, meta.treatment_samples, meta.control_samples]}
        .groupTuple()        
        .combine(peakcaller_xls
                    .map{it -> it[1]}
                    .collect()
                    .toList()) 
        .map{group_name, tr_name, ctrl_name, tr_samples, ctrl_samples,l -> 
                def xls = l.findAll{it =~ "${tr_samples[0]}|${ctrl_samples[0]}"}
                return [group_name, tr_name[0], ctrl_name[0], tr_samples[0], ctrl_samples[0], peakcaller_name, xls]
            } 
         
    pairs_to_compare | quality_pca_correlation  
    // quality_pca_correlation.out.pdfs.collect() | combine_pdf
    combine_pdf(quality_pca_correlation.out.pdfs.collect(), peakcaller_name)

    // emit:
    // // input_params = quality_pca_correlation.out 
    // combined_pdf = combine_pdf.out 

}

process quality_pca_correlation {
    tag "${group_name}"

    executor 'local'

    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/summary/pileup_data_for_pca/", mode: "copy", pattern: "*.xlsx", overwrite: true

    input:
    tuple val(group_name), val(tr_name), val(ctrl_name), val(tr_samples), val(ctrl_samples), val(peakcaller_name), path(xls_files)
    
    output:
    path("*.pdf"), emit: pdfs
    path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    quality_pca.R --treatment_name $tr_name \
        --control_name $ctrl_name \
        --treatment_samples '$tr_samples'\
        --control_samples '$ctrl_samples'\
        --peakcaller $peakcaller_name

    quality_pca.R --treatment_name $tr_name \
        --control_name $ctrl_name \
        --treatment_samples '$tr_samples'\
        --control_samples '$ctrl_samples' \
        --remove_chrXY \
        --peakcaller $peakcaller_name
    """
}

process combine_pdf {

    executor "local"
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "Summary*.pdf", overwrite: true

    input:
    path(pdfs)
    val(peakcaller_name)

    output:
    path("Summary*.pdf")

    script:

    """
    module load poppler
    all_pdfs=(\$(find . -name '*.pdf' | sort))
    pdfunite \${all_pdfs[@]} "Summary_${peakcaller_name}_PCA_correlation_plots.pdf"

    """
}
