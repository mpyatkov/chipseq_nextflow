#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process diffreps {
    tag "${NUM}"
    executor 'sge'
    cpus 16
    //cache false
    time '2h'
    // echo true

    beforeScript 'source $HOME/.bashrc'
    
    input:
    tuple val(NUM), val(TREATMENT_NAME), val(CONTROL_NAME), val(TREATMENT_SAMPLES), val(CONTROL_SAMPLES),
        val(NORMALIZATION), val(WINDOW_SIZE), path(TREATMENT_FILES), path(CONTROL_FILES), path(NORM_FILE)
    output:
    tuple val(NUM), val(report_name), path("diffReps_*")

    shell:
    report_name="diffReps_${TREATMENT_NAME}.vs.${CONTROL_NAME}_${NORMALIZATION}_${WINDOW_SIZE}"
    
    template 'diffreps.sh'
}

process diffreps_summary {
    tag "${NUM}"

    executor 'local'
    // cpus 16
    // echo true
    
    publishDir path: "${params.output_dir}/diffreps_output/", mode: "copy", pattern: "${output_dir}/*", overwrite: true
    
    input:
    tuple val(NUM), val(TREATMENT_NAME), val(CONTROL_NAME), val(TREATMENT_SAMPLES), val(CONTROL_SAMPLES),
        val(NORMALIZATION), val(WINDOW_SIZE), val(report_name), path(diffout), path(xls), path(SAMPLE_LABELS)
    path(mm9_chrom_sizes)

    
    output:
    tuple val(NUM), path("${output_dir}/*")
    tuple val(report_name), path("Hist*.pdf"), emit: hist_pdf
    tuple val(report_name), path("FDR*.pdf"), emit: fdr_pdf
    tuple val(report_name), path("Bar*.pdf"), emit: bar_pdf
    tuple val(NUM), val(output_dir), path("*.bb"), emit: diffreps_track
    path("${report_name}"), emit: diffreps_output    

    shell:
    peakcaller="${params.peakcaller}"
    output_dir="${NUM}_${report_name}_${peakcaller}"
    
    dataset_label="${params.dataset_label}"
    template 'diffreps_summary.sh'
}

process aggregate_diffreps_pdf {
    tag("${group_name}")
    executor 'local'
    publishDir path: "${params.output_dir}/diffreps_output/aggregated_pdfs/${group_name}", mode: "copy", pattern: "${group_name}_*.pdf", overwrite: true
    
    input:
    tuple val(group_name), path(hist), path(fdr), path(bar)
    output:
    path("${group_name}_*.pdf")
    script:
    
    """
    module load poppler
    pdfunite $hist "${group_name}_Histogram_Barcharts_${params.peakcaller}.pdf"
    pdfunite $fdr "${group_name}_FDR_Barcharts_${params.peakcaller}.pdf"
    pdfunite $bar "${group_name}_Barcharts_${params.peakcaller}.pdf"
    """
}

process collect_diffreps_norm_factors {

    executor 'local'
    publishDir path: "${params.output_dir}/diffreps_output/", mode: "copy", pattern: "*.xlsx", overwrite: true
    echo true

    input:
    tuple path(diffreps_output), path(sample_stats)
    output:
    path("*.xlsx")

    script:
    """
    module load R
    diffreps_output_parser.R --path "." \
        --rippm_report ${sample_stats}
    """
}

workflow DIFFREPS {
    take:
    diffreps_config      //parse_configuration_xls.out.diffreps_config
    sample_labels_config //parse_configuration_xls.out.sample_labels_config
    fragments_bed6       //bam_count.out.fragments_bed6
    norm_factors         //calc_norm_factors.out
    peakcaller_xls       //macs2_callpeak.out.xls
    mm9_chrom_sizes      //mm9_chrom_sizes

    main:
    diffreps_params = diffreps_config.splitCsv()
        .combine(
            fragments_bed6
                .map{it->it[1]}
                .collect().toList())
        .map{num,tn,cn,tsmp,csmp,norm,w,l->
            def tl = l.findAll{it =~ tsmp.trim()} 
            def cl = l.findAll{it =~ csmp.trim()}
            return [num.trim(),tn.trim(),cn.trim(),tsmp.trim(),csmp.trim(),norm.trim(),w.trim(),tl,cl]}
        
    diffreps_norm_params = diffreps_params
        .map{it->[it[0], it[3].trim(), it[4].trim()]}

    diffreps_params_withnorm = diffreps_params.combine(norm_factors)
    diffreps_params_withnorm | diffreps

    dsummary_ch = diffreps_params.map{it -> it[0..-3]} // exclude bed6 files
        .join(diffreps.out)
        .combine(peakcaller_xls
                 .map{it -> it[1]}.collect().toList())
        .map{it ->
            //filter xls only participating in comparison
            def new_xls = it[-1].findAll{it =~ "${it[3]}|${it[4]}"} 
            def new_out = it[0..-2] << new_xls
            return new_out}
        .combine(sample_labels_config)
    
    diffreps_summary(dsummary_ch, mm9_chrom_sizes)
    
    //aggregate diffpres pdfs
    hist_pdfs = diffreps_summary.out.hist_pdf
        .groupTuple()
    fdr_pdfs = diffreps_summary.out.fdr_pdf
        .groupTuple()
    bar_pdfs = diffreps_summary.out.bar_pdf
        .groupTuple()
    
    all_pdfs = hist_pdfs
        .join(fdr_pdfs)
        .join(bar_pdfs) | aggregate_diffreps_pdf
    
    // -- collect diffreps norm factors (publishDir them)
    diffreps_summary.out.diffreps_output
        .collect()
        .toList()
        .combine(norm_factors) | collect_diffreps_norm_factors
    
    
    emit:
    diffreps_track = diffreps_summary.out.diffreps_track
}
