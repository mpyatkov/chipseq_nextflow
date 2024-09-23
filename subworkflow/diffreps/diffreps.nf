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

process diffreps {
    tag "${meta.num}"
    executor 'sge'
    cpus 16
    //cache false
    time '2h'
    // echo true

    beforeScript 'source $HOME/.bashrc'
    
    input:
    tuple val(meta), path(TREATMENT_FILES), path(CONTROL_FILES), path(NORM_FILE)
    
    output:
    tuple val(meta) , path("diffReps_*")

    shell:
    template 'diffreps.sh'
}

process diffreps_summary {
    tag "${meta.num}"

    executor 'local'
    // cpus 16
    // echo true
    errorStrategy 'retry'
    maxRetries 2

    
    publishDir path: "${params.output_dir}/diffreps_output/${meta.group_name}/", mode: "copy", pattern: "${output_dir}/*", overwrite: true
    
    input:
    tuple val(meta), path(diffout), path(xls), path(SAMPLE_LABELS)
    path(mm9_chrom_sizes)
    
    output:
    tuple val(meta.num), path("${output_dir}/*")
    tuple val(meta), path("*.xlsx"), emit: full_report
    tuple val(meta.group_name), path("Histograms_AllChr*.pdf"), path("Histograms_noXY*.pdf"), emit: hist_pdf
    tuple val(meta.group_name), path("FDR*.pdf"), emit: fdr_pdf
    tuple val(meta.group_name), path("Bar*.pdf"), emit: bar_pdf
    tuple val(meta.num), val(output_dir), path("*.bb"), emit: diffreps_track
    path("${meta.report_name}"), emit: diffreps_output    

    shell:
    peakcaller="${params.peakcaller}"
    output_dir="${meta.num}_${meta.report_name}_${peakcaller}"
    dataset_label="${params.dataset_label}"
    
    template 'diffreps_summary.sh'
}

process aggregate_diffreps_pdf {
    tag("${group_name}")
    executor 'local'
    publishDir path: "${params.output_dir}/summary/FDR_Barcharts/", mode: "copy", pattern: "${group_name}_*FDR*.pdf", overwrite: true
    publishDir path: "${params.output_dir}/summary/Feature_Barcharts/", mode: "copy", pattern: "${group_name}_*Feature*.pdf", overwrite: true
    
    input:
    tuple val(group_name), path(hist_allchr), path(hist_noxy), path(fdr), path(bar)

    output:
    path("${group_name}_*.pdf")
    tuple val(group_name), path("${group_name}_AllChr_Histogram.pdf"), emit: diffreps_aggregated_histogram_allchr
    tuple val(group_name), path("${group_name}_noXY_Histogram.pdf"), emit: diffreps_aggregated_histogram_noxy
    
    script:
    """
    module load poppler
    allchr_hist_pdf=(\$(find . -name '*Histograms_AllChr*.pdf' | sort))
    pdfunite \${allchr_hist_pdf[@]} "${group_name}_AllChr_Histogram.pdf"

    noXY_hist_pdf=(\$(find . -name '*Histograms_noXY*.pdf' | sort))
    pdfunite \${noXY_hist_pdf[@]} "${group_name}_noXY_Histogram.pdf"

    fdr_pdfs=(\$(find . -name '*FDR*.pdf' | sort))
    pdfunite \${fdr_pdfs[@]} "${group_name}_FDR_Barcharts_${params.peakcaller}.pdf"

    features_pdfs=(\$(find . -name 'Barchart_*.pdf' | sort))
    pdfunite \${features_pdfs[@]} "${group_name}_Feature_Barcharts_${params.peakcaller}.pdf"
    """
}

process collect_diffreps_norm_factors {

    executor 'local'
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "*.xlsx", overwrite: true
    // echo true

    input:
    tuple path(diffreps_output), path(sample_stats), path(fq_num_reads)
    output:
    path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    diffreps_output_parser.R --path "." \
        --rippm_report ${sample_stats} \
        --fastq_num_reads ${fq_num_reads}
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
    fq_num_reads         //table with number of reads in R1.fq files for each sample 

    main:

    input_params = diffreps_config
        .splitCsv() 
        .map{it->create_diffreps_channel(it)}
        .branch {
            diffreps: it.normalization =~"DIFFREPS|RIPPM"
            manorm2: true
        }
    
    diffreps_params = input_params.diffreps
        .combine(
            fragments_bed6
                .map{it->it[1]}
                .collect().toList())
        .map{meta,l ->
            def tl = l.findAll{it =~ meta.treatment_samples}
            def cl = l.findAll{it =~ meta.control_samples}
            return [meta, tl, cl]
        }
    
    diffreps_params_withnorm = diffreps_params.combine(norm_factors)
    diffreps_params_withnorm | diffreps

    dsummary_ch = diffreps_params.map{it -> it[0..-3]} // exclude bed6 files
        .join(diffreps.out)
        .combine(peakcaller_xls
                 .map{it -> it[1]}.collect().toSortedList())
        .map{meta,diffout,xls ->
            //filter xls only participating in comparison
            def new_xls = xls.findAll{it =~ "${meta.treatment_samples}|${meta.control_samples}"} 
            return [meta, diffout, new_xls]}
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
        .combine(norm_factors) 
        .combine(fq_num_reads) | collect_diffreps_norm_factors
    
    emit:
    diffreps_track = diffreps_summary.out.diffreps_track
    full_report = diffreps_summary.out.full_report
    aggregated_histograms_diffreps_allchr = aggregate_diffreps_pdf.out.diffreps_aggregated_histogram_allchr    
    aggregated_histograms_diffreps_noxy = aggregate_diffreps_pdf.out.diffreps_aggregated_histogram_noxy    
}
