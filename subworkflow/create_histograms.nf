workflow CREATE_HISTOGRAMS {
    take:
    histograms_ch


    main:

    histograms_ch | create_individual_histograms

    hist_by_comparisons_ch = create_individual_histograms.out.histogram
        .groupTuple() 

    hist_by_comparisons_ch | combine_histograms
    
}

process create_individual_histograms {
    tag "${meta.num}_${meta.method}_${meta.window_size}"

    container params.rlang_sif

    executor 'local'
    cpus 1
    
    input:
    tuple val(meta), path(report)
    
    output:
    tuple val(meta.group_name), path(hist_name), emit: histogram

    script:
    hist_name="${meta.report_name}_histogram.pdf"
    """
    
    create_histogram.R --report_path ${report} \
        --control_name ${meta.control_name} \
        --treatment_name ${meta.treatment_name} \
        --method ${meta.method} \
        --window_size ${meta.window_size} \
        --control_samples "${meta.control_samples}" \
        --treatment_samples "${meta.treatment_samples}" \
        --compar_number ${meta.num} \
        --output_name ${hist_name}
    """
}


process combine_histograms {
    tag "${group_name}"

    executor "local"
    container params.containers.poppler
    publishDir path: "${params.output_dir}/summary/Histograms/", mode: "copy", pattern: "*combined_histograms.pdf", overwrite: true
    // echo true

    input:
    tuple val(group_name), path(diffreps_pdfs)

    output:
    path("*combined_histograms.pdf")

    script:

    """
    sorted_pdfs=(\$(find . -name '*.pdf' | sort))
    pdfunite \${sorted_pdfs[@]} "${group_name}_combined_histograms.pdf"
    """
}
