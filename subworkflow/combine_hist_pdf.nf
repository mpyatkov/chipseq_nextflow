workflow COMBINE_HIST_PDF {
	take:

	diffreps_noxy
	diffreps_allchr

	main:

    dr_noxy = diffreps_noxy
        .groupTuple() 

    all_noxy = dr_noxy
    combine_noxy(all_noxy)

    // Same but combining pdfs that have all chromosomes reports
    dr_allchr = diffreps_allchr
        .groupTuple() 

    all_allchr= dr_allchr
    combine_allchr(all_allchr)
}



process combine_noxy {
    tag "${group_name}"

    executor "local"
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/summary/Histograms/", mode: "copy", pattern: "*.pdf", overwrite: true
    // echo true

    input:
    tuple val(group_name), path(diffreps_pdfs)

    output:
    path("*.pdf")

    script:

    """
    ls -l
    module load poppler
    pdfunite $diffreps_pdfs "${group_name}_Histogram_noXY_${params.peakcaller}.pdf"
    """
}

process combine_allchr {
    tag "${group_name}"

    executor "local"
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/summary/Histograms/", mode: "copy", pattern: "*.pdf", overwrite: true
    // echo true

    input:
    tuple val(group_name), path(diffreps_pdfs)

    output:
    path("*.pdf")

    script:

    """
    ls -l
    module load poppler
    pdfunite $diffreps_pdfs "${group_name}_Histogram_AllChrom_${params.peakcaller}.pdf"
    """
}
