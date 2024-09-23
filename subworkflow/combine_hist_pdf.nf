workflow COMBINE_HIST_PDF {
	take:

	diffreps_noxy
	diffreps_allchr
	manorm2_noxy
	manorm2_allchr

	main:

	  // Combine pdf reports where filtered out X and Y chromosomes for
	  // Manorm2 and diffReps
    mn2_noxy = manorm2_noxy
        .map{meta,rest -> [meta.group_name, rest]}
        .groupTuple() 

    dr_noxy = diffreps_noxy
        .groupTuple() 

    all_noxy = dr_noxy.join(mn2_noxy)  
    combine_noxy(all_noxy)

    // Same but combining pdfs that have all chromosomes reports
    mn2_allchr = manorm2_allchr
        .map{meta,rest -> [meta.group_name, rest]}
        .groupTuple() 

    dr_allchr = diffreps_allchr
        .groupTuple() 

    all_allchr= dr_allchr.join(mn2_allchr)
    combine_allchr(all_allchr)

}



process combine_noxy {
    tag "${group_name}"

    executor "local"
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/summary/Histograms/", mode: "copy", pattern: "*.pdf", overwrite: true
    // echo true

    input:
    tuple val(group_name), path(diffreps_pdfs), path(manorm2_pdf)

    output:
    path("*.pdf")

    script:

    """
    ls -l
    module load poppler
    pdfunite $manorm2_pdf $diffreps_pdfs "${group_name}_Histogram_noXY_${params.peakcaller}.pdf"
    """
}

process combine_allchr {
    tag "${group_name}"

    executor "local"
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/summary/Histograms/", mode: "copy", pattern: "*.pdf", overwrite: true
    // echo true

    input:
    tuple val(group_name), path(diffreps_pdfs), path(manorm2_pdf)

    output:
    path("*.pdf")

    script:

    """
    ls -l
    module load poppler
    pdfunite $manorm2_pdf $diffreps_pdfs "${group_name}_Histogram_AllChrom_${params.peakcaller}.pdf"
    """
}
