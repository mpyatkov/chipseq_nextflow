workflow CREATE_BYCOMPARISON_TRACKS {
    take:
    xlsx_reports_ch
    bed_reports_ch
    chrom_sizes
    
    main:
    xlsx2track(xlsx_reports_ch, chrom_sizes)
    bed2track(bed_reports_ch, chrom_sizes)
    
    emit:
    tracks = xlsx2track.out.track
        .mix(bed2track.out.track)

}

process bed2track {
    tag "${group_name}"

    beforeScript 'source $HOME/.bashrc'

    executor 'local'
    cpus 1

    input:
    tuple val(group_name), path(bedfile)
    path(chrom_sizes)
    
    output:
    tuple val(group_name), path("*.bb"), emit: track

    script:
    
    """
    module load bedtools
    
    # ignore chrM,random and header
    if [[ \$(wc -l <*.bed) -ge 2 ]]; then
       for bed in `find . -name "*.bed"`; do
           cat \${bed} | grep -vE "track|chrM|random|chrUn" > tmp.bed
           bedtools sort -i tmp.bed > tmp.sorted.bed
           bedClip tmp.sorted.bed ${chrom_sizes} tmp.sorted.clipped.bed
           bedToBigBed -allow1bpOverlap tmp.sorted.clipped.bed ${chrom_sizes} "${group_name}_MUMERGE_UCSC_track.bb"
           rm tmp*.bed
       done
    fi

    """
}

process xlsx2track {
    tag "${meta.num}_${meta.method}_${meta.window_size}"

    container params.rlang_sif

    executor 'local'
    cpus 1

    input:
    tuple val(meta), path(report)
    path(chrom_sizes)
    
    output:
    tuple val(meta.group_name), path("*.bb"), emit: track

    script:
    
    """
    
    create_by_comparison_track.R --report_path ${report} \
        --track_name "${meta.report_name}_UCSC_track.bed"
    
    # ignore chrM,random and header
    if [[ \$(wc -l <*.bed) -ge 2 ]]; then
       for bed in `find . -name "*_UCSC_track*.bed"`; do
           cat \${bed} | grep -vE "track|chrM|random|chrUn" > tmp.bed
           bedtools sort -i tmp.bed > tmp.sorted.bed
           bedClip tmp.sorted.bed ${chrom_sizes} tmp.sorted.clipped.bed
           bedToBigBed -allow1bpOverlap tmp.sorted.clipped.bed ${chrom_sizes} "${meta.report_name}_UCSC_track.bb"
           rm tmp*.bed
       done
    fi

    """

    
}


