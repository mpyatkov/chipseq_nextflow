process fragments_union_overlap_stats {
    tag "${sample_id}"
    
    executor "sge"
    cpus 1
    memory '8 GB'
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/macs2/", mode: "copy", pattern: "${sample_id}_fragments_union_coverage.bed", overwrite: true
    
    input:
    tuple val(sample_id), path(fragments)
    path(union)
    
    output:
    path("${sample_id}_output.txt")
    
    script:
    """
    module load bedtools
    bedtools coverage -a $union -b $fragments > "${sample_id}_fragments_union_coverage.bed"
    
    FRAGMENT_COUNT=\$(cat $fragments | wc -l)
    FRAGMENT_IN_PEAK_COUNT=\$(awk '{n+=\$4;} END {print n;}' ${sample_id}_fragments_union_coverage.bed)
    #FRAGMENT_IN_PEAK_RATIO=\$(awk "BEGIN {printf \"%.4f\", \$FRAGMENT_IN_PEAK_COUNT/\$FRAGMENT_COUNT}")
    FRAGMENT_IN_PEAK_RATIO=\$(awk -v count=\$FRAGMENT_IN_PEAK_COUNT -v total=\$FRAGMENT_COUNT 'BEGIN {printf "%.4f", count/total}')    
    echo "sample_id,nfragments,nfragments_in_peak,ratio" > ${sample_id}_output.txt
    echo "${sample_id},\$FRAGMENT_COUNT,\$FRAGMENT_IN_PEAK_COUNT,\$FRAGMENT_IN_PEAK_RATIO" >> ${sample_id}_output.txt
    """
}

process calc_norm_factors {

    executor "local"
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "*.tsv", overwrite: true
        
    input:
    path(sample_stats)
    
    output:
    path("Norm_Factors.tsv")
    
    script:
    """
    module load R/${params.rversion}
    calculate_norm_factors.R ${sample_stats} "Norm_Factors.tsv"
    """
}

workflow RIPPM_NORM_FACTORS {
    take:
    MUMERGE_UNION
    FRAGMENTS
    
    main:
    fragments_union_overlap_stats(FRAGMENTS, MUMERGE_UNION)

    // TODO: make the output from fragments_union_overlap_stats the coverage files
    // combine them and calc PCA for each sample against mumerge union
    
    sample_stats = fragments_union_overlap_stats.out
        .collectFile(name: "file.csv", keepHeader: true)
    
    sample_stats | calc_norm_factors
    
    emit:
    norm_factors = calc_norm_factors.out
}
