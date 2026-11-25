workflow CREATE_BIGWIG_FILES {
    take:
    norm_factors
    fragments
    labels
    chrom_sizes
    
    main:

    // individual BigWig files
    sid_normfact_ch = norm_factors.splitCsv(sep: "\t")
        .map{it->[it[0],it[4]]}
    
    sid_fr_norm = fragments
        .join(sid_normfact_ch)

    create_bigwig_files(sid_fr_norm, chrom_sizes)


    // Combine VigWig files from one group into one track
    bigwig_group_ch = labels
        .splitCsv()
        .join(create_bigwig_files.out)
        .map{sid,n1,id1,n2,r,g,b,pth -> 
             [n2, sid, "${r}__${g}__${b}", pth]
        }
        .groupTuple() 
        .map{group_name, sids, colors, paths ->
            def new_sids = sids.collect{it -> it.toString()}.sort().join("__") 
            return [group_name, new_sids, colors[1], paths]
        } 

    combine_bigwig_tracks(bigwig_group_ch, chrom_sizes)

    emit:
    individual = create_bigwig_files.out
    grouped = combine_bigwig_tracks.out
}


process create_bigwig_files {

    tag "${sample_id}"
    executor 'sge'
    cpus 1
    time '1h'
    memory '16 GB'

    container params.tools_sif

    storeDir "${params.chipseq_bam_cache}/${sample_id}/bigwig/"
    
    input:
    tuple val(sample_id), path(fragments), val(norm_factor)
    path(chrom_sizes)
    
    output:
    tuple val(sample_id), path("*.bw")
    
    script:
    // Create the storeDir folder here, on the Nextflow node: if the task creates it,
    // this node may not see it yet when checking the outputs ("Missing output file").
    file("${params.chipseq_bam_cache}/${sample_id}/bigwig").mkdirs()
    """

    bedtools sort -i ${fragments} > tmp.sorted.bed
    bedtools genomecov -bg -i tmp.sorted.bed -g ${chrom_sizes} > tmp.bedGraph
    sed -i '/track/d;/random/d;/chrM/d;/chrUn/d' tmp.bedGraph
    awk -v OFS='\t' -v norm_factor="${norm_factor}" '{print \$1, \$2, \$3, \$4 * norm_factor}' "tmp.bedGraph" > "tmp.norm.bedGraph"
    bedGraphToBigWig "tmp.norm.bedGraph" ${chrom_sizes} "${sample_id}_RiPPM_norm.bw"
    """
}

process combine_bigwig_tracks {
    tag "${group_name}"

    executor 'sge'
    cpus 8
    memory '16 GB'
    time '2h'

    container params.containers.deeptools
    // echo true

    input:
    tuple val(group_name), val(samples_string), val(color), path(sample_paths)
    path(chrom_sizes)
    
    output:
    tuple val(group_name), val(samples_string), val(color), path("${group_name}.bw")

    script:
    """
    ## for wiggletools
    bigwigAverage -p $task.cpus -b *.bw -o "${group_name}.bw"
    """
}
