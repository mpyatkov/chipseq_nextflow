workflow CREATE_TRACKLINES {
    take:
    individual_ch
    group_ch
    comparison_ch
    sample_labels_config_ch
    
    main:

    // INDIVIDUAL (SAMPLE SPECIFIC)
    individual_ch
        .map{it->[it[0], it[1].getName()]}
        .collectFile{item -> item.join(",")+'\n'}
        .combine(sample_labels_config_ch) |
        sample_specific_tracks

    
    // GROUPED BIGWIG 
    group_tracklines_ch = group_ch
        .map{group_name,samples_str,color,pth -> [group_name, samples_str, color, pth.getName()]}
        .collectFile{item -> item.join(",")+'\n'} |
        group_combined_tracks
    
    // group_combined_tracks.out | view

    // COMPARISON TRACKS - DIFFREPS,...
    comparison_ch
        .map{it->[it[0],it[1],it[2].getName()]}
        .collectFile{item -> item.join(",")+'\n'} |
        comparison_specific_tracks

    emit:
    individual_tracklines = sample_specific_tracks.out.sid_tracks
    group_combined_tracklines = group_combined_tracks.out
    comparison_specific_tracklines = comparison_specific_tracks.out
}

process comparison_specific_tracks {
    executor 'local'
    //publishDir copy to server
    input:
    path(diffreps_tracks)

    output:
    path("diffreps_tracks.txt")

    script:
    data_path="${workflow.userName}/${params.dataset_label}"
    
    """
    module load R/${params.rversion}
    create_by_comparison_trackline.R \
        --diffreps_tracks ${diffreps_tracks} \
        --data_path ${data_path} \
        --output_name "diffreps_tracks.txt"
    """
}

process combined_track_lines {
    executor 'local'
    
    input:
    tuple path(sid_tracks), path(sample_labels)//, path(files)

    output:
    path("sid_tracks.txt"), emit: sid_tracks
    path("all_bigwig_group_autoscale_hub.txt"), emit: bigwig_hub

    script:
    data_path="${workflow.userName}/${params.dataset_label}"
    
    """
    module load R/${params.rversion}
    generate_sid_tracks.R \
        --sample_labels ${sample_labels} \
        --sid_tracks ${sid_tracks} \
        --data_path ${data_path} \
        --output_name "sid_tracks.txt" \
        --peakcaller ${params.peakcaller}

    ## create all_bigwig_group_autoscale_hub.txt
    convert_tohub.py ${params.dataset_label} ./bigwig_for_hub.csv
    """
}

process group_combined_tracks {
    executor 'local'

    input:
    path(group_tracks)

    output:
    path("bw_combined_tracks.txt")

    script:
    data_path="${workflow.userName}/${params.dataset_label}"

    """
    module load R/${params.rversion}
    generate_bw_combined_tracks.R \
        --combined_tracks ${group_tracks} \
        --data_path ${data_path} \
        --output_name "bw_combined_tracks.txt"
    """
}

process sample_specific_tracks {
    executor 'local'
    
    input:
    tuple path(sid_tracks), path(sample_labels)//, path(files)

    output:
    path("sid_tracks.txt"), emit: sid_tracks
    path("all_bigwig_group_autoscale_hub.txt"), emit: bigwig_hub

    script:
    data_path="${workflow.userName}/${params.dataset_label}"
    
    """
    module load R/${params.rversion}
    generate_sid_tracks.R \
        --sample_labels ${sample_labels} \
        --sid_tracks ${sid_tracks} \
        --data_path ${data_path} \
        --output_name "sid_tracks.txt" \
        --peakcaller ${params.peakcaller}

    ## create all_bigwig_group_autoscale_hub.txt
    convert_tohub.py ${params.dataset_label} ./bigwig_for_hub.csv
    """
}
