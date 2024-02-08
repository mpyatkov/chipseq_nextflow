#!/usr/bin/env nextflow
nextflow.enable.dsl=2

mm9_chrom_sizes  = file("$projectDir/assets/mm9.chrom.sizes", checkIfExists: true)
mm9_black_complement  = file("$projectDir/assets/mm9-blacklist_complement", checkIfExists: true)
default_tracks  = file("$projectDir/assets/default_tracks.txt", checkIfExists: true)


params.bowtie2_index="/projectnb/wax-es/aramp10/Bowtie2/Bowtie2Index/genome"
params.output_dir="./RESULTS"
params.peakcaller="MACS2"
params.dataset_label="TEST1"

// params.xlsx_config = file("$projectDir/hnf6_samples.xlsx", checkIfExists: true)
// params.xlsx_config = file("./hnf6_samples.xlsx", checkIfExists: true)
params.xlsx_config = file(params.input_config, checkIfExists: true)

include {DIFFREPS} from './subworkflow/diffreps/diffreps.nf'

process bowtie2_align {

    tag "${sample_id}"
    
    //echo true
    // executor "local"
    cpus 16
    memory '32 GB'
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "${sample_id}_sorted.bam*", overwrite: true
    // publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "library.txt", overwrite: true
    
    input:
    tuple val(sample_id), val(r1), val(r2)
    file(anti_blacklist)
    
    output:
    tuple val(sample_id), val(library), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"),  emit: bam
    
    script:

    reads_args = null
    library = null
    
    if (r2 == "NO") {
        reads_args = "-U ${r1}"
        library = "single-end"
    } else {
        reads_args = "-1 ${r1} -2 ${r2}"
        library = "paired-end"
    }

    // TODO: remove prefixes in script and use usual sam and bam, rename later in publishdir
    
    """
    #echo "${sample_id} library: ${library}"
    module load bowtie2
    module load samtools
    bowtie2 -p $task.cpus -x $params.bowtie2_index $reads_args -S ${sample_id}_bowtie2.sam
    samtools view -bS ${sample_id}_bowtie2.sam > ${sample_id}_alignments.bam
    samtools sort -@ $task.cpus ${sample_id}_alignments.bam -o "sorted.bam"
    samtools index -@ $task.cpus "sorted.bam"

    rm -rf ${sample_id}_alignments.bam ${sample_id}_bowtie2.sam
    
    ## removing blacklist regions from bam
    samtools view -@ $task.cpus -L ${anti_blacklist} -O BAM -o "${sample_id}_sorted.bam" sorted.bam
    samtools index -@ $task.cpus "${sample_id}_sorted.bam"
    """

    stub:
    """
    touch ${sample_id}_sorted.bam
    touch ${sample_id}_sorted.bam.bai
    touch ${sample_id}_library.report
    """
}

process bam_count {
    tag "${sample_id}"
    
    // echo true
    cpus 4
    memory '32 GB'
    // executor 'local'

    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "*fragments*.bed*", overwrite: true
    publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "${sample_id}_sorted_filtered.bam*", overwrite: true
    
    input:
    tuple val(sample_id), val(library), path(bam), path(bai)
    
    output:
    
    tuple val(sample_id), path("${sample_id}_fragments.bed"), emit: fragments
    tuple val(sample_id), path("${sample_id}_fragments_bed6.bed"), emit: fragments_bed6
    tuple val(sample_id), val(library), path("${sample_id}_sorted_filtered.bam"), path("${sample_id}_sorted_filtered.bam.bai"), emit: final_bam

    script:
    filter_bedpe=library == "paired-end" ? "and proper_pair" : ""
    bedtools_bedpe=library == "paired-end" ? "-bedpe" : ""
    """
    module load bedtools
    #echo "LIBRARY: ${library}"
    
    sambamba view -h -t $task.cpus -f bam -F "[XS] == null and not unmapped ${filter_bedpe}" $bam > 1.bam
    sambamba sort -n -t $task.cpus 1.bam
    bedtools bamtobed ${bedtools_bedpe} -i 1.sorted.bam > ${sample_id}_fragments.bed 2> /dev/null

    awk -v OFS='\t' '{print \$1,\$2,\$3,".","0","."}' ${sample_id}_fragments.bed > ${sample_id}_fragments_bed6.bed

    #gzip ${sample_id}_fragments.bed

    mv 1.sorted.bam 2.bam
    sambamba sort -t $task.cpus 2.bam
    mv 2.sorted.bam ${sample_id}_sorted_filtered.bam    
    sambamba index ${sample_id}_sorted_filtered.bam
    """
    
    stub:
    x = "paired"
    """
    touch ${sample_id}_fragments.bed.gz
    touch ${sample_id}_sorted_mapped.bam
    """
    
}

process macs2_callpeak {

    tag "${sample_id}"
    
    // echo true
    // executor "local"
    cpus 1
    memory '8 GB'
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/${sample_id}/macs2/", mode: "copy", pattern: "*.{narrowPeak,broadPeak,xls,bed}", overwrite: true
        
    executor 'sge'
    beforeScript 'source $HOME/.bashrc'

    input:
    tuple val(sample_id), val(library), path(bam), path(bai)
    path(mm9_chrom_sizes)
    
    output:
    tuple val(sample_id), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(sample_id), path("*narrow_MACS2_peaks.xls")                   , emit: xls
    tuple val(sample_id), path("*narrowPeak.bed")       , emit: narrow_bed
    tuple val(sample_id), path("*broadPeak.bed")       , emit: broad_bed
    tuple val(sample_id), path("*narrow*.bb"), emit: narrow_bb
    tuple val(sample_id), path("*broad*.bb"), emit: broad_bb


    shell:
    lib=library == "paired-end" ? "BAMPE" : "BAM"
    template 'macs2_callpeak.sh'

    stub:
    """
    touch ${sample_id}.narrowPeak
    touch ${sample_id}.broadPeak
    touch ${sample_id}.narrow.xls
    touch ${sample_id}.broad.xls
    touch ${sample_id}_narrow.bb
    touch ${sample_id}_broad.bb
    echo "${sample_id}_narrow" > ${sample_id}.narrowPeak.bed
    echo "${sample_id}_broad"  > ${sample_id}.broadPeak.bed
    """
}

process peak_union {

    executor "local"
    echo true
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/peak_union/", mode: "copy", pattern: "peak_union.bed", overwrite: true

    input:
    path("*")
    
    output:
    path("peak_union.bed"), emit: union

    script:
    """
    module load bedtools
    set -x
    cat * > tmp.bed
    sort -k1,1 -k2,2n tmp.bed > tmp1.bed
    bedtools merge -i tmp1.bed > peak_union.bed
    """
    
    stub:
    """
    cat * > peak_union.bed
    """
}

process fragments_union_overlap {
    tag "${sample_id}"
    
    executor "sge"
    cpus 1
    memory '8 GB'
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/${sample_id}/macs2/", mode: "copy", pattern: "${sample_id}_fragments_union_coverage.bed", overwrite: true
        
    input:
    tuple val(sample_id), path(fragments)
    path(union)
    
    output:
    tuple val(sample_id), path("${sample_id}_fragments_union_coverage.bed"), emit: fr_union_overlap

    script:
    """
    module load bedtools
    bedtools coverage -a $union -b $fragments > "${sample_id}_fragments_union_coverage.bed"
    """

}

process calc_sample_stats {
    // executor "sge"
    // cpus 1
    executor 'local'
    echo true
    input:
    tuple val(sample_id), path(fragments), path(fragments_union_coverage), path(union)
    
    output:
    path("${sample_id}_output.txt")
    script:

    // ## FRAGMENT_COUNT=\$(grep 'in treatment:' $xls | awk '{print \$NF}') ## macs2
    // ## FRAGMENT_COUNT=\$(grep -A 1 'Line count of fragments BED file:' *\${Sample_ID}'.o'* | awk 'FNR==2{print \$0}') ## sicer
    // ## echo "${sample_id} - fr count: \${FRAGMENT_COUNT} - \${FRAGMENT_IN_PEAK_COUNT}"    
    // ## echo "${sample_id},\${FRAGMENT_COUNT},\${FRAGMENT_IN_PEAK_COUNT},\${FRAGMENT_IN_PEAK_RATIO}" 

    """
    FRAGMENT_COUNT=\$(cat $fragments | wc -l )
    FRAGMENT_IN_PEAK_COUNT=\$(awk '{n+=\$4;} ; END {print n;}' ${fragments_union_coverage})
    FRAGMENT_IN_PEAK_RATIO=\$(echo "scale=4;\${FRAGMENT_IN_PEAK_COUNT}/\${FRAGMENT_COUNT}" | bc)
    echo "${sample_id},\${FRAGMENT_COUNT},\${FRAGMENT_IN_PEAK_COUNT},\${FRAGMENT_IN_PEAK_RATIO}" > ${sample_id}_output.txt
    """
    
    stub:
    """
    fragment_count=${sample_id}
    fragment_in_peak_count=${sample_id}
    fragment_in_peak_ratio=${sample_id}
    echo "${sample_id},\${fragment_count},\${fragment_in_peak_count},\${fragment_in_peak_ratio}" > ${sample_id}_output.txt
    """
}

process calc_norm_factors {

    executor "local"
    beforeScript 'source $HOME/.bashrc'
    publishDir path: "${params.output_dir}/peak_union/", mode: "copy", pattern: "*.tsv", overwrite: true
        
    input:
    path(sample_stats)
    
    output:
    path("Norm_Factors.tsv")
    
    script:
    """
    module load R
    calculate_norm_factors.R ${sample_stats} "Norm_Factors.tsv"
    """
}

process create_bigwig_files {

    tag "${sample_id}"
    executor 'sge'
    cpus 1
    time '1h'
    memory '16 GB'

    beforeScript 'source $HOME/.bashrc'
    
    input:
    tuple val(sample_id), path(fragments), val(norm_factor)
    path(mm9_chrom_sizes)
    
    output:
    tuple val(sample_id), path("*.bw")
    
    script:
    """
    module load bedtools
    module load ucscutils

    bedtools sort -i ${fragments} > tmp.sorted.bed
    bedtools genomecov -bg -i tmp.sorted.bed -g ${mm9_chrom_sizes} > tmp.bedGraph
    sed -i '/track/d;/random/d;/chrM/d' tmp.bedGraph
    awk -v OFS='\t' -v norm_factor="${norm_factor}" '{print \$1, \$2, \$3, \$4 * norm_factor}' "tmp.bedGraph" > "tmp.norm.bedGraph"
    bedGraphToBigWig "tmp.norm.bedGraph" ${mm9_chrom_sizes} "${sample_id}_RiPPM_norm.bw"
    """
}

process create_sample_specific_tracks {
    executor 'local'
    
    input:
    tuple path(sid_tracks), path(sample_labels)//, path(files)

    output:
    path("sid_tracks.txt")

    script:
    data_path="${workflow.userName}/${params.dataset_label}"
    
    """
    module load R
    generate_sid_tracks.R \
        --sample_labels ${sample_labels} \
        --sid_tracks ${sid_tracks} \
        --data_path ${data_path} \
        --output_name "sid_tracks.txt"
    """
}

process create_diffreps_tracks {
    executor 'local'
    //publishDir copy to server
    input:
    tuple path(diffreps_tracks), path(diffreps_config)

    output:
    path("diffreps_tracks.txt")

    script:
    data_path="${workflow.userName}/${params.dataset_label}"
    
    """
    module load R
    generate_diffreps_tracks.R \
        --diffreps_config ${diffreps_config} \
        --diffreps_tracks ${diffreps_tracks} \
        --data_path ${data_path} \
        --output_name "diffreps_tracks.txt"
    """
}


process parse_configuration_xls {
    executor 'local'
    input:
    path(xlsx)

    output:
    path("sample_labels.csv"), emit: sample_labels_config
    path("diffreps_config.csv"), emit: diffreps_config
    path("fastq_config.csv"), emit: fastq_config

    script:
    """
    module load R
    config_parser.R --input_xlsx ${xlsx}
    """
}

process epic2_callpeak {
    tag "${sample_id}"
    
    // echo true
    // executor "local"
    
    cpus 1
    memory '8 GB'
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/${sample_id}/epic2/", mode: "copy", pattern: "*bed", overwrite: true
        
    executor 'sge'
    beforeScript 'source $HOME/.bashrc'

    input:
    tuple val(sample_id), path(fragments)
    path(mm9_chrom_sizes)
    
    output:
    tuple val(sample_id), path("${sample_id}_epic_bed6.bed"), emit: bed6
    tuple val(sample_id), path("${sample_id}_epic_bed3.bed"), emit: bed3
    tuple val(sample_id), path("*epic2*.bb"), emit: epic2_bb

    script:
    //species=mm9
    window_size=400
    fragment_size=200
    effective_genome_fraction=0.80
    gap_size=2400 // 6
    e_value=100
    
    """
    module load epic2
    module load bedtools
    epic2 --output "${sample_id}_epic_bed6.bed" -egf ${effective_genome_fraction} -cs ${mm9_chrom_sizes} -t ${fragments} -e ${e_value} -fs ${fragment_size} -bin ${window_size} -g 6
    awk 'OFS="\t" {print \$1,\$2,\$3}' "${sample_id}_epic_bed6.bed" > "${sample_id}_epic_bed3.bed"
    
    cat "${sample_id}_epic_bed3.bed" | grep -vE "track|chrM|random" > tmp.bed
    bedtools sort -i tmp.bed > tmp.sorted.bed
    bedToBigBed -allow1bpOverlap tmp.sorted.bed ${mm9_chrom_sizes} "${sample_id}_epic2.bb"
    """
}

process read12_tester {
    executor 'local'
    echo true
    
    input:
    tuple val(sample_id), val(r1), val(r2)

    output:
    stdout

    script:
    library = null
    if (r2 == "NO") {
        library = "single-end"
    } else {
        library = "paired-end"
    }
    // println(Objects.equals(r2, new String("NA")))
    // println(r2.class)
    // println(na.class)
    // println(r2)
    // println(na)
    """
    #echo $r2
    echo "library: $library"
    """
}

process copy_files_to_server {
    executor 'local'
    // publishDir path: "${data_path}", mode: "copy", pattern: "*.{bw,bam*,bb}", overwrite: false
    // publishDir path: "${data_path}/TRACK_LINES", mode: "copy", pattern: "*.txt", overwrite: true
    
    input:
    path(files)
    path(track_lines)

    // output:
    // path("*.{bw,bam,bai,bb,txt}")
    // stdout

    script:
    data_path="/net/waxman-server/mnt/data/waxmanlabvm_home/${workflow.userName}/${params.dataset_label}"
    """
    mkdir -p ${data_path}/TRACK_LINES
    cp ${files} ${data_path}
    cp ${track_lines} ${data_path}/TRACK_LINES/
    """
}

// process test1 {
//     executor 'local'
//     echo true

//     input:
//     val(arr)

//     output:
//     // path("result.txt")
//     stdout
//     script:
//     def reports = arr.join(" ")
    
//     """
//     echo $reports
//     for x in ${reports}; do
//         echo \$x
//     done
//     """
// }



workflow {

    parse_configuration_xls(params.xlsx_config)
    // parse_configuration_xls.out.sample_labels_config | view
    // parse_configuration_xls.out.diffreps_config | view
    // parse_configuration_xls.out.fastq_config.splitCsv() | view
    // parse_configuration_xls.out.fastq_config.splitCsv() | read12_tester

    
    // fastq = Channel.fromPath( './hnf6.csv' ).splitCsv()
    bowtie2_align(parse_configuration_xls.out.fastq_config.splitCsv(),
                  mm9_black_complement)
    
    bam_count(bowtie2_align.output.bam)
    
    macs2_callpeak(bam_count.output.final_bam, mm9_chrom_sizes)
    epic2_callpeak(bam_count.output.fragments, mm9_chrom_sizes)
    
    
    //peak union (MACS2 / EPIC2(SICER))
    if (params.peakcaller == "MACS2") {
        peakcaller_bed_files = macs2_callpeak.out.narrow_bed
            .map{it->it[1]}.collect()
    } else {
        peakcaller_bed_files = epic2_callpeak.out.bed3
            .map{it->it[1]}.collect()
    }

    peak_union(peakcaller_bed_files)


    //overlap peakcaller union and fragments files
    fragments_union_overlap(bam_count.out.fragments, peak_union.out.union)

    
    peaks_for_stats_ch = bam_count.out.fragments
        .join(fragments_union_overlap.out.fr_union_overlap)
        .combine(peak_union.out.union)
    
    calc_sample_stats(peaks_for_stats_ch)
    
    sample_stats = calc_sample_stats.out
        .collectFile{item -> item.text}

    sample_stats | calc_norm_factors

    // log.info("params.peakcaller: $params.peakcaller")

    if (params.peakcaller == "MACS2") {
        extra_columns = macs2_callpeak.out.xls
    } else {
        extra_columns = epic2_callpeak.out.bed6
    }

    calc_norm_factors.out

    DIFFREPS(
        parse_configuration_xls.out.diffreps_config,
        parse_configuration_xls.out.sample_labels_config,
        bam_count.out.fragments_bed6,
        calc_norm_factors.out,
        extra_columns, //macs2_callpeak.out.xls
        mm9_chrom_sizes
    )
    
    // TRACKS (copying, generating track lines)
    // create_bigwig_files 
    sid_normfact_ch = calc_norm_factors.out.splitCsv(sep: "\t")
        .map{it->[it[0],it[4]]}
    
    sid_fr_norm = bam_count.out.fragments
        .join(sid_normfact_ch)
    
    create_bigwig_files(sid_fr_norm, mm9_chrom_sizes)


    // create track lines
    sid_specific_ch = create_bigwig_files.out
        .mix(macs2_callpeak.out.broad_bb)
        .mix(macs2_callpeak.out.narrow_bb)
        .mix(epic2_callpeak.out.epic2_bb)
        .mix(bam_count.out.final_bam.map{it->[it[0],it[2]]})
        .map{it->[it[0], it[1].getName()]}
        .collectFile{item -> item.join(",")+'\n'}
        .combine(parse_configuration_xls.out.sample_labels_config)


    group_specific_ch = DIFFREPS.out.diffreps_track 
        .map{it->[it[0], it[2].getName()]}
        .collectFile{item -> item.join(",")+'\n'}
        .combine(parse_configuration_xls.out.diffreps_config)


    create_sample_specific_tracks(sid_specific_ch)
    // create_sample_specific_tracks.out |view

    create_diffreps_tracks(group_specific_ch)

    track_lines = Channel.from(default_tracks).concat(create_diffreps_tracks.out,create_sample_specific_tracks.out)
        .collectFile(name: 'autolimit_tracks.txt', sort: 'index'){item -> item.text} | view

    // copy bb,bam,bw files to server
    bw_files = create_bigwig_files.out.map{it -> it[1]} 
    narrow_files= macs2_callpeak.out.narrow_bb.map{it -> it[1]}
    broad_files= macs2_callpeak.out.broad_bb.map{it -> it[1]}
    broad_epic_files=epic2_callpeak.out.epic2_bb.map{it -> it[1]}
    bam_files=bam_count.out.final_bam.map{it->[it[2],it[3]]}
    diffreps_files=DIFFREPS.out.diffreps_track.map{it->it[2]}
    
    track_files_to_server = bw_files.concat(bam_files, broad_epic_files, broad_files, narrow_files, diffreps_files).collect()
    copy_files_to_server(track_files_to_server, track_lines)

    
    // read12_tester

    // log.info """\
    //      --
    // run as       : ${workflow.commandLine}
    // user       : ${workflow.userName}
         
    //      config files : ${workflow.configFiles}
         
    //      """
    //      .stripIndent()

}

