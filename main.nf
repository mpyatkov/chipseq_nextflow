#!/usr/bin/env nextflow
nextflow.enable.dsl=2

mm9_chrom_sizes  = file("$projectDir/assets/mm9.chrom.sizes", checkIfExists: true)
mm9_black_complement  = file("$projectDir/assets/mm9-blacklist_complement", checkIfExists: true)
default_tracks  = file("$projectDir/assets/default_tracks.txt", checkIfExists: true)

params.bowtie2_index="/projectnb/wax-es/aramp10/Bowtie2/Bowtie2Index/genome"
params.output_dir="./RESULTS"
params.peakcaller="MACS2"
params.dataset_label="TEST1"
params.copy_to_server_bool=false
params.fastq_config = file("$projectDir/${params.input_configs}/fastq_config.csv", checkIfExists: true)
params.sample_labels_config = file("$projectDir/${params.input_configs}/sample_labels.csv", checkIfExists: true)
params.diffreps_config = file("$projectDir/${params.input_configs}/diffreps_config.csv")

// need_diffexpr = is_empty_file(params.diffreps_config.toString()) ? false : true

// println(params.input_configs)
// println(params.diffreps_config)
// println(params.fastq_config)
// println(params.sample_labels_config)

include {DIFFREPS} from './subworkflow/diffreps/diffreps.nf'
include {MANORM2} from './subworkflow/diffreps/manorm2.nf'

process bowtie2_align {

    tag "${sample_id}"
    
    //echo true
    // executor "local"
    cpus 16
    memory '32 GB'
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "${sample_id}_sorted.bam*", overwrite: true
    publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "*.log", overwrite: true
    publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "library.txt", overwrite: true
    
    input:
    tuple val(sample_id), val(r1), val(r2)
    file(anti_blacklist)
    
    output:
    tuple val(sample_id), val(library), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"),  emit: bam
    tuple val(sample_id), path("*.log"), emit: log
    
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
    echo "${sample_id}: ${library}" > library.txt
    module load bowtie2
    module load samtools
    bowtie2 -p $task.cpus -x $params.bowtie2_index $reads_args -S ${sample_id}_bowtie2.sam 2> ${sample_id}.bowtie2.log
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
    tuple val(sample_id), path("*stats"), emit: stats

    script:
    filter_bedpe=library == "paired-end" ? "and proper_pair" : ""
    bedtools_bedpe=library == "paired-end" ? "-bedpe" : ""
    """
    module load bedtools
    module load samtools
    
    #echo "LIBRARY: ${library}"
    
    sambamba view -h -t $task.cpus -f bam -F "[XS] == null and not unmapped ${filter_bedpe}" $bam > 1.bam
    sambamba sort -n -t $task.cpus 1.bam

    ## extracting fragments
    ## for single-end it will be read themselfs
    ## for paired-end (start1,end1 - start2,end2) it will be difference (start1,end2)
    if [[ ${library} == "single-end" ]]; then
        bedtools bamtobed ${bedtools_bedpe} -i 1.sorted.bam > ${sample_id}_fragments.bed 2> /dev/null
    else
        bedtools bamtobed ${bedtools_bedpe} -i 1.sorted.bam | cut -f 1,2,6 > ${sample_id}_fragments.bed 2> /dev/null
    fi

    awk -v OFS='\t' '{print \$1,\$2,\$3,".","0","."}' ${sample_id}_fragments.bed > ${sample_id}_fragments_bed6.bed

    #gzip ${sample_id}_fragments.bed

    mv 1.sorted.bam 2.bam
    sambamba sort -t $task.cpus 2.bam
    mv 2.sorted.bam ${sample_id}_sorted_filtered.bam    
    sambamba index ${sample_id}_sorted_filtered.bam

    samtools stats ${sample_id}_sorted_filtered.bam > ${sample_id}_sorted_filtered.stats
    samtools flagstats ${sample_id}_sorted_filtered.bam > ${sample_id}_sorted_filtered.flagstats
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
    tuple val(sample_id), path("*narrow_MACS2_peaks.xls"), emit: xls
    tuple val(sample_id), path("*broad_MACS2_peaks.xls"), emit: broad_xls
    tuple val(sample_id), path("*narrowPeak.bed"), emit: narrow_bed
    tuple val(sample_id), path("*broadPeak.bed"), emit: broad_bed
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
    #set -x
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


// process parse_configuration_xls {
//     executor 'local'
//     input:
//     path(xlsx)

//     output:
//     path("sample_labels.csv"), emit: sample_labels_config
//     path("diffreps_config.csv"), emit: diffreps_config, optional: true
//     path("fastq_config.csv"), emit: fastq_config

//     script:
//     """
//     module load R
//     config_parser.R --input_xlsx ${xlsx}
//     """
// }

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


process collect_metrics {
    tag "${sample_id}"
    cpus 1
    memory '16 GB'
    publishDir path: "${params.output_dir}/${sample_id}/metrics/", mode: "copy", overwrite: true

    beforeScript 'source $HOME/.bashrc'
    
    // echo true
    input:
    tuple val(sample_id), val(library), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("*_metrics"), emit: metrics
    tuple val(sample_id), path("*.pdf"), emit: metrics_pdf
    
    script:
    """
    module load picard
    module load R

    # set -x
    
    picard \
        CollectMultipleMetrics \
        --INPUT $bam \
        --OUTPUT ${sample_id}.CollectMultipleMetrics

    ## only for paired-end
    if [[ $library == "paired-end" ]];
    then
    picard \
        CollectInsertSizeMetrics \
        --Histogram_FILE ${sample_id}.InsertSizeMetrics.pdf \
        --INPUT ${bam} \
        --OUTPUT ${sample_id}.InsertSizeMetrics
    fi
    """
}

process fastqc {
    tag "${sample_id}"
    cpus 4
    time '3h'
    memory '32 GB'
    errorStrategy 'retry'
    maxRetries 3
    
    publishDir path: "${params.output_dir}/${sample_id}/metrics/fastqc/", mode: "copy", overwrite: true
    
    beforeScript 'source $HOME/.bashrc'
    
    // echo true
    input:
    tuple val(sample_id), val(r1), val(r2)
    
    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip") , emit: zip

    script:
    if (r2 == "NO") {
        """
        module load fastqc
        mkdir -p ${sample_id}_out
        ln -s $r1 ${sample_id}_1.fq.gz
        fastqc -o ${sample_id}_out --threads $task.cpus ${sample_id}_1.fq.gz
        mv ${sample_id}_out/* ./
        """
    } else {
        """
        module load fastqc
        mkdir -p ${sample_id}_out
        ln -s $r1 ${sample_id}_1.fq.gz
        ln -s $r2 ${sample_id}_2.fq.gz
        fastqc --threads $task.cpus $r1 $r2
        fastqc -o ${sample_id}_out --threads $task.cpus ${sample_id}_1.fq.gz ${sample_id}_2.fq.gz
        mv ${sample_id}_out/* ./
        """
    }
}

process multiqc {

    cpus 1
    publishDir path: "${params.output_dir}/multiqc/", mode: "copy", pattern: "multiqc_report.html", overwrite: true
    publishDir path: "/net/waxman-server/mnt/data/waxmanlabvm_home/${workflow.userName}/${params.dataset_label}/multiqc/", mode: "copy", pattern: "multiqc_report.html", overwrite: true
    
    beforeScript 'source $HOME/.bashrc'
    
    input:
    path(fastqc)          // fastq
    path(aligner)         // bowtie logs
    path(bam_count_stats) // samtools stats
    path(macs2_xls)       // MACS2 xls files
    path(picard)          // picard files
    // echo true
    
    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    

    script:
    """
    module load multiqc
    multiqc -f .
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

process checkifempty {
    executor "local"
    echo true
    input:
    path(f)
    output:
    stdout

    script:
    println "In process: " + f.getClass()
    """
    echo "FILE: $f"
    """
    
}

def is_empty_file(fp) {
    File file = new File(fp);
    return !file.exists() || file.length() == 0
}

workflow {
    // parse_configuration_xls(params.xlsx_config)
    // parse_configuration_xls.out.sample_labels_config | view
    // parse_configuration_xls.out.diffreps_config | view
    // test = parse_configuration_xls.out.diffreps_config.ifEmpty(false)
    // .ifEmpty{exit 1, "Cannot find any input RDS files for Module 3"}

    fastq_config_ch = Channel.from(params.fastq_config)
    sample_labels_config_ch = Channel.from(params.sample_labels_config)
    diffreps_config_ch = Channel.from(params.diffreps_config)
    
    // if (is_empty_file(params.diffreps_config.toString())) {
    //     println(params.diffreps_config.toString())
    //     println("File is empty")
    // } else {
    //     println("File is not empty")
    // }

    // parse_configuration_xls.out.diffreps_config.splitCsv() 
    // parse_configuration_xls.out.fastq_config.splitCsv() | read12_tester
    // checkifempty(parse_configuration_xls.out.diffreps_config)
    // parse_configuration_xls.out.diffreps_config | view
    
    // // fastq = Channel.fromPath( './hnf6.csv' ).splitCsv()
    // fastq_reads_ch = parse_configuration_xls.out.fastq_config.splitCsv()
    // fastq_reads_ch = params.fastq_config.splitCsv() | view()

    bowtie2_align(fastq_config_ch.splitCsv(), mm9_black_complement)
    
    bam_count(bowtie2_align.output.bam)
    
    macs2_callpeak(bam_count.output.final_bam, mm9_chrom_sizes)
    epic2_callpeak(bam_count.output.fragments_bed6, mm9_chrom_sizes)
    
    
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
        extra_columns = macs2_callpeak.out.xls // for diffreps summary
        peaks_for_manorm2 = macs2_callpeak.out.narrow_bed
    } else {
        extra_columns = epic2_callpeak.out.bed6 // for diffreps summary
        peaks_for_manorm2 = epic2_callpeak.out.bed6
    }

    DIFFREPS(
        // parse_configuration_xls.out.diffreps_config,
        diffreps_config_ch,
        // parse_configuration_xls.out.sample_labels_config,
        sample_labels_config_ch,
        bam_count.out.fragments_bed6,
        calc_norm_factors.out,
        extra_columns, //macs2_callpeak.out.xls
        mm9_chrom_sizes
    )

    MANORM2(
        diffreps_config_ch,
        bam_count.out.fragments,
        peaks_for_manorm2,
        mm9_chrom_sizes
    )

    // MANORM2.out.profile | view
    manorm2_group_report_ch = MANORM2.out.diff_table
        .map{meta, rest -> [meta.group_name, rest]}
    
    diffreps_group_report_ch = DIFFREPS.out.full_report
        .map{meta, rest -> [meta.group_name, rest]}
        .groupTuple() 

    combined_manorm2_diffreps_ch = diffreps_group_report_ch.join(manorm2_group_report_ch)
    diffreps_manorm2_overlap(combined_manorm2_diffreps_ch)
    
    // TRACKS (copying, generating track lines)
    // create_bigwig_files 
    sid_normfact_ch = calc_norm_factors.out.splitCsv(sep: "\t")
        .map{it->[it[0],it[4]]}
    
    sid_fr_norm = bam_count.out.fragments
        .join(sid_normfact_ch)
    
    create_bigwig_files(sid_fr_norm, mm9_chrom_sizes)


    // create track lines (sample specific)
    sid_specific_ch = create_bigwig_files.out
        .mix(macs2_callpeak.out.broad_bb)
        .mix(macs2_callpeak.out.narrow_bb)
        .mix(epic2_callpeak.out.epic2_bb)
        //.mix(bam_count.out.final_bam.map{it->[it[0],it[2]]}) //excluded bam track lines
        .map{it->[it[0], it[1].getName()]}
        .collectFile{item -> item.join(",")+'\n'}
        .combine(sample_labels_config_ch)

    // create track lines (group specific)    
    group_specific_ch = DIFFREPS.out.diffreps_track
        .mix(MANORM2.out.manorm2_track)
        .map{it->[it[0], it[2].getName()]}
        .collectFile{item -> item.join(",")+'\n'}
        .combine(diffreps_config_ch) | view

    create_sample_specific_tracks(sid_specific_ch)

    create_diffreps_tracks(group_specific_ch)

    track_lines = Channel.from(default_tracks).concat(create_diffreps_tracks.out,create_sample_specific_tracks.out)
        .collectFile(name: 'autolimit_tracks.txt', sort: 'index'){item -> item.text}

    // copy bb,bam,bw files to server
    bw_files = create_bigwig_files.out.map{it -> it[1]} 
    narrow_files= macs2_callpeak.out.narrow_bb.map{it -> it[1]}
    broad_files= macs2_callpeak.out.broad_bb.map{it -> it[1]}
    broad_epic_files=epic2_callpeak.out.epic2_bb.map{it -> it[1]}
    bam_files=bam_count.out.final_bam.map{it->[it[2],it[3]]}
    diffreps_files=DIFFREPS.out.diffreps_track.map{it->it[2]}
    manorm2_files=MANORM2.out.manorm2_track.map{it->it[2]}
    
    // track_files_to_server = bw_files.concat(bam_files, broad_epic_files, broad_files, narrow_files, diffreps_files).collect() //excluded bam track files
    track_files_to_server = bw_files.concat(broad_epic_files, broad_files, narrow_files, diffreps_files, manorm2_files).collect()
    
    if (params.copy_to_server_bool){
        copy_files_to_server(track_files_to_server, track_lines)
    }

    // calculate overlaps for all MACS2 narrow/broad peaks
    macs2_narrow_broad_ch = macs2_callpeak.out.xls
        .mix(macs2_callpeak.out.broad_xls)
        .map{sid,bedfile -> bedfile}
        .collect()

    extradetailed_stats_for_macs2(macs2_narrow_broad_ch, params.sample_labels_config)
    
    //picard
    collect_metrics(bam_count.out.final_bam)
    
    fastqc(fastq_config_ch.splitCsv())

    // bowtie2_align.out.log | view
    // bam_count.out.stats | view
    multiqc(
        fastqc.out.zip.map{it -> it[1]}.collect(),
        bowtie2_align.out.log.map{it -> it[1]}.collect(),
        bam_count.out.stats.map{it -> it[1]}.collect(),
        macs2_callpeak.out.xls.map{it -> it[1]}.collect(),
        collect_metrics.out.metrics.map{it -> it[1]}.collect()
    )

    multiqc.out.report | view
    // read12_tester

    // log.info """\
    //      --
    // run as       : ${workflow.commandLine}
    // user       : ${workflow.userName}
         
    //      config files : ${workflow.configFiles}
         
    //      """
    //      .stripIndent()

}

//TODO: To avoid multiple if conditions that diffreps config exists
// it is required to pack creating tracks to separate workflow
// if (diffreps) {
//     foo(input1)
// } else {
//     foo(input2)
// }


process diffreps_manorm2_overlap {
    tag "${group_name}"
    executor 'local'
    publishDir path: "${params.output_dir}/manorm2_diffreps_overlap/", mode: "copy", pattern: "*.xlsx", overwrite: true
    
    beforeScript 'source $HOME/.bashrc'
    input:
    tuple val(group_name), path(diffreps_reports), path(manorm2_report)
    
    output:
    tuple val(group_name), path("*.xlsx")

    script:
    """
    module load R
    diffreps_manorm2_overlap.R --output_prefix ${group_name}
    """
}

process extradetailed_stats_for_macs2 {
    
    executor "local"
    // echo true
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/MACS2_peaks_overlaps/", mode: "copy", pattern: "*.xlsx", overwrite: false

    input:
    path(peaks)
    path(sample_labels)
    
    output:
    path("*.xlsx")

    script:
    """
    module load R
    detailed_macs2_peaks_overlap.R --peak_prefix "broad" --sample_labels ./sample_labels.csv
    detailed_macs2_peaks_overlap.R --peak_prefix "narrow" --sample_labels ./sample_labels.csv
    """
}
