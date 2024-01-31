#!/usr/bin/env nextflow
nextflow.enable.dsl=2

mm9_chrom_sizes  = file("$projectDir/assets/mm9.chrom.sizes", checkIfExists: true)
mm9_black_complement  = file("$projectDir/assets/mm9-blacklist_complement", checkIfExists: true)


params.bowtie2_index="/projectnb/wax-es/aramp10/Bowtie2/Bowtie2Index/genome"
params.output_dir="./RESULTS"
params.peakcaller="MACS2"
params.dataset_label="TEST1"

// params.xlsx_config = file("$projectDir/hnf6_samples.xlsx", checkIfExists: true)
// params.xlsx_config = file("./hnf6_samples.xlsx", checkIfExists: true)
params.xlsx_config = file(params.input_config, checkIfExists: true)

process bowtie2_align {

    tag "${sample_id}"
    
    echo true
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
    echo "${sample_id} library: ${library}"
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
    
    echo true
    cpus 4
    memory '32 GB'
    // executor 'local'

    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/${sample_id}/bam/", mode: "copy", pattern: "*.gz*", overwrite: true
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
    echo "LIBRARY: ${library}"
    
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

process macs2_union {

    executor "local"
    echo true
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/macs2_union/", mode: "copy", pattern: "macs2_union.bed", overwrite: true

    input:
    path("*")
    
    output:
    path("macs2_union.bed"), emit: union

    script:
    """
    module load bedtools
    set -x
    cat * > tmp.bed
    sort -k1,1 -k2,2n tmp.bed > tmp1.bed
    bedtools merge -i tmp1.bed > macs2_union.bed
    """
    
    stub:
    """
    cat * > macs2_union.bed
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
    tuple val(sample_id), path(xls), path(fragments_union_coverage), path(union)
    
    output:
    path("${sample_id}_output.txt")
    script:

    """
    FRAGMENT_COUNT=\$(grep 'in treatment:' $xls | awk '{print \$NF}') ## macs2
    ## FRAGMENT_COUNT=\$(grep -A 1 'Line count of fragments BED file:' *\${Sample_ID}'.o'* | awk 'FNR==2{print \$0}') ## sicer
    FRAGMENT_IN_PEAK_COUNT=\$(awk '{n+=\$4;} ; END {print n;}' ${fragments_union_coverage})
    ## echo "${sample_id} - fr count: \${FRAGMENT_COUNT} - \${FRAGMENT_IN_PEAK_COUNT}"
    FRAGMENT_IN_PEAK_RATIO=\$(echo "scale=4;\${FRAGMENT_IN_PEAK_COUNT}/\${FRAGMENT_COUNT}" | bc)
    ## echo "${sample_id},\${FRAGMENT_COUNT},\${FRAGMENT_IN_PEAK_COUNT},\${FRAGMENT_IN_PEAK_RATIO}" 
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
    publishDir path: "${params.output_dir}/macs2_union/", mode: "copy", pattern: "*.tsv", overwrite: true
        
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

process diffreps {
    tag "${NUM}"
    executor 'sge'
    cpus 16
    //cache false
    time '2h'
    // echo true

    beforeScript 'source $HOME/.bashrc'
    
    input:
    tuple val(NUM), val(TREATMENT_NAME), val(CONTROL_NAME), val(TREATMENT_SAMPLES), val(CONTROL_SAMPLES),
        val(NORMALIZATION), val(WINDOW_SIZE), path(TREATMENT_FILES), path(CONTROL_FILES), path(NORM_FILE)
    output:
    tuple val(NUM), val(report_name), path("diffReps_*")

    shell:
    report_name="diffReps_${TREATMENT_NAME}.vs.${CONTROL_NAME}_${NORMALIZATION}_${WINDOW_SIZE}"
    
    template 'diffreps.sh'
}

process diffreps_summary {
    tag "${NUM}"

    executor 'local'
    // cpus 16
    // echo true
    
    publishDir path: "${params.output_dir}/diffreps_output/", mode: "copy", pattern: "${output_dir}/*", overwrite: true
    
    input:
    tuple val(NUM), val(TREATMENT_NAME), val(CONTROL_NAME), val(TREATMENT_SAMPLES), val(CONTROL_SAMPLES),
        val(NORMALIZATION), val(WINDOW_SIZE), val(report_name), path(diffout), path(xls), path(SAMPLE_LABELS)
    path(mm9_chrom_sizes)

    
    output:
    tuple val(NUM), path("${output_dir}/*")
    tuple val(report_name), path("Hist*.pdf"), emit: hist_pdf
    tuple val(report_name), path("FDR*.pdf"), emit: fdr_pdf
    tuple val(report_name), path("Bar*.pdf"), emit: bar_pdf
    tuple val(NUM), val(output_dir), path("*.bb"), emit: diffreps_track
    path("${report_name}"), emit: diffreps_output    

    shell:
    peakcaller="${params.peakcaller}"
    output_dir="${NUM}_${report_name}_${peakcaller}"
    
    dataset_label="${params.dataset_label}"
    template 'diffreps_summary.sh'
}

process aggregate_diffreps_pdf {
    tag("${group_name}")
    executor 'local'
    publishDir path: "${params.output_dir}/diffreps_output/aggregated_pdfs/${group_name}", mode: "copy", pattern: "${group_name}_*.pdf", overwrite: true
    
    input:
    tuple val(group_name), path(hist), path(fdr), path(bar)
    output:
    path("${group_name}_*.pdf")
    script:
    
    """
    module load poppler
    pdfunite $hist "${group_name}_Histogram_Barcharts_${params.peakcaller}.pdf"
    pdfunite $fdr "${group_name}_FDR_Barcharts_${params.peakcaller}.pdf"
    pdfunite $bar "${group_name}_Barcharts_${params.peakcaller}.pdf"
    """
}

process collect_diffreps_norm_factors {

    executor 'local'
    publishDir path: "${params.output_dir}/diffreps_output/", mode: "copy", pattern: "*.xlsx", overwrite: true
    echo true

    input:
    tuple path(diffreps_output), path(sample_stats)
    output:
    path("*.xlsx")

    script:
    """
    module load R
    diffreps_output_parser.R --path "." \
        --rippm_report ${sample_stats}
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
    //publishDir copy to server
    input:
    tuple path(sid_tracks), path(sample_labels)

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

    all_macs2 = macs2_callpeak.out.narrow_bed
        .map{it->it[1]}.collect()
    
    macs2_union(all_macs2)

    fragments_union_overlap(bam_count.output.fragments, macs2_union.output.union)
    // fragments_union_overlap.out.fr_union_overlap 

    peaks_for_stats = macs2_callpeak.out.xls
        .join(fragments_union_overlap.out.fr_union_overlap)
        .combine(macs2_union.out.union)
    
    calc_sample_stats(peaks_for_stats)
    sample_stats = calc_sample_stats.out
        .collectFile{item -> item.text}

    sample_stats | calc_norm_factors


    //diffreps_config = Channel.fromPath( './diffreps_config.csv' ).splitCsv()
    // diffreps_params = diffreps_config
    diffreps_params = parse_configuration_xls.out.diffreps_config.splitCsv()
        .combine(
            //bed3_to_bed6.out.bed6
            bam_count.out.fragments_bed6
                .map{it->it[1]}
                .collect().toList())
        .map{num,tn,cn,tsmp,csmp,norm,w,l->
            def tl = l.findAll{it =~ tsmp.trim()} 
            def cl = l.findAll{it =~ csmp.trim()}
            return [num.trim(),tn.trim(),cn.trim(),tsmp.trim(),csmp.trim(),norm.trim(),w.trim(),tl,cl]}
        
    // diffreps(diffreps_params) | view

    diffreps_norm_params = diffreps_params
        .map{it->[it[0], it[3].trim(), it[4].trim()]}

    // diffreps_norm_params.combine(calc_norm_factors.out) | diffreps_rippm_norm
    // diffreps_rippm_norm.out

    diffreps_params_withnorm = diffreps_params.combine(calc_norm_factors.out)
    diffreps_params_withnorm | diffreps

    // diffreps summary
    //sample_labels = Channel.from("$projectDir/Sample_Labels.txt")

    // diffreps.out | view
    
    dsummary_ch = diffreps_params.map{it -> it[0..-3]} // exclude bed6 files
        .join(diffreps.out)
        .combine(macs2_callpeak.out.xls
                 .map{it -> it[1]}.collect().toList())
        .map{it ->
            //filter xls only participating in comparison
            def new_xls = it[-1].findAll{it =~ "${it[3]}|${it[4]}"} 
            def new_out = it[0..-2] << new_xls
            return new_out}
    // .combine(sample_labels)
        .combine(parse_configuration_xls.out.sample_labels_config)
    
    diffreps_summary(dsummary_ch, mm9_chrom_sizes)
    
    // diffreps_summary.out.diffreps_output | view

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
    
    // -- collect diffreps norm factors
    diffreps_summary.out.diffreps_output
        .collect()
        .toList()
        .combine(calc_norm_factors.out) | collect_diffreps_norm_factors

    // collect_diffreps_norm_factors.out | view

    //creating tracks
    // create_bigwig_files 
    sid_normfact_ch = calc_norm_factors.out.splitCsv(sep: "\t")
        .map{it->[it[0],it[4]]}
    
    sid_fr_norm = bam_count.out.fragments
        .join(sid_normfact_ch)
    
    create_bigwig_files(sid_fr_norm, mm9_chrom_sizes)

    // create_bigwig_files.out | view
    // macs2_callpeak.out.broad_bb | view
    // macs2_callpeak.out.narrow_bb | view
    // diffreps_summary.out.diffreps_track | view

    bw_files = create_bigwig_files.out.map{it -> it[1]} 
    narrow_files= macs2_callpeak.out.narrow_bb.map{it -> it[1]}
    broad_files= macs2_callpeak.out.broad_bb.map{it -> it[1]}
    diffreps_files=diffreps_summary.out.diffreps_track.map{it->it[2]}
    // bam_files=bam_count.out.final_bam.map{it->[it[2],it[3]]} | view

    sid_specific_ch = create_bigwig_files.out
        .mix(macs2_callpeak.out.broad_bb)
        .mix(macs2_callpeak.out.narrow_bb)
        .map{it->[it[0], it[1].getName()]}
        .collectFile{item -> item.join(",")+'\n'}
    // .combine(sample_labels)
        .combine(parse_configuration_xls.out.sample_labels_config)

    create_sample_specific_tracks(sid_specific_ch)
    create_sample_specific_tracks.out |view



    // test3 = Channel.fromPath( './test3.csv' ).splitCsv() | read12_tester | view
    //     // .mix(diffreps_summary.out.diffreps_track)| view

    // read12_tester

    // log.info """\
    //      --
    // run as       : ${workflow.commandLine}
    // user       : ${workflow.userName}
         
    //      config files : ${workflow.configFiles}
         
    //      """
    //      .stripIndent()

    // parse_configuration_xls(params.xlsx_config)
    // parse_configuration_xls.out.sample_labels_config | view
    // parse_configuration_xls.out.diffreps_config | view
    // parse_configuration_xls.out.fastq_config | view

}

