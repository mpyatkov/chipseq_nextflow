#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// TODO: use cleanup = true inside nextflow.config to remove the work directory after run 
//>>> PARAMETERS TO CHANGE
params.copy_to_server_bool=false
params.trim_adapters = true // does not have any effect, trimming by default // TODO: deprecated
params.peakcaller="MACS2"    // by default going to use MACS2 (SICER/EPIC2 alternative)
// ChIPSEQ TF                                --> MACS2 model will be used
// ATAC (DAR), CutNRun, Histone modification --> MACS2 model is not required
params.macs2_model = false
//<<< PARAMETERS TO CHANGE

mm9_chrom_sizes  = file("$projectDir/assets/mm9.chrom.sizes", checkIfExists: true)
mm9_black_complement  = file("$projectDir/assets/mm9-blacklist_complement", checkIfExists: true)
default_tracks  = file("$projectDir/assets/default_tracks.txt", checkIfExists: true)

params.bowtie2_index="/projectnb/wax-es/aramp10/Bowtie2/Bowtie2Index/genome"
params.chipseq_bam_cache="/projectnb/wax-es/CHIPSEQ_BAM_CACHE"
params.dataset_label="TEST1" // default dataset label if not provided
params.rversion="4.2"
params.output_dir="./RESULTS_${params.dataset_label}"
params.fastq_config = file("$projectDir/${params.input_configs}/fastq_config.csv", checkIfExists: true)
params.sample_labels_config = file("$projectDir/${params.input_configs}/sample_labels.csv", checkIfExists: true)
params.diffreps_config = file("$projectDir/${params.input_configs}/diffreps_config.csv")
params.overwrite_outputs = true

include {MUMERGE} from './subworkflow/mumerge.nf'
include {DIFFREPS} from './subworkflow/diffreps/diffreps.nf'
include {MANORM2} from './subworkflow/diffreps/manorm2.nf'
include {QUALITY_PCA} from './subworkflow/quality_pca.nf'
include {COMBINE_HIST_PDF} from './subworkflow/combine_hist_pdf.nf'

process trim_adapters {
    tag "${sample_id}"

    cpus 8
    memory '8 GB'

    beforeScript 'source $HOME/.bashrc'

    input:
    tuple val(sample_id), val(library), val(downsample), val(r1), val(r2)

    output:
    tuple val(sample_id), val(library), val(downsample), path("${sample_id}_val_1.fq.gz"), path("${sample_id}_val_2.fq.gz")

    script:

    lib = library.toString().toLowerCase()
    if ( lib == "paired") {
        """
        module load trimgalore
        module load cutadapt

        ## remove autodetected adapters
        trim_galore --gzip --stringency 13 --trim1 --length 30 --quality 0 --paired $r1 $r2 -j 4 --basename ${sample_id}
        """
    } else {
        """
        module load trimgalore  
        module load cutadapt

        trim_galore --gzip --stringency 13 --trim1 --length 30 --quality 0 $r1 -j 4 --basename ${sample_id}
        # create dummy file for output consistency
        mv ${sample_id}_trimmed.fq.gz ${sample_id}_val_1.fq.gz
        touch ${sample_id}_val_2.fq.gz
        """
    }
}

process bowtie2_align {

    tag "${sample_id}"
    
    //echo true
    // executor "local"
    cpus 16
    memory '32 GB'
    debug true
    beforeScript 'source $HOME/.bashrc'

    // put to cache without overwrite

    publishDir path: "${params.chipseq_bam_cache}/${sample_id}/bam/", mode: 'copy', pattern: "${sample_id}_sorted_filtered.bam*", overwrite: false
    publishDir path: "${params.chipseq_bam_cache}/${sample_id}/bam/", mode: 'copy', pattern: "${sample_id}.bowtie2.log", overwrite: false
    // publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/bam/", mode: "symlink", pattern: "${sample_id}_sorted_filtered.bam*", overwrite: params.overwrite_outputs 
    // publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/bam/", mode: "copy", pattern: "*.log", overwrite: true
    
    input:
    tuple val(sample_id), val(library), val(downsample), path(r1), path(r2)
    file(anti_blacklist)
    
    output:
    tuple val(sample_id), val(lib), path("${sample_id}_sorted_filtered.bam"), path("${sample_id}_sorted_filtered.bam.bai"),  emit: bam
    tuple val(sample_id), path("*.log"), emit: log
    
    script:

    reads_args = null
    filter_bedpe = null
    lib = library.toString().toLowerCase()
    if ( lib == "paired") {
        reads_args = "-1 ${r1} -2 ${r2}"
        filter_bedpe = "and proper_pair"
    } else {
        reads_args = "-U ${r1}"
        filter_bedpe = ""
    }

    """
    module load bowtie2
    module load samtools

    set -x
    bowtie2 -p $task.cpus -x $params.bowtie2_index $reads_args -S ${sample_id}_bowtie2.sam 2> ${sample_id}.bowtie2.log
    samtools view -@ $task.cpus -b ${sample_id}_bowtie2.sam > ${sample_id}_alignments.bam
    samtools sort -@ $task.cpus ${sample_id}_alignments.bam -o "sorted.bam"
    samtools index -@ $task.cpus "sorted.bam"

    rm -rf ${sample_id}_alignments.bam ${sample_id}_bowtie2.sam

    ## removing blacklist regions from bam
    samtools view -@ $task.cpus -L ${anti_blacklist} -O BAM -o "${sample_id}_sorted.bam" sorted.bam
    samtools index -@ $task.cpus "${sample_id}_sorted.bam"

    ## filtering from low quality, not mapped, not paired
    sambamba view -h -t $task.cpus -f bam -F "[XS] == null and not unmapped ${filter_bedpe}" "${sample_id}_sorted.bam" > 1.bam
    
    sambamba sort -t $task.cpus 1.bam
    mv 1.sorted.bam ${sample_id}_sorted_filtered.bam    
    sambamba index ${sample_id}_sorted_filtered.bam

    ## downsampling if it is required
    if [[ "${downsample}" != "NO" ]]; then
       mv ${sample_id}_sorted_filtered.bam 1.bam
       mv ${sample_id}_sorted_filtered.bam.bai 1.bam.bai

       ## ## remove symlinks 
       ## rm ${sample_id}_sorted_filtered.bam
       ## rm ${sample_id}_sorted_filtered.bam.bai

       ## before downsampling sorting by name
       sambamba sort -n -t $task.cpus 1.bam
       samtools view -b -s ${downsample} 1.sorted.bam > tmp.bam

       ## sorting back by coordinates
       sambamba sort -t $task.cpus tmp.bam
       mv tmp.sorted.bam ${sample_id}_sorted_filtered.bam
       sambamba index ${sample_id}_sorted_filtered.bam
    fi
    
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
    
    // debug true
    cpus 4
    memory '32 GB'
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/bam/", mode: "symlink", pattern: "*fragments*.bed*", overwrite: params.overwrite_outputs 
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/bam/", mode: "symlink", pattern: "${sample_id}_sorted_filtered.bam*", overwrite: params.overwrite_outputs 

    storeDir "${params.chipseq_bam_cache}/${sample_id}/bam_count/"
    
    input:
    tuple val(sample_id), val(lib), path(bam), path(bai)
    
    output:
    
    tuple val(sample_id), path("${sample_id}_fragments.bed"), emit: fragments
    tuple val(sample_id), path("${sample_id}_fragments_bed6.bed"), emit: fragments_bed6
    tuple val(sample_id), path("*stats"), emit: stats

    script:
    bedtools_bedpe = lib == "paired" ? "-bedpe" : ""

    """
    module load bedtools
    module load samtools
    
    ## extracting fragments
    ## for paired-end (start1,end1 - start2,end2) it will be difference (start1,end2)
    ## for single-end it will be the read itself
    if [[ ${lib} == "paired" ]]; then

        ## it is required to sort bam file by name for paired-end
        ## bedtools does not work properly with coordinate sorted bam files
        sambamba sort -n -t $task.cpus -o tmp.bam ${bam} 
        bedtools bamtobed ${bedtools_bedpe} -i tmp.bam 2>/dev/null | cut -f 1,2,6 > ${sample_id}_fragments.bed
        rm -rf tmp.bam
    else
        bedtools bamtobed ${bedtools_bedpe} -i ${bam} 2>/dev/null > ${sample_id}_fragments.bed 2> /dev/null
    fi

    awk -v OFS='\t' '{print \$1,\$2,\$3,".","0","."}' ${sample_id}_fragments.bed > ${sample_id}_fragments_bed6.bed

    ## calc some stats
    samtools stats -@ $task.cpus ${bam} > ${sample_id}_sorted_filtered.stats
    samtools flagstats -@ $task.cpus ${bam} > ${sample_id}_sorted_filtered.flagstats

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
    
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/macs2/", mode: "copy", pattern: "*.{narrowPeak,broadPeak,xls,bed}", overwrite: true

    storeDir "${params.chipseq_bam_cache}/${sample_id}/macs2/"
    
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
    lib=library == "paired" ? "BAMPE" : "BAM"
    model=params.macs2_model ? "" : "--nomodel"
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
    // echo true
    
    beforeScript 'source $HOME/.bashrc'
    
    // publishDir path: "${params.output_dir}/peak_union/", mode: "copy", pattern: "peak_union.bed", overwrite: true

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
    // cache false // TODO: remove
    
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/macs2/", mode: "copy", pattern: "${sample_id}_fragments_union_coverage.bed", overwrite: true
        
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
    // echo true
    
    input:
    tuple val(sample_id), path(fragments), path(fragments_union_coverage)
    
    output:
    path("${sample_id}_output.txt")
    script:

    """
    FRAGMENT_COUNT=\$(cat $fragments | wc -l )
    FRAGMENT_IN_PEAK_COUNT=\$(awk '{n+=\$4;} ; END {print n;}' ${fragments_union_coverage})
    FRAGMENT_IN_PEAK_RATIO=\$(echo "scale=4;\${FRAGMENT_IN_PEAK_COUNT}/\${FRAGMENT_COUNT}" | bc)
    echo "sample_id,nfragments,nfragments_in_peak,ratio" > ${sample_id}_output.txt

    echo "${sample_id},\${FRAGMENT_COUNT},\${FRAGMENT_IN_PEAK_COUNT},\${FRAGMENT_IN_PEAK_RATIO}" >> ${sample_id}_output.txt
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

process create_group_combined_tracks {
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

process create_sample_specific_tracks {
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
    module load R/${params.rversion}
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
//     module load R/${params.rversion}
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
    
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/epic2/", mode: "copy", pattern: "*bed", overwrite: true

    storeDir "${params.chipseq_bam_cache}/${sample_id}/epic2/"
    
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
    
    ## removing header
    sed -i '1d' "${sample_id}_epic_bed6.bed" 

    ## convert to bed3
    awk 'OFS="\t" {print \$1,\$2,\$3}' "${sample_id}_epic_bed6.bed" > "${sample_id}_epic_bed3.bed"

    ## filtering out junk and creating bigbed track
    cat "${sample_id}_epic_bed3.bed" | grep -vE "track|chrM|random" > tmp.bed
    bedtools sort -i tmp.bed > tmp.sorted.bed
    bedClip tmp.sorted.bed ${mm9_chrom_sizes} tmp.sorted.clipped.bed
    bedToBigBed -allow1bpOverlap tmp.sorted.clipped.bed ${mm9_chrom_sizes} "${sample_id}_epic2.bb"
    """
}

process copy_files_to_server {
    executor 'local'
    // publishDir path: "${data_path}", mode: "copy", pattern: "*.{bw,bam*,bb}", overwrite: false
    // publishDir path: "${data_path}/TRACK_LINES", mode: "copy", pattern: "*.txt", overwrite: true
    
    input:
    path(files)
    path(track_lines)
    path(combined_track_lines)
    path(track_hub)

    // output:
    // path("*.{bw,bam,bai,bb,txt}")
    // stdout
    
    script:
    data_path="/net/waxman-server/mnt/data/waxmanlabvm_home/${workflow.userName}/${params.dataset_label}"
    """
    mkdir -p ${data_path}/TRACK_LINES
    cp ${files} ${data_path}
    cp ${track_lines} ${track_hub} ${combined_track_lines} ${data_path}/TRACK_LINES/
    """
}


process collect_metrics {
    tag "${sample_id}"
    cpus 1
    memory '16 GB'
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/metrics/", mode: "copy", overwrite: true

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
    module load R/${params.rversion}

    # set -x
    
    picard \
        CollectMultipleMetrics \
        --INPUT $bam \
        --OUTPUT ${sample_id}.CollectMultipleMetrics

    ## only for paired-end
    if [[ $library == "paired" ]];
    then
    picard \
        CollectInsertSizeMetrics \
        --Histogram_FILE ${sample_id}.InsertSizeMetrics.pdf \
        --INPUT ${bam} \
        --OUTPUT ${sample_id}.InsertSizeMetrics
    fi
    """
}

process fastq_num_reads {
    tag "${sample_id}"
    cpus 1
    time '1h'

    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/metrics/fastqc/", mode: "copy", overwrite: params.overwrite_outputs 

    input:
    tuple val(sample_id), val(library), val(downsample), val(r1), val(r2)

    output:
    tuple val(sample_id), path("*_raw_reads.txt") , emit: raw_reads

    script:
    """
    ## calculate number of raw reads
    NUMREADS=`echo \$(zcat $r1 | wc -l)/4 | bc`
    echo "sample_id,num_raw_reads" >> ${sample_id}_raw_reads.txt
    echo "${sample_id},\${NUMREADS}" >> ${sample_id}_raw_reads.txt
    """
}

process fastqc {
    tag "${sample_id}"
    cpus 4
    time '3h'
    memory '32 GB'
    errorStrategy 'retry'
    maxRetries 3
    
    // publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/metrics/fastqc/", mode: "symlink", overwrite: params.overwrite_outputs 
    publishDir path: "${params.output_dir}/SAMPLES/${sample_id}/metrics/fastqc/", mode: "symlink", overwrite: params.overwrite_outputs 
    storeDir "${params.chipseq_bam_cache}/${sample_id}/fastqc/"

    
    beforeScript 'source $HOME/.bashrc'
    
    // echo true
    input:
    tuple val(sample_id), val(library), val(downsample), val(r1), val(r2)
    
    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip") , emit: zip
    tuple val(sample_id), path("*_raw_reads.txt") , emit: raw_reads
    
    script:
    lib = library.toString().toLowerCase()

    """
    module load fastqc

    if [[ "${lib}" == "paired" ]]; then
        
        mkdir -p ${sample_id}_out
        ln -s $r1 ${sample_id}_1.fq.gz
        ln -s $r2 ${sample_id}_2.fq.gz
        fastqc --threads $task.cpus $r1 $r2
        fastqc -o ${sample_id}_out --threads $task.cpus ${sample_id}_1.fq.gz ${sample_id}_2.fq.gz
        mv ${sample_id}_out/* ./
    else
        mkdir -p ${sample_id}_out
        ln -s $r1 ${sample_id}_1.fq.gz
        fastqc -o ${sample_id}_out --threads $task.cpus ${sample_id}_1.fq.gz
        mv ${sample_id}_out/* ./
    fi

    ## calculate number of raw reads
    NUMREADS=`echo \$(zcat $r1 | wc -l)/4 | bc`
    echo "sample_id,num_raw_reads" >> ${sample_id}_raw_reads.txt
    echo "${sample_id},\${NUMREADS}" >> ${sample_id}_raw_reads.txt
    """
}

process multiqc {

    cpus 1
    publishDir path: "${params.output_dir}/summary/multiqc/", mode: "copy", pattern: "multiqc_report.html", overwrite: params.overwrite_outputs 
    publishDir path: "/net/waxman-server/mnt/data/waxmanlabvm_home/${workflow.userName}/${params.dataset_label}/multiqc/", mode: "copy", pattern: "multiqc_report.html", overwrite: params.overwrite_outputs 
    
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

def is_empty_file(fp) {
    File file = new File(fp);
    return !file.exists() || file.length() == 0
}

// Define the processes for cache operations
process check_bam_cache {

    executor "local"
    // cache false
    input:
    tuple val(sample_id), val(library), val(downsample), val(r1), val(r2)
    
    
    output:
    tuple val(sample_id), val(library), val(downsample), val(r1), val(r2), val(exists)
    
    exec:
    // Define the path to the cached BAM based on sample ID
    def cachePath = "${params.chipseq_bam_cache}/${sample_id}/bam/${sample_id}_sorted_filtered.bam"
    def cacheFile = new File(cachePath)
    exists = cacheFile.exists()
    
    """
    echo "Checking cache for ${sample_id}: ${exists ? 'FOUND' : 'NOT FOUND'}"
    """
}

process retrieve_cached_bams {

    executor 'local'
    // cache false
    // debug true
    input:
    tuple val(sample_id), val(library), val(downsample), val(r1), val(r2)
    
    output:
    tuple val(sample_id),
        val(lib),
        path("${sample_id}_sorted_filtered.bam"),
        path("${sample_id}_sorted_filtered.bam.bai"),  emit: bam

    tuple val(sample_id), path("${sample_id}.bowtie2.log"),  emit: log
    
    script:
    lib = library.toString().toLowerCase()
    
    """
    ln -s ${params.chipseq_bam_cache}/${sample_id}/bam/${sample_id}_sorted_filtered.bam .
    ln -s ${params.chipseq_bam_cache}/${sample_id}/bam/${sample_id}_sorted_filtered.bam.bai .
    ln -s ${params.chipseq_bam_cache}/${sample_id}/bam/${sample_id}.bowtie2.log .
    """
}

//aggregate all diffreps and manorm2 reports
process diffreps_manorm2_overlap_general {
    tag "diffreps_manorm2_overlap_general"

    executor 'sge'
    cpus 4
    memory '32 GB'
    time '1h'
    // executor 'local'
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "*overlap_manorm2_vs_diffreps.xlsx", overwrite: true
    beforeScript 'source $HOME/.bashrc'

    input:
    path(diffreps_reports)
    
    output:
    path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    diffreps_manorm2_overlap.R --output_prefix "all_groups_together"
    """
}


process extradetailed_peaks_overlaps {
    
    executor "local"
    // echo true
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/summary/extradetailed_peaks_overlaps/", mode: "copy", pattern: "*.xlsx", overwrite: true

    input:
    path(peaks)
    path(sample_labels)
    
    output:
    path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    detailed_macs2_peaks_overlap.R --peak_prefix "broad" --sample_labels ./sample_labels.csv
    detailed_macs2_peaks_overlap.R --peak_prefix "narrow" --sample_labels ./sample_labels.csv
    detailed_epic2_peaks_overlap.R --peak_prefix "epic" --sample_labels ./sample_labels.csv
    """
}

process combine_mn2_dr_pdfs {
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
    pdfunite $manorm2_pdf $diffreps_pdfs "${group_name}_Histogram_Barcharts_${params.peakcaller}.pdf"
    """
}

process combine_bigwig_tracks {
    tag "${group_name}"

    executor 'sge'
    cpus 8
    memory '16 GB'
    time '2h'

    beforeScript 'source $HOME/.bashrc'
    // echo true

    input:
    tuple val(group_name), val(samples_string), val(color), path(sample_paths)
    path(mm9_chrom_sizes)
    
    output:
    tuple val(group_name), val(samples_string), val(color), path("${group_name}.bw")

    script:
    """
    ## for wiggletools
    module load miniconda
    conda activate /projectnb/wax-es/routines/condaenv/deeptools
    bigwigAverage -p $task.cpus -b *.bw -o "${group_name}.bw"
    """
}

process diffreps_manorm2_overlap {
    tag "${group_name}"
    executor 'local'
    publishDir path: "${params.output_dir}/summary/manorm2_diffreps_overlap_individual_reports_for_groups/", mode: "copy", pattern: "*overlap_manorm2_vs_diffreps.xlsx", overwrite: true
    publishDir path: "${params.output_dir}/summary/manorm2_diffreps_overlap_top25/", mode: "copy", pattern: "*top25*.xlsx", overwrite: true
    
    beforeScript 'source $HOME/.bashrc'
    input:
    tuple val(group_name), path(diffreps_reports), path(manorm2_report)
    
    output:
    tuple val(group_name), path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    diffreps_manorm2_overlap.R --output_prefix ${group_name}
    """
}

// process cache_bams {

//     executor 'local'
//     publishDir path: "${params.chipseq_bam_cache}/${sample_id}/bam/", mode: 'copy', pattern: "${sample_id}_sorted_filtered.bam*", overwrite: false

//     input:
//     tuple val(sample_id), val(library), path(bam), path(bai)
    
//     output:
//     tuple val(sample_id), path("${sample_id}_sorted_filtered.bam"), path("${sample_id}_sorted_filtered.bam.bai"),  emit: bam
    
//     script:
//     """
//     cp ${bam} ${sample_id}_sorted_filtered.bam
//     cp ${bai} ${sample_id}_sorted_filtered.bam.bai
//     """
// }

workflow {
    // parse_configuration_xls(params.xlsx_config)
    // parse_configuration_xls.out.sample_labels_config | view
    // parse_configuration_xls.out.diffreps_config | view
    // test = parse_configuration_xls.out.diffreps_config.ifEmpty(false)
    // .ifEmpty{exit 1, "Cannot find any input RDS files for Module 3"}

    fastq_config_ch = Channel.from(params.fastq_config)
    sample_labels_config_ch = Channel.from(params.sample_labels_config)
    diffreps_config_ch = Channel.from(params.diffreps_config)

    // [sid, library, downsample, r1, r2]
    fastq_records = fastq_config_ch.splitCsv()

    // [sid, library, downsample, r1, r2, cached:true|false]
    // add "cached" column
    cached_bams = check_bam_cache(fastq_records)

    // split bams by cached and not cached
    bam_cache_status = cached_bams.branch {
        cached: it[5] == true
        uncached: it[5] == false
    }

    // Process only uncached samples through the trimming step
    // TODO: revised step, I suppose we are going to use only trimmed fastq files, see below
    // if (params.trim_adapters) {
    //     fastq_for_mapping = trim_adapters(bam_cache_status.uncached.map{it-> it[0..4]})
    // } else {
    //     fastq_for_mapping = bam_cache_status.uncached.map{it -> it[0..4]}
    // }

    fastq_for_mapping = trim_adapters(bam_cache_status.uncached.map{it-> it[0..4]}) 
    
    // Make QC analysis (only for uncached, it will work because of storeDir)
    // fastqc(fastq_records)  
    fastqc(fastq_for_mapping) // calculate fastqc only for new samples with trimmed adapters 

    // Calculate number of reads in fastq files  (only for uncached)
    // (required for downstream report)
    // fastq_num_reads(fastq_for_mapping)

    fq_num_reads = fastqc.out.raw_reads  
        .map{it -> it[1]}
        .collectFile(name: "numreads.csv", keepHeader: true)
        
    bowtie2_align(fastq_for_mapping, mm9_black_complement)

    // TODO: make caching optional because for public data probably it will not
    // be required to cache bam files
    // if (params.need_cache) {
    //     cache_bams(bowtie2_align.out.bam)
    // } 
    
    // Combine newly generated BAMs with cached BAMs
    retrieve_cached_bams(bam_cache_status.cached.map { it[0..4] })
    all_bams = bowtie2_align.output.bam.mix(retrieve_cached_bams.out.bam)
    all_bowtie2_logs = bowtie2_align.output.log.mix(retrieve_cached_bams.out.log)
    
    bam_count(all_bams)
    
    macs2_callpeak(all_bams, mm9_chrom_sizes)
    epic2_callpeak(bam_count.out.fragments_bed6, mm9_chrom_sizes)
    
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

    // TODO: need to remove .combine(union) because we do not use it
    peaks_for_stats_ch = bam_count.out.fragments
        .join(fragments_union_overlap.out.fr_union_overlap)
    
    calc_sample_stats(peaks_for_stats_ch)
    
    sample_stats = calc_sample_stats.out
        .collectFile(name: "file.csv", keepHeader: true)
        // .collectFile{item -> item.text, keepHeader: true}
    sample_stats | calc_norm_factors

    // log.info("params.peakcaller: $params.peakcaller")

    if (params.peakcaller == "MACS2") {
        extra_columns = macs2_callpeak.out.xls // for diffreps summary
        peaks_for_manorm2 = macs2_callpeak.out.narrow_bed
        quality_control_peaks = macs2_callpeak.out.xls
        mumerge_peaks = macs2_callpeak.out.narrow_bed
    } else {
        extra_columns = epic2_callpeak.out.bed6 // for diffreps summary
        peaks_for_manorm2 = epic2_callpeak.out.bed6
        quality_control_peaks = epic2_callpeak.out.bed6
        mumerge_peaks = epic2_callpeak.out.bed6
    }

    MUMERGE(diffreps_config_ch, mumerge_peaks) 
    DIFFREPS(
        // parse_configuration_xls.out.diffreps_config,
        diffreps_config_ch,
        // parse_configuration_xls.out.sample_labels_config,
        sample_labels_config_ch,
        bam_count.out.fragments_bed6,
        calc_norm_factors.out,
        extra_columns, //macs2_callpeak.out.xls
        mm9_chrom_sizes,
        fq_num_reads, // table with number of reads in R1.fq files for each sample
        MUMERGE.out.mumerge_peaks // mumerge peaks instead of MACS2 union 
    )

    MANORM2(
        diffreps_config_ch,
        bam_count.out.fragments,
        peaks_for_manorm2,
        mm9_chrom_sizes,
        MUMERGE.out.mumerge_peaks // mumerge peaks instead of MACS2 union 
    )
    
    QUALITY_PCA(
        diffreps_config_ch,
        quality_control_peaks,
        params.peakcaller
    )


    // Aggregate DIFFREPS and MANORM2 reports for each INDIVIDUAL group separately
    manorm2_group_report_ch = MANORM2.out.diff_table
        .map{meta, rest -> [meta.group_name, rest]}
    
    diffreps_group_report_ch = DIFFREPS.out.full_report
        .map{meta, rest -> [meta.group_name, rest]}
        .groupTuple() 

    combined_manorm2_diffreps_ch = diffreps_group_report_ch.join(manorm2_group_report_ch)
    diffreps_manorm2_overlap(combined_manorm2_diffreps_ch)


    // Aggregated report which contains all DIFFREPS and MANORM2 reports together
    only_reports_ch  = DIFFREPS.out.full_report
        .mix(MANORM2.out.diff_table)
        .map{meta,rest -> rest}
        .collect()
        
    diffreps_manorm2_overlap_general(only_reports_ch)

    // Combine MAnorm2 and diffReps pdfs
    COMBINE_HIST_PDF (
        DIFFREPS.out.aggregated_histograms_diffreps_noxy,    
        DIFFREPS.out.aggregated_histograms_diffreps_allchr,
        MANORM2.out.manorm2_histogram_noxy,
        MANORM2.out.manorm2_histogram_allchr
    )
  
    // TRACKS (copying, generating track lines)
    // create_bigwig_files 
    sid_normfact_ch = calc_norm_factors.out.splitCsv(sep: "\t")
        .map{it->[it[0],it[4]]}
    
    sid_fr_norm = bam_count.out.fragments
        .join(sid_normfact_ch)

    // Combine bigWig files for one group into one track
    create_bigwig_files(sid_fr_norm, mm9_chrom_sizes)
    bigwig_group_ch = sample_labels_config_ch
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

    combine_bigwig_tracks(bigwig_group_ch, mm9_chrom_sizes)
    combined_bigwig_ch = combine_bigwig_tracks.out
        .map{group_name,samples_str,color,pth -> [group_name, samples_str, color, pth.getName()]}
        .collectFile{item -> item.join(",")+'\n'}    

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
        .combine(diffreps_config_ch) 

    create_sample_specific_tracks(sid_specific_ch)
    create_diffreps_tracks(group_specific_ch)
    create_group_combined_tracks(combined_bigwig_ch)

    track_lines = Channel.from(default_tracks)
        .concat(create_diffreps_tracks.out,create_sample_specific_tracks.out.sid_tracks)
        .collectFile(name: 'autolimit_tracks.txt', sort: 'index'){item -> item.text}

    combined_bw_track_lines = Channel.from(default_tracks)
        .concat(create_diffreps_tracks.out,create_group_combined_tracks.out)
        .collectFile(name: 'autolimit_combined_tracks.txt', sort: 'index'){item -> item.text}

    // copy bb,bam,bw files to server
    bw_files = create_bigwig_files.out.map{it -> it[1]} 
    narrow_files= macs2_callpeak.out.narrow_bb.map{it -> it[1]}
    broad_files= macs2_callpeak.out.broad_bb.map{it -> it[1]}
    broad_epic_files=epic2_callpeak.out.epic2_bb.map{it -> it[1]}
    bam_files=all_bams.map{it->[it[2],it[3]]}
    diffreps_files=DIFFREPS.out.diffreps_track.map{it->it[2]}
    manorm2_files=MANORM2.out.manorm2_track.map{it->it[2]}
    combined_bwfiles=combine_bigwig_tracks.out.map{it->it[3]} 

    // track_files_to_server = bw_files.concat(bam_files, broad_epic_files, broad_files, narrow_files, diffreps_files).collect() //excluded bam track files
    track_files_to_server = bw_files.concat(broad_epic_files, broad_files, narrow_files, diffreps_files, manorm2_files, combined_bwfiles).collect()

    if (params.copy_to_server_bool){
        copy_files_to_server(track_files_to_server,
                             track_lines,
                             combined_bw_track_lines,
                             create_sample_specific_tracks.out.bigwig_hub)
    }

    // calculate overlaps for all MACS2 narrow/broad and SICER peaks
    // for each type of peaks create separate xlsx report
    peaks_for_aggregation = macs2_callpeak.out.xls
        .mix(macs2_callpeak.out.broad_xls)
        .mix(epic2_callpeak.out.bed6)
        .map{sid,bedfile -> bedfile}
        .collect()

    extradetailed_peaks_overlaps(peaks_for_aggregation, params.sample_labels_config)
    
    //picard
    collect_metrics(all_bams)
    
    multiqc(
        fastqc.out.zip.map{it -> it[1]}.collect(),
        all_bowtie2_logs.map{it -> it[1]}.collect(),
        bam_count.out.stats.map{it -> it[1]}.collect(),
        macs2_callpeak.out.xls.map{it -> it[1]}.collect(),
        collect_metrics.out.metrics.map{it -> it[1]}.collect()
    )

    multiqc.out.report
}

