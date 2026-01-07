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

// method, window
// for CSAW and DESEQ2 no window parameters
methods = Channel.from(["DESEQ2"],
                       ["CSAW"],
                       ["DIFFREPS", 200],
                       ["DIFFREPS", 1000],
                       ["RIPPM", 200],
                       ["RIPPM", 1000])

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
include {DESEQ2_MUMERGE} from './subworkflow/deseq2_mumerge.nf'
include {CSAW_MUMERGE} from './subworkflow/csaw_mumerge.nf'
include {DIFFREPS} from './subworkflow/diffreps_mumerge.nf'
include {QUALITY_PCA} from './subworkflow/quality_pca.nf'
include {COMBINE_HIST_PDF} from './subworkflow/combine_hist_pdf.nf'
include {CREATE_HISTOGRAMS} from './subworkflow/create_histograms.nf'
include {CREATE_BYCOMPARISON_TRACKS} from './subworkflow/create_bycomparison_tracks.nf'
include {RIPPM_NORM_FACTORS} from './subworkflow/rippm_norm_factors.nf'
include {CREATE_BIGWIG_FILES} from './subworkflow/create_bigwig_files.nf'
include {CREATE_TRACKLINES} from './subworkflow/create_tracklines.nf'

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
    
    cpus 16
    memory '32 GB'
    beforeScript 'source $HOME/.bashrc'

    // put to cache without overwrite
    publishDir path: "${params.chipseq_bam_cache}/${sample_id}/bam/", mode: 'copy', pattern: "${sample_id}_sorted_filtered.bam*", overwrite: false
    publishDir path: "${params.chipseq_bam_cache}/${sample_id}/bam/", mode: 'copy', pattern: "${sample_id}.bowtie2.log", overwrite: false
    
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
   
    input:
    path(files)
    path(track_lines)
    path(combined_track_lines)

    script:
    data_path="/net/waxman-server/mnt/data/waxmanlabvm_home/${workflow.userName}/${params.dataset_label}"
    """
    mkdir -p ${data_path}/TRACK_LINES
    cp ${files} ${data_path}
    cp ${track_lines} ${combined_track_lines} ${data_path}/TRACK_LINES/
    """
}

process collect_metrics {
    tag "${sample_id}"
    cpus 1
    memory '16 GB'
    
    storeDir "${params.chipseq_bam_cache}/${sample_id}/metrics/"

    beforeScript 'source $HOME/.bashrc'
    
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
    
    input:
    tuple val(sample_id), 
          val(library), 
          val(downsample), 
          val(r1), 
          val(r2),
          path(bam_file),      
          path(bai_file),
          path(log_file)
    
    output:
    tuple val(sample_id), val(lib), path(bam_file), path(bai_file), emit: bam
    tuple val(sample_id), path(log_file), emit: log
    
    script:
    lib = library.toString().toLowerCase()
    """
    """
}

// aggregate individual group reports
process by_comparison_overlap {
    tag "${group_name}"
    executor 'local'
    publishDir path: "${params.output_dir}/summary/by_comparison_overlap/", mode: "copy", pattern: "*_overlap.xlsx", overwrite: true
    publishDir path: "${params.output_dir}/summary/by_comparison_overlap/${group_name}_TOP25", mode: "copy", pattern: "*TOP25*.xlsx", overwrite: true
    beforeScript 'source $HOME/.bashrc'
    input:
    tuple val(group_name), path(diffreps_reports)
    
    output:
    tuple val(group_name), path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    methods_overlap.R --output_prefix ${group_name}
    """
}

//aggregate all diffreps and manorm2 reports
process all_comparisons_overlap {
    tag "all_comparisons_overlap"

    executor 'sge'
    cpus 4
    memory '32 GB'
    time '1h'
    // executor 'local'
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "*_overlap.xlsx", overwrite: true
    publishDir path: "${params.output_dir}/summary/all_groups_together_overlap_TOP25/", mode: "copy", pattern: "*TOP25*.xlsx", overwrite: true

    beforeScript 'source $HOME/.bashrc'

    input:
    path(diffreps_reports)
    
    output:
    path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    methods_overlap.R --output_prefix "all_groups_together"
    """
}


process extradetailed_macs2_epic_peaks_overlaps {
    
    executor 'sge'
    cpus 4
    memory '32 GB'
    errorStrategy 'retry'
    maxRetries 3
    
    beforeScript 'source $HOME/.bashrc'
    
    publishDir path: "${params.output_dir}/summary/extradetailed_macs2_epic_peaks_overlaps/", mode: "copy", pattern: "*.xlsx", overwrite: true

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


def parse_diffreps_configuration(row) {
    def meta = [:]
    // 1, Hnf6_Male, Hnf6_Female, G73M03|G73M04|G76M12|G76M13, G73M01|G73M02|G76M10|G76M11, RIPPM, 1000
    meta.num               = row[0].trim()
    meta.treatment_name    = row[1].trim()
    meta.control_name      = row[2].trim()
    meta.treatment_samples = row[3].trim()
    meta.control_samples   = row[4].trim()
    meta.method            = row[5].trim()
    meta.window_size = (6 < row.size()) ? row.get(6) : null

    window_size_str = meta.window_size == null ? "" : "_${meta.window_size}"
    meta.report_name = "${meta.num}_${meta.method}${window_size_str}_${meta.treatment_name}_vs_${meta.control_name}"
    meta.group_name = "${meta.num}_${meta.treatment_name}_vs_${meta.control_name}"
    return meta
}

def parse_comparisons_configuration(row) {
    def meta = [:]
    // 1, Hnf6_Male, Hnf6_Female, G73M03|G73M04|G76M12|G76M13, G73M01|G73M02|G76M10|G76M11, RIPPM, 1000
    meta.num               = row[0].trim()
    meta.treatment_name    = row[1].trim()
    meta.control_name      = row[2].trim()
    meta.treatment_samples = row[3].trim()
    meta.control_samples   = row[4].trim()
    meta.group_name = "${meta.num}_${meta.treatment_name}_vs_${meta.control_name}"
    return meta
}


workflow {

    fastq_config_ch = Channel.from(params.fastq_config)
    sample_labels_config_ch = Channel.from(params.sample_labels_config)
    diffreps_config_ch = Channel.from(params.diffreps_config)

    comparisons = diffreps_config_ch
        .splitCsv()
        .map{it->parse_comparisons_configuration(it)}

    diffreps_final = diffreps_config_ch
        .splitCsv()
        .combine(methods)
        .map{it->parse_diffreps_configuration(it)}
        .branch {
            diffreps: it.method =~"DIFFREPS|RIPPM"
            deseq2: it.method =~"DESEQ2"
            csaw: it.method =~"CSAW"
        }

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

    fastq_for_mapping = trim_adapters(bam_cache_status.uncached.map{it-> it[0..4]}) 
    
    // Make QC analysis (only for uncached, it will work because of storeDir)
    // fastqc(fastq_for_mapping) // calculate fastqc only for new samples with trimmed adapters
    fastqc(fastq_records)

    // Calculate number of reads in fastq files  (only for uncached)
    fq_num_reads = fastqc.out.raw_reads  
        .map{it -> it[1]}
        .collectFile(name: "numreads.csv", keepHeader: true)

    bowtie2_align(fastq_for_mapping, mm9_black_complement)
    
    // Combine newly generated BAMs with cached BAMs
    cached_bams_ch = bam_cache_status.cached
        .map { sample_id, library, downsample, r1, r2, exist  ->
            tuple(
                sample_id, 
                library, 
                downsample, 
                r1, 
                r2,
                file("${params.chipseq_bam_cache}/${sample_id}/bam/${sample_id}_sorted_filtered.bam"),
                file("${params.chipseq_bam_cache}/${sample_id}/bam/${sample_id}_sorted_filtered.bam.bai"),
                file("${params.chipseq_bam_cache}/${sample_id}/bam/${sample_id}.bowtie2.log")
            )
        } 
    
    retrieve_cached_bams(cached_bams_ch)
    all_bams = bowtie2_align.output.bam.mix(retrieve_cached_bams.out.bam)
    all_bowtie2_logs = bowtie2_align.output.log.mix(retrieve_cached_bams.out.log)
    
    bam_count(all_bams)
    
    macs2_callpeak(all_bams, mm9_chrom_sizes)
    epic2_callpeak(bam_count.out.fragments_bed6, mm9_chrom_sizes)

    //picard
    collect_metrics(all_bams)
    
    // calculate overlaps for all MACS2 narrow/broad and SICER peaks
    // for each type of peaks create separate xlsx report
    peaks_for_aggregation = macs2_callpeak.out.xls
        .mix(macs2_callpeak.out.broad_xls)
        .mix(epic2_callpeak.out.bed6)
        .map{sid,bedfile -> bedfile}
        .collect()

    extradetailed_macs2_epic_peaks_overlaps(peaks_for_aggregation, params.sample_labels_config)
    
    // log.info("params.peakcaller: $params.peakcaller")

    if (params.peakcaller == "MACS2") {
        extra_columns = macs2_callpeak.out.xls // for diffreps summary
        quality_control_peaks = macs2_callpeak.out.xls
        mumerge_peaks = macs2_callpeak.out.narrow_bed
    } else {
        extra_columns = epic2_callpeak.out.bed6 // for diffreps summary
        quality_control_peaks = epic2_callpeak.out.bed6
        mumerge_peaks = epic2_callpeak.out.bed6
    }

    MUMERGE(comparisons, mumerge_peaks)

    DESEQ2_MUMERGE(
        diffreps_final.deseq2,
        all_bams,
        MUMERGE.out.mumerge_peaks // mumerge peaks instead of MACS2 union 
    )

    CSAW_MUMERGE(
        diffreps_final.csaw,
        all_bams,
        MUMERGE.out.mumerge_peaks // mumerge peaks instead of MACS2 union 
    )
        
    RIPPM_NORM_FACTORS(MUMERGE.out.mumerge_overlap, bam_count.out.fragments)

    DIFFREPS(
        diffreps_final.diffreps,
        bam_count.out.fragments_bed6,
        RIPPM_NORM_FACTORS.out.norm_factors,
        fq_num_reads, // table with number of reads in R1.fq files for each sample
        MUMERGE.out.mumerge_peaks // mumerge peaks instead of MACS2 union 
    )

    QUALITY_PCA(
        comparisons,
        quality_control_peaks,
        params.peakcaller
    )

    // individual by comparison reports
    by_comparisons_report_ch = DIFFREPS.out.full_report
        .mix(DESEQ2_MUMERGE.out.full_report)
        .mix(CSAW_MUMERGE.out.full_report)
        .map{meta, rest -> [meta.group_name, rest]}
        .groupTuple()

    by_comparison_overlap(by_comparisons_report_ch)

    // Aggregated report which contains all comparisons together
    only_reports_ch  = DIFFREPS.out.full_report
        .mix(DESEQ2_MUMERGE.out.full_report)
        .mix(CSAW_MUMERGE.out.full_report)
        .map{meta,rest -> rest}
        .collect()

    all_comparisons_overlap(only_reports_ch)

    // Create histograms
    combined_diffpeak_reports_ch = DIFFREPS.out.full_report
        .mix(DESEQ2_MUMERGE.out.full_report)
        .mix(CSAW_MUMERGE.out.full_report)

    CREATE_HISTOGRAMS(combined_diffpeak_reports_ch)

    CREATE_BIGWIG_FILES(
        RIPPM_NORM_FACTORS.out.norm_factors,
        bam_count.out.fragments,
        sample_labels_config_ch,
        mm9_chrom_sizes
    )

    CREATE_BYCOMPARISON_TRACKS(
        combined_diffpeak_reports_ch,
        MUMERGE.out.mumerge_peaks,
        mm9_chrom_sizes)
    
    // CREATE_BYCOMPARISON_TRACKS.out.tracks | view
    //     .mix(MUMERGE.out.mumerge_peaks)
    //     .groupTuple() | view

    CREATE_TRACKLINES(
        CREATE_BIGWIG_FILES.out.individual
            .mix(macs2_callpeak.out.broad_bb)
            .mix(macs2_callpeak.out.narrow_bb)
            .mix(epic2_callpeak.out.epic2_bb), // Individual tracks by sample_id

        CREATE_BIGWIG_FILES.out.grouped, //Bigwig combined for each group 
        
        CREATE_BYCOMPARISON_TRACKS.out.tracks, //DIFFREPS,CSAW,..,MUMERGE by comparison

        sample_labels_config_ch // for proper labels
    )
       
    track_lines = Channel.from(default_tracks)
        .concat(CREATE_TRACKLINES.out.comparison_specific_tracklines,
                CREATE_TRACKLINES.out.individual_tracklines)
        .collectFile(name: 'autolimit_tracks.txt', sort: 'index'){item -> item.text}

    combined_bw_track_lines = Channel.from(default_tracks)
        .concat(CREATE_TRACKLINES.out.comparison_specific_tracklines,
                CREATE_TRACKLINES.out.group_combined_tracklines)
        .collectFile(name: 'autolimit_combined_tracks.txt', sort: 'index'){item -> item.text}

    track_files_to_server = CREATE_BIGWIG_FILES.out.individual.map{it -> it[1]}
        .concat(
            epic2_callpeak.out.epic2_bb.map{it -> it[1]},
            macs2_callpeak.out.narrow_bb.map{it -> it[1]},
            macs2_callpeak.out.broad_bb.map{it -> it[1]},
            CREATE_BIGWIG_FILES.out.grouped.map{it->it[3]},
            CREATE_BYCOMPARISON_TRACKS.out.tracks.map{it->it[2]}
        )
        .collect()

    if (params.copy_to_server_bool){
        copy_files_to_server(track_files_to_server,
                             track_lines,
                             combined_bw_track_lines)
    }


    multiqc(
        fastqc.out.zip.map{it -> it[1]}.collect(),
        all_bowtie2_logs.map{it -> it[1]}.collect(),
        bam_count.out.stats.map{it -> it[1]}.collect(),
        macs2_callpeak.out.xls.map{it -> it[1]}.collect(),
        collect_metrics.out.metrics.map{it -> it[1]}.collect()
    )
}

