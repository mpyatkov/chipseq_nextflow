#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// process diffreps {
//     tag "${meta.num}_${meta.method}_${meta.window}"
    
//     executor 'sge'
//     cpus 16
//     //cache false
//     time '2h'
//     // echo true

//     beforeScript 'source $HOME/.bashrc'
    
//     input:
//     tuple val(meta), path(TREATMENT_FILES), path(CONTROL_FILES), path(NORM_FILE)
    
//     output:
//     tuple val(meta) , path("*_vs_*")

//     shell:
//     template 'diffreps.sh'
// }

process diffreps {
    tag "${meta.num}_${meta.method}_${meta.window_size}"
    
    executor 'sge'
    cpus 16
    //cache false
    time '2h'
    // echo true

    // debug true
    
    beforeScript 'source $HOME/.bashrc'
    
    input:
    tuple val(meta), path(TREATMENT_FILES), path(CONTROL_FILES), path(NORM_FILE)
    
    output:
    tuple val(meta) , path("*_vs_*")

    script:
    """
    set +eu
    module load bedtools
    module load miniconda
    conda activate /projectnb/wax-es/routines/condaenv/perlenv
    set -ue

    frag_size=0

    echo "${meta.method}_${meta.window_size}_${meta.report_name}"
    
    if [ "${meta.method}" == "RIPPM" ]; then

    ##echo "Making normalization for RIPPM"
    norm=\$(cat $NORM_FILE | \
        grep -E "${meta.control_samples}|${meta.treatment_samples}" | \
        awk '{ sum += \$3 } END { if (NR > 0) print sum / NR }')

    control_norm_line=\$(cat ${NORM_FILE} | \
	grep -v frag | \
	cut -f1,3 | \
        sort -t "M" -k2 -n | \
        awk -v mnorm="\$norm" '{printf "%s %.2f\\n", \$1, \$2/mnorm}'| \
        grep -E "${meta.control_samples}" | \
        cut -d " " -f2 | \
        paste -s -d " ")

    treatment_norm_line=\$(cat ${NORM_FILE} | \
	grep -v frag | \
	cut -f1,3 | \
        sort -t "M" -k2 -n | \
        awk -v mnorm="\$norm" '{printf "%s %.2f\\n", \$1, \$2/mnorm}'| \
        grep -E "${meta.treatment_samples}" | \
        cut -d " " -f2 | \
        paste -s -d " ")

    echo "treatment \${treatment_norm_line}" > norm.txt
    echo "control \${control_norm_line}" >> norm.txt
    fi

    num_ctrl_reps=\$(ls -1 *.bed | grep -E "${meta.control_samples}"  | sort -t "M" -k2 -n | wc -l)
    ctrl_reps=\$(ls -1 *.bed | grep -E "${meta.control_samples}" | sort -t "M" -k2 -n | tr "\n" " " )

    num_treatment_reps=\$(ls -1 *.bed | grep -E "${meta.treatment_samples}"  | sort -t "M" -k2 -n | wc -l)
    treatment_reps=\$(ls -1 *.bed | grep -E "${meta.treatment_samples}" | sort -t "M" -k2 -n | tr "\n" " " )

    norm_line=""
    if [ "${meta.method}" == "RIPPM" ]; then
    norm_line="--norm norm.txt"
    fi

    if [ "\${num_ctrl_reps}" == "1" ] || [ "\${num_treatment_reps}" == "1" ]; then
    (set -x; time diffReps.pl --treatment \${treatment_reps} \
        --control \${ctrl_reps} --report ${meta.report_name} \
        --gname mm9 --window ${meta.window_size} --nsd broad\
        --frag \${frag_size} --nproc $task.cpus --meth gt \${norm_line})
    else

    (set -x; time diffReps.pl --treatment \${treatment_reps} \
        --control \${ctrl_reps} --report ${meta.report_name} \
        --gname mm9 --window ${meta.window_size} --nsd broad\
        --frag \${frag_size} --nproc $task.cpus \${norm_line})
    fi
    
    """
}

process diffreps_summary {
    tag "${meta.num}_${meta.method}_${meta.window_size}"

    beforeScript 'source $HOME/.bashrc'
    
    executor 'local'
    // cpus 1
    errorStrategy 'retry'
    maxRetries 2
    publishDir path: "${params.output_dir}/diffreps_output/${meta.group_name}/", mode: "copy", pattern: "*.xlsx", overwrite: true
    
    input:
    tuple val(meta), path(diffout), path(MUMERGE)

    output:
    tuple val(meta), path("*.xlsx"), emit: full_report
    path("${meta.report_name}"), emit: diffreps_output 

    script:
    output_prefix="${meta.group_name}_MuMerge_${meta.method}_${meta.window_size}"

    """
    module load R/${params.rversion}
    diffreps_diffpeak_analysis.R --annotated_path ${meta.report_name}.annotated \
        --min_avg_count 20 \
        --log2fc_cutoff 1 \
        --mumerge_path ${MUMERGE} \
        --output_prefix ${output_prefix}

    """
}

process collect_diffreps_norm_factors {

    executor 'local'
    publishDir path: "${params.output_dir}/summary/", mode: "copy", pattern: "*.xlsx", overwrite: true
    // echo true

    input:
    tuple path(diffreps_output), path(sample_stats), path(fq_num_reads)
    output:
    path("*.xlsx")

    script:
    """
    module load R/${params.rversion}
    diffreps_output_parser.R --path "." \
        --rippm_report ${sample_stats} \
        --fastq_num_reads ${fq_num_reads}
    """
}

workflow DIFFREPS {
    take:
    diffreps_config      //parse_configuration_xls.out.diffreps_config
    fragments_bed6       //bam_count.out.fragments_bed6
    norm_factors         //calc_norm_factors.out
    fq_num_reads         //table with number of reads in R1.fq files for each sample
    mumerge_peaks              //mumerge peaks instead of MACS2 union

    main:

    diffreps_params = diffreps_config
        .combine(
            fragments_bed6
                .map{it->it[1]}
                .collect().toList())
        .map{meta,l ->
            def tl = l.findAll{it =~ meta.treatment_samples}
            def cl = l.findAll{it =~ meta.control_samples}
            return [meta, tl, cl]
        }
    
    diffreps_params_withnorm = diffreps_params.combine(norm_factors)
    diffreps_params_withnorm | diffreps


    prelim_ch = diffreps_params.map{it -> it[0..-3]} // exclude bed6 files
        .join(diffreps.out)
        .map{meta,diffout ->
            [meta.group_name,meta,diffout]} 
        .combine(mumerge_peaks, by: 0) 
        .map{g,m,d,mu -> [m,d,mu]}

    prelim_ch | diffreps_summary

    diffreps_summary.out.diffreps_output
        .collect()
        .toList()
        .combine(norm_factors) 
        .combine(fq_num_reads) | collect_diffreps_norm_factors
    
    emit:
    full_report = diffreps_summary.out.full_report
}
