#!/bin/bash

module load macs2

macs2 callpeak -t ${bam} -f ${lib} -g mm -n "${sample_id}_narrow_MACS2" --keep-dup all --nomodel --extsize 200
macs2 callpeak -t ${bam} -f ${lib} -g mm -n "${sample_id}_broad_MACS2" --keep-dup all --nomodel --extsize 200 --broad

function peak_to_bed() {
    local input=$1
    local output=$2

    awk '{print $1"\t"$2"\t"$3"\t"$4}' $input > temp1.bed

    #  Filtering peaks: Omit chrM, random, end position is greater than start
    #  start position is greater than zero, end position is greater than zero
    awk 'BEGIN {OFS="\t"} {if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0) print $0 > "temp2.bed"}' temp1.bed
    
    # echo 'Filtering peaks less than 100bp peak width from BED file'
    awk '{if($3-$2>99){print $0}}' temp2.bed > $output

    rm temp1.bed temp2.bed
}

peak_to_bed "${sample_id}_narrow_MACS2_peaks.narrowPeak" "${sample_id}_narrow_MACS2_peaks.narrowPeak.bed"
peak_to_bed "${sample_id}_broad_MACS2_peaks.broadPeak" "${sample_id}_broad_MACS2_peaks.broadPeak.bed"

