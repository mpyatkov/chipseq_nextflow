#!/bin/bash

module load macs2
module load bedtools

sample_id="!{sample_id}"
lib="!{lib}"
bam="!{bam}"
model="!{model}"
mm9_chrom_sizes="!{mm9_chrom_sizes}"

# By default MACS2 uses model to find 'shift' and 'extsize 'parameters
# if model are using then MACS2 automatically calculate 'shift' and 'extsize' params:
# in case when MACS2 uses '--nomodel' parameter, 'shift' and 'extsize' parameters the following:
# --shift=0
# --extsize=200

(
    set -x
    macs2 callpeak -t ${bam} -f ${lib} -g mm -n "${sample_id}_narrow_MACS2" --keep-dup all ${model}
    macs2 callpeak -t ${bam} -f ${lib} -g mm -n "${sample_id}_broad_MACS2" --keep-dup all ${model} --broad
)

function peak_to_bed() {
    local input=$1
    local output=$2

    awk '{print $1"\t"$2"\t"$3"\t"$4}' $input >temp1.bed

    #  Filtering peaks: Omit chrM, random, end position is greater than start
    #  start position is greater than zero, end position is greater than zero
    awk 'BEGIN {OFS="\t"} {if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0) print $0 > "temp2.bed"}' temp1.bed

    # echo 'Filtering peaks less than 100bp peak width from BED file'
    awk '{if($3-$2>99){print $0}}' temp2.bed >$output

    rm temp1.bed temp2.bed
}

peak_to_bed "${sample_id}_narrow_MACS2_peaks.narrowPeak" "${sample_id}_narrow_MACS2_peaks.narrowPeak.bed"
peak_to_bed "${sample_id}_broad_MACS2_peaks.broadPeak" "${sample_id}_broad_MACS2_peaks.broadPeak.bed"

## replace in all xls -log10(pvalue/qvalue) to minus_log10_pvalue/qvalue
sed -i -E 's/-log10\(/minus_log10_/g;s/value\)/value/g' *.xls

## create bb tracks
# ignore chrM,random and header
cat "${sample_id}_narrow_MACS2_peaks.narrowPeak.bed" | grep -vE "track|chrM|random" >tmp.bed
bedtools sort -i tmp.bed >tmp.sorted.bed
bedToBigBed -allow1bpOverlap tmp.sorted.bed ${mm9_chrom_sizes} "${sample_id}_narrow_MACS2.bb"

rm tmp*

# ignore chrM,random and header
cat "${sample_id}_broad_MACS2_peaks.broadPeak.bed" | grep -vE "track|chrM|random" >tmp.bed
bedtools sort -i tmp.bed >tmp.sorted.bed
bedToBigBed -allow1bpOverlap tmp.sorted.bed ${mm9_chrom_sizes} "${sample_id}_broad_MACS2.bb"

rm tmp*
