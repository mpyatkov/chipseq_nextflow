#!/bin/bash

NUM="!{NUM}"
TREATMENT_NAME="!{TREATMENT_NAME}"
CONTROL_NAME="!{CONTROL_NAME}"
TREATMENT_SAMPLES="!{TREATMENT_SAMPLES}"
CONTROL_SAMPLES="!{CONTROL_SAMPLES}"
NORMALIZATION="!{NORMALIZATION}"
WINDOW_SIZE="!{WINDOW_SIZE}"
report_name="!{report_name}"
output_dir="!{output_dir}"

PEAKCALLER="!{peakcaller}"
DATASET_LABEL="!{dataset_label}"
SAMPLE_LABELS="!{SAMPLE_LABELS}"

mm9_chrom_sizes="!{mm9_chrom_sizes}"

mkdir -p XLSfiles && mv *.xls ./XLSfiles

module load R
module load bedtools

(set -x; diffreps_genomicRanges.R \
                 --annotated_path ${report_name}.annotated \
                 --hotspot_path ${report_name}.hotspot \
                 --macs2_xls_dir_path ./XLSfiles/ \
                 --sample_labels_path ${SAMPLE_LABELS} \
                 --min_avg_count 20 \
                 --log2fc_cutoff 1 \
                 --control_name ${CONTROL_NAME} \
                 --treatment_name ${TREATMENT_NAME} \
                 --peak_caller ${PEAKCALLER} \
                 --histone_mark ${DATASET_LABEL} \
                 --normalization_caller "${NORMALIZATION}_${WINDOW_SIZE}" \
                 --treatment_samples ${TREATMENT_SAMPLES}  \
                 --control_samples ${CONTROL_SAMPLES})

mkdir -p plots && cp *.pdf ./plots
mkdir -p fullreport && cp *.xlsx ./fullreport
mkdir -p ucsc_tracks && cp *.bed ./ucsc_tracks

# ignore chrM,random and header
for bed in `find . -name "*_FILTERED*.bed"`; do
    cat ${bed} | grep -vE "track|chrM|random" > tmp.bed
    bedtools sort -i tmp.bed > tmp.sorted.bed
    bedToBigBed -allow1bpOverlap tmp.sorted.bed ${mm9_chrom_sizes} "${output_dir}.bb"
    rm tmp.bed
done

mkdir -p ./${output_dir}/
mv plots ./${output_dir}/
mv fullreport ./${output_dir}/
mv ucsc_tracks ./${output_dir}/
cp ${report_name}* ./${output_dir}/

rm -rf XLSfiles
rm -rf ${SAMPLE_LABELS}



