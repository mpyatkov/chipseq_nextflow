#!/bin/bash

set -eu
xlsconfig=$1

# G222_G207_H3K27ac
dataset_label=${2:-$(pwd | xargs basename)}
config_dir="raw_configs"

RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
NORMAL=$(tput sgr0)

printf ' %s%s%s\n' "$NORMAL" "XLS config file:" "$GREEN" "${xlsconfig}"
printf ' %s%s%s\n' "$NORMAL" "Directory where raw config files will be located" "$GREEN" "${config_dir}"
printf ' %s%s%s\n' "$NORMAL" "Dataset label (will be a directory on remote server)" "$GREEN" "${dataset_label}" 
printf '%s\n' "$NORMAL"

ready=0
while true; do
    read -p "Are you ready to start pipeline and all params ok? [yN]. Press Enter to Y. " yn
    case $yn in
        [Yy]* ) ready=1; break;;
        [Nn]* ) break;;
	* ) ready=1;break;;
    esac
done

if [[ ${ready} == 0 ]]; then
    echo "Ok, check params first"
    exit 0
fi

module load nextflow

## split xlsx config to individual config files and put to $config_dir 
module load R
#set -x
if [ -d "${config_dir}" ]; then
    rm -rf _tmp && mkdir -p _tmp
    echo "Extracting config files from xls to temporary dir..."
    Rscript ./bin/config_parser.R --output_dir _tmp --input_xlsx ${xlsconfig}
    for f in $(find _tmp -name "*" -type f); do
        # echo "Processing file: $f"
        fname=$(basename $f)
        #DIFF="$(diff -q $f ${config_dir}/$fname)"
        if ! cmp -s "$f" "${config_dir}/${fname}" ; then
        #if [[ -z "$DIFF" ]]; then
            echo "Updating: $f"
            cp "$f" "${config_dir}"
        fi
    done
    rm -rf _tmp
else
    echo "Extracting config files from xls..."
    mkdir -p ${config_dir}
    Rscript ./bin/config_parser.R --output_dir ${config_dir} --input_xlsx ${xlsconfig}    
fi
module unload R

NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf --input_configs ${config_dir} --dataset_label ${dataset_label} -resume

