#!/bin/bash

set -eu
RVERSION=4.2
xlsconfig=$1

# G222_G207_H3K27ac
dataset_label=${2:-$(pwd | xargs basename)}
config_dir="raw_configs"
output_dir="RESULTS_${dataset_label}"

echo "Copying xls config to output directory"
mkdir -p "${output_dir}/summary/" && cp ${xlsconfig} "${output_dir}/summary/"

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

# We need to preinstall packages before first run.
# If "work" directory does not exist we suppose that this is first run
if [ ! -d "work" ]; then
    # ./install_packages.sh ${RVERSION}
    # version=${1:-4.4}

    module load R/${RVERSION}
    export R_LIBS="${HOME}/R/x86_64-pc-linux-gnu-library/${RVERSION}"
    mkdir -p ${R_LIBS}

    Rscript bin/install_packages.R
fi

## split xlsx config to individual config files and put to $config_dir 
if [ -d "${config_dir}" ]; then
    module load R/${RVERSION}

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

module load java/21.0.4
module load nextflow/24.04.4

NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf --input_configs ${config_dir} --dataset_label ${dataset_label} --output_dir ${output_dir} --rversion ${RVERSION} -resume 
