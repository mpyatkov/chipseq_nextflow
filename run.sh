#!/bin/bash
module load nextflow
# config=$1 # --xlsx_config ${config}
NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf  -resume
#NXF_OPTS='-Xms500M -Xmx2G' nextflow run main.nf -resume -stub-run
