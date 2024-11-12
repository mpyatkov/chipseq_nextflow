#!/bin/bash
set -eu
## USAGE
## ./run.sh session1 [session2] [session3] ...
## support multiple sessions separated by space
## if session name contains spaces enclose in quotes:  "some session name"

### cleanup after Ctrl-C
function cleanup() {
    echo "Caught Ctrl+C! Running cleanup..."
    JOB_IDS=$(qstat -u $USER | grep "UCSC" | awk '{print $1}')
    if [ -z "${JOB_IDS}" ]; then
        echo "No jobs found"
    else
        for JOB_ID in $JOB_IDS; do
            qdel "$JOB_ID"
        done
    fi
}

## catch the Ctrl-C here and run cleanup
trap cleanup SIGINT
### end cleanup

## wait jobs with specific prefix
function wait_for_jobs() {
    local prefix="$1"
    local job_ids=$(qstat | grep "$prefix" | awk '{print $1}')

    while [[ -n "$job_ids" ]]; do
        sleep 60 ## wait 1 minute
        job_ids=$(qstat | grep "$prefix" | awk '{print $1}')
    done
}

### check if any sessions provided
if [ $# -eq 0 ]; then
    echo "Error: you did not provide any UCSC session as argument..."
    echo "./run.sh session1 [session2] [session3] ..."
    exit 1
fi

## Take the R version from pipelines run.sh script
## I supposed that you already started pipeline and pipeline have already installed
## all the required packages
RVERSION=$(cat ../run.sh | grep RVERSION= | cut -d "=" -f2)

## Check only the last RESULTS directory if we have multiple of them
TOP25DIR=$(ls -t ../ | grep RESULTS | head -1)
TOP25DIR_PATH="../${TOP25DIR}"

module load "R/${RVERSION}"

for session in "$@"; do
    for top25file in $(find ${TOP25DIR_PATH} -name "*top25*.xlsx"); do
        # echo "processing $top25file"
        fname=$(basename $top25file)
        qsub -j y -o "${fname}.log" -N "UCSC_${fname}" download.qsub ${session} ${top25file} ${RVERSION}
    done
done

## Wait until all jobs will be complete
echo "You can check current status of jobs using qstat -u $USER | grep UCSC"
echo "Please wait until all downloads are complete..."
wait_for_jobs "UCSC"

## Check if some log files does not contain IAMOK as the last line
## which means something happen during download
NUMBER_OF_BAD_LOGS=$(grep -rL IAMOK *.log | wc -l)
if [ "${NUMBER_OF_BAD_LOGS}" -ne 0 ]; then
    echo "Some logs contain errors: "
    grep -rL IAMOK *.log
    exit 1
else
    echo "All downloads are complete"
    rm -rf *.log
fi
