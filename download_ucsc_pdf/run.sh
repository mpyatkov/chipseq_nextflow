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

# Prompt for UCSC credentials (one set for all jobs)
read -p "UCSC login: " UCSC_LOGIN
read -s -p "UCSC password: " UCSC_PASSWORD
echo
if [ -z "${UCSC_LOGIN}" ] || [ -z "${UCSC_PASSWORD}" ]; then
    echo "Error: UCSC login and password are required."
    exit 1
fi

# Defaults for DB and zoom (override by exporting UCSC_DB/UCSC_ZOOM before running)
UCSC_DB="${UCSC_DB:-mm9}"
UCSC_ZOOM="${UCSC_ZOOM:-3}"

# Export for qsub -V (do NOT echo credentials)
export UCSC_LOGIN UCSC_PASSWORD UCSC_DB UCSC_ZOOM

## Check only the last RESULTS directory if we have multiple of them
TOP25DIR=$(ls -t ../ | grep RESULTS | head -1)
TOP25DIR_PATH="../${TOP25DIR}"

for session in "$@"; do
    for top25file in $(find ${TOP25DIR_PATH} -name "*TOP25*.xlsx" | grep by_comparison | grep NOXY); do
        fname=$(basename $top25file)
        fname_noext=${fname%.*}
        top25pdf="${session}__${fname_noext}_ucsc.pdf"
        
        qsub -V -j y -o "${fname}.log" -N "UCSC_${fname}" download.qsub "${session}" "${top25file}" "${top25pdf}"
    done
    echo "Starting jobs for session: ${session}..."
    echo "Please wait until all downloads are complete..."
    wait_for_jobs "UCSC"
done

## Wait until all jobs will be complete
# echo "You can check current status of jobs using qstat -u $USER | grep UCSC"
# echo "Please wait until all downloads are complete..."
# wait_for_jobs "UCSC"

## Check if some log files does not contain IAMOK as the last line
## which means something happen during download
NUMBER_OF_BAD_LOGS=$(grep -rL IAMOK *.log | wc -l)
if [ "${NUMBER_OF_BAD_LOGS}" -ne 0 ]; then
    echo "Some logs contain errors: "
    grep -rL IAMOK *.log
    exit 1
else
    PDF_DIR="${TOP25DIR_PATH}/summary/UCSC_TOP25_PDFS" 
    mkdir -p ${PDF_DIR}
    mv *.pdf "${PDF_DIR}" 
    rm -rf *.log
    echo "All downloads are complete. Check ${PDF_DIR} for the output"
fi
