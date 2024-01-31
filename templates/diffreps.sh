#!/bin/bash

NUM="!{NUM}"
TREATMENT_NAME="!{TREATMENT_NAME}"
CONTROL_NAME="!{CONTROL_NAME}"
TREATMENT_SAMPLES="!{TREATMENT_SAMPLES}"
CONTROL_SAMPLES="!{CONTROL_SAMPLES}"
NORMALIZATION="!{NORMALIZATION}"
WINDOW_SIZE="!{WINDOW_SIZE}"
TREATMENT_FILES="!{TREATMENT_FILES}"
CONTROL_FILES="!{CONTROL_FILES}"
NORM_FILE="!{NORM_FILE}"
report_name="!{report_name}"
NSLOT="!{task.cpu}"
FRAG_SIZE=0

module load bedtools
module load miniconda

echo "${report_name}"
## load conda environment for diffReps
echo "Activating conda environment"
conda activate /projectnb/wax-es/routines/condaenv/perlenv

if [ "${NORMALIZATION}" == "RIPPM" ]; then

    echo "Making normalization for RIPPM"
    
    norm=$(cat $NORM_FILE | \
               grep -E "${CONTROL_SAMPLES}|${TREATMENT_SAMPLES}" | \
               awk '{ sum += $3 } END { if (NR > 0) print sum / NR }')

    control_norm_line=$(cat ${NORM_FILE} | \
	                    grep -v frag | \
	                    cut -f1,3 | \
                            sort -t "M" -k2 -n |\
                            awk -v mnorm="$norm" '{printf "%s %.2f\n", $1, $2/mnorm}'| \
                            grep -E "${CONTROL_SAMPLES}" | \
                            cut -d " " -f2 | \
                            paste -s -d " ")

    treatment_norm_line=$(cat ${NORM_FILE} | \
	                      grep -v frag | \
	                      cut -f1,3 | \
                              sort -t "M" -k2 -n |\
                              awk -v mnorm="$norm" '{printf "%s %.2f\n", $1, $2/mnorm}'| \
                              grep -E "${TREATMENT_SAMPLES}" | \
                              cut -d " " -f2 | \
                              paste -s -d " ")

    echo "treatment ${treatment_norm_line}" > norm.txt
    echo "control ${control_norm_line}" >> norm.txt
fi


num_ctrl_reps=$(ls -1 *.bed | grep -E ${CONTROL_SAMPLES}  | sort -t "M" -k2 -n | wc -l)
ctrl_reps=$(ls -1 *.bed | grep -E ${CONTROL_SAMPLES} | sort -t "M" -k2 -n | tr "\n" " " )

num_treatment_reps=$(ls -1 *.bed | grep -E ${TREATMENT_SAMPLES}  | sort -t "M" -k2 -n | wc -l)
treatment_reps=$(ls -1 *.bed | grep -E ${TREATMENT_SAMPLES} | sort -t "M" -k2 -n | tr "\n" " " )


# report_name="diffReps_${TREATMENT_NAME}.vs.${CONTROL_NAME}"

norm_line=""
if [ "${NORMALIZATION}" == "RIPPM" ]; then
    norm_line="--norm norm.txt"
fi

if [ "${num_ctrl_reps}" == "1" ] || [ "${num_treatment_reps}" == "1" ]; then
    (set -x; time diffReps.pl --treatment ${treatment_reps} \
                  --control ${ctrl_reps} --report ${report_name} \
                  --gname mm9 --window ${WINDOW_SIZE} --nsd broad\
                  --frag ${FRAG_SIZE} --nproc $NSLOTS --meth gt ${norm_line})
else
    echo $NORMALIZATION
    (set -x; time diffReps.pl --treatment ${treatment_reps} \
                  --control ${ctrl_reps} --report ${report_name} \
                  --gname mm9 --window ${WINDOW_SIZE} --nsd broad\
                  --frag ${FRAG_SIZE} --nproc $NSLOTS ${norm_line})
fi

# ls -l
# diffReps_Hnf6_Male.vs.Hnf6_Female
# diffReps_Hnf6_Male.vs.Hnf6_Female.annotated
# diffReps_Hnf6_Male.vs.Hnf6_Female.hotspot
