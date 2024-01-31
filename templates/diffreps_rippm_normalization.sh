#!/bin/bash

samples_rippm_norm=!{samples_rippm_norm}
CONTROL_SAMPLES="!{CONTROL_SAMPLES}"
TREATMENT_SAMPLES="!{TREATMENT_SAMPLES}"

norm=$(cat $samples_rippm_norm | \
           grep -E "${CONTROL_SAMPLES}|${TREATMENT_SAMPLES}" | \
           awk '{ sum += $3 } END { if (NR > 0) print sum / NR }')

control_norm_line=$(cat ${samples_rippm_norm} | \
	                grep -v frag | \
	                cut -f1,3 | \
                        sort -t "M" -k2 -n |\
                        awk -v mnorm="$norm" '{printf "%s %.2f\n", $1, $2/mnorm}'| \
                        grep -E "${CONTROL_SAMPLES}" | \
                        cut -d " " -f2 | \
                        paste -s -d " ")

treatment_norm_line=$(cat ${samples_rippm_norm} | \
	                  grep -v frag | \
	                  cut -f1,3 | \
                          sort -t "M" -k2 -n |\
                          awk -v mnorm="$norm" '{printf "%s %.2f\n", $1, $2/mnorm}'| \
                          grep -E "${TREATMENT_SAMPLES}" | \
                          cut -d " " -f2 | \
                          paste -s -d " ")

echo "treatment ${treatment_norm_line}" > norm.txt
echo "control ${control_norm_line}" >> norm.txt
