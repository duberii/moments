#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 DATASET N_BOOTSTRAPS MAXL1 ... MAXLN"
    exit 1
fi

dataset="${1}"
n_bootstraps="${2}"
all_L_to_compare="${@:3}"
echo $all_L_to_compare
max_L_in_list="${3}"
maxL_vector="{"
counter=0
for i in ${all_L_to_compare}; do
    if [[ ${counter} -eq 0 ]]; then
        maxL_vector="${maxL_vector}${i}"
    else 
        maxL_vector="${maxL_vector},${i}"
    fi
    if [[ ${i} -gt ${max_L_in_list} ]]; then
        max_L_in_list=${i}
    fi
    counter+=1
done
maxL_vector="${maxL_vector}}"
output_directory="/N/u/rdube/Quartz/work/moments_workflow/results"

echo "Comparing maxL=${all_L_to_compare}"

root -b -q -x -e ".L /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/aggregate_results.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/comparison_plots.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/maxL_comparison_plots.C+"|| exit 1

maxL_job_id=$(/N/u/rdube/Quartz/work/moments_workflow/moments.sh "${dataset}" "${n_bootstraps}" "${max_L_in_list}")
jobs_to_wait_for_to_compare="${maxL_job_id}"

for i in ${all_L_to_compare}; do
    if [[ ${i} -eq ${max_L_in_list} ]]; then
        continue
    fi
    job_id=$(/N/u/rdube/Quartz/work/moments_workflow/moments_parasite.sh "${dataset}" "${n_bootstraps}" "${i}" "${max_L_in_list}" "${maxL_job_id}")
    jobs_to_wait_for_to_compare="${jobs_to_wait_for_to_compare}:${job_id}"
done
sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${jobs_to_wait_for_to_compare} --output="${output_directory}/${dataset}/maxL_compare.out" /N/u/rdube/Quartz/work/moments_workflow/maxL_comparison_plots.slurm "${dataset}" "${maxL_vector}"
