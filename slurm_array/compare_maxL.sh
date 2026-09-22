#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 DATASET N_BOOTSTRAPS MAXL1 ... MAXLN"
    exit 1
fi

dataset="${1}"
n_bootstraps="${2}"
all_L_to_compare="${@:3}"
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
echo "MaxL is ${max_L_in_list}"

root -b -q -x -e ".L /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/aggregate_results.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/comparison_plots.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/maxL_comparison_plots.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/slurm_array/single_step_calculate_moments.C" || exit 1

array_job_id=$(sbatch --array=0-${n_bootstraps} --output="/dev/null" --parsable /N/u/rdube/Quartz/work/moments_workflow/slurm_array/compare_maxL.slurm "${dataset}" "${max_L_in_list}" ${all_L_to_compare})

sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${array_job_id} --output="${output_directory}/${dataset}/aggregation.out" /N/u/rdube/Quartz/work/moments_workflow/slurm_array/aggregation.slurm "${dataset}" "${n_bootstraps}" "${maxL_vector}" ${all_L_to_compare}
