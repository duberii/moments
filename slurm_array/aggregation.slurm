#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -A r00039

dataset="${1}"
n_bootstraps="${2}"
maxL_vector="${3}"
all_L_to_compare="${@:4}"
output_directory="/N/u/rdube/Quartz/work/moments_workflow/results"

for i in ${all_L_to_compare}; do
    if [[ ${n_bootstraps} -gt 0 ]]; then
        root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/aggregate_results.C+(\"${dataset}\",${i},${n_bootstraps})" > "${output_directory}/${dataset}/${i}/full/logs/aggregation.out" 2>&1
    fi
    root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/comparison_plots.C+(\"${dataset}\", ${i})" > "${output_directory}/${dataset}/${i}/full/logs/comparison_plots.out" 2>&1
done

root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/maxL_comparison_plots.C+(\"${dataset}\",std::vector<int>(${maxL_vector}))" > "${output_directory}/${dataset}/maxL_compare.out" 2>&1
