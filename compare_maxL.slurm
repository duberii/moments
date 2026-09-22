#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH -A r00039


dataset="${1}"
bootstrap_number=$((${SLURM_ARRAY_TASK_ID}-1))
maxL="${2}"
all_L_to_compare="${@:3}"
echo ${all_L_to_compare}
output_directory="/N/u/rdube/Quartz/work/moments_workflow/results"

if [[ ${bootstrap_number} -eq -1 ]]; then
    dir_name="full"
else
    dir_name="bootstraps/${bootstrap_number}"
fi

mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/logs"
mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/data"
mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/accMC"
mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/genMC"

# Calculate the moments and acceptance matrix for all datasets with maximum L
root -b -q -x -e ".x /N/u/rdube/Quartz/work/moments_workflow/single_step_calculate_moments.C(\"${dataset}\",${maxL}, ${bootstrap_number})" > "${output_directory}/${dataset}/${maxL}/${dir_name}/logs/moments_and_acceptance_correcting.out" 2>&1


for i in ${all_L_to_compare}; do
    if [[ ${i} -eq ${maxL} ]]; then
        continue
    fi
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/logs"
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/data"
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/accMC"
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/genMC"
    # Get the moments and acceptance (sub)matrix from the previous jobs
    root -b -q -x -e ".x /N/u/rdube/Quartz/work/moments_workflow/single_step_calculate_moments_parasite.C(\"${dataset}\",${i}, ${bootstrap_number},${maxL})" > "${output_directory}/${dataset}/${i}/${dir_name}/logs/moments_and_acceptance_correcting.out" 2>&1

done
