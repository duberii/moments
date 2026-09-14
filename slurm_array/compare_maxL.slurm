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
    dir_name="bootstrap_${bootstrap_number}"
fi

mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/logs"
mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/data"
mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/accMC"
mkdir -p "${output_directory}/${dataset}/${maxL}/${dir_name}/genMC"

# Calculate the moments and acceptance matrix for all datasets with maximum L
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments.C+(\"${dataset}\",\"data\",${maxL}, ${bootstrap_number})" > "${output_directory}/${dataset}/${maxL}/${dir_name}/logs/moments_data.out" 2>&1
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments.C+(\"${dataset}\",\"accMC\",${maxL}, ${bootstrap_number})" > "${output_directory}/${dataset}/${maxL}/${dir_name}/logs/moments_accMC.out" 2>&1
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments.C+(\"${dataset}\",\"genMC\",${maxL}, ${bootstrap_number})" > "${output_directory}/${dataset}/${maxL}/${dir_name}/logs/moments_genMC.out" 2>&1
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix.C+(\"${dataset}\",${maxL},${bootstrap_number})" > "${output_directory}/${dataset}/${maxL}/${dir_name}/logs/acceptance_matrix.out" 2>&1

#acceptance correct accMC and data with maximum L
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.C+(\"${dataset}\",\"data\",${maxL},${bootstrap_number})" > "${output_directory}/${dataset}/${maxL}/${dir_name}/logs/acceptance_correct_data.out" 2>&1
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.C+(\"${dataset}\",\"accMC\",${maxL},${bootstrap_number})" > "${output_directory}/${dataset}/${maxL}/${dir_name}/logs/acceptance_correct_accMC.out" 2>&1


for i in ${all_L_to_compare}; do
    if [[ ${i} -eq ${max_L_in_list} ]]; then
        continue
    fi
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/logs"
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/data"
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/accMC"
    mkdir -p "${output_directory}/${dataset}/${i}/${dir_name}/genMC"
    # Get the moments and acceptance (sub)matrix from the previous jobs
    root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.C+(\"${dataset}\",\"data\",${i}, ${bootstrap_number},${maxL})" > "${output_directory}/${dataset}/${i}/${dir_name}/logs/moments_data.out" 2>&1
    root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.C+(\"${dataset}\",\"accMC\",${i}, ${bootstrap_number},${maxL})" > "${output_directory}/${dataset}/${i}/${dir_name}/logs/moments_accMC.out" 2>&1
    root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.C+(\"${dataset}\",\"genMC\",${i}, ${bootstrap_number},${maxL})" > "${output_directory}/${dataset}/${i}/${dir_name}/logs/moments_genMC.out" 2>&1
    root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.C+(\"${dataset}\",${i},${bootstrap_number},${maxL})" > "${output_directory}/${dataset}/${i}/${dir_name}/logs/acceptance_matrix.out" 2>&1

    #acceptance correct accMC and data using the acceptance (sub)matrix from previous step
    root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.C+(\"${dataset}\", \"data\", ${i}, ${bootstrap_number})" > "${output_directory}/${dataset}/${i}/${dir_name}/logs/acceptance_correct_data.out" 2>&1
    root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.C+(\"${dataset}\", \"accMC\", ${i}, ${bootstrap_number})" > "${output_directory}/${dataset}/${i}/${dir_name}/logs/acceptance_correct_accMC.out" 2>&1
done
