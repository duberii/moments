#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=05:00:00
#SBATCH -A r00039

dataset_name=${1}
maxL=${2}
bootstrap_number=${3}

mkdir -p "/N/u/rdube/Quartz/work/moments_workflow/results/${dataset_name}/${maxL}/full/logs"
mkdir -p "/N/u/rdube/Quartz/work/moments_workflow/results/${dataset_name}/${maxL}/full/data"
mkdir -p "/N/u/rdube/Quartz/work/moments_workflow/results/${dataset_name}/${maxL}/full/accMC"
mkdir -p "/N/u/rdube/Quartz/work/moments_workflow/results/${dataset_name}/${maxL}/full/genMC"

root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments.C+(\"${dataset_name}\",\"genMC\",${maxL}, ${bootstrap_number})"
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments.C+(\"${dataset_name}\",\"accMC\",${maxL}, ${bootstrap_number})"
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments.C+(\"${dataset_name}\",\"data\",${maxL}, ${bootstrap_number})"

root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix.C+(\"${dataset_name}\",${maxL},${bootstrap_number})"

root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.C+(\"${dataset_name}\",\"data\",${maxL},${bootstrap_number})"
root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.C+(\"${dataset_name}\",\"accMC\",${maxL},${bootstrap_number})"

