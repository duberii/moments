#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 DATASET N_BOOTSTRAPS SMALLER_L BIGGER_L JOB_ID_TO_WAIT_FOR"
    exit 1
fi

dataset="${1}"
n_bootstraps="${2}"
maxL="${3}"
largerMaxL="${4}"
largerMaxL_job_id=${5}
output_directory="/N/u/rdube/Quartz/work/moments_workflow/results"

echo "Calculating acceptance corrected moments for ${dataset} with ${n_bootstraps} bootstrapping instances. maxL=${maxL}. This job will be a parasite for maxL=${largerMaxL}, which has (final) pid= ${largerMaxL_job_id}." >&2
mkdir -p "${output_directory}/${dataset}/${maxL}/full/logs"
mkdir -p "${output_directory}/${dataset}/${maxL}/full/data"
mkdir -p "${output_directory}/${dataset}/${maxL}/full/accMC"
mkdir -p "${output_directory}/${dataset}/${maxL}/full/genMC"

root -b -q -x -e ".L /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/aggregate_results.C+" -e ".L /N/u/rdube/Quartz/work/moments_workflow/comparison_plots.C+" &> "${output_directory}/${dataset}/${maxL}/compilation.out"|| exit 1

#submit unsampled data job
data_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/full/logs/data.out" /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.slurm ${dataset} "data" ${maxL} "-1" ${largerMaxL})
#submit unsampled accMC job
accMC_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/full/logs/accMC.out" /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.slurm ${dataset} "accMC" ${maxL} "-1" ${largerMaxL})
#submit unsampled genMC job
genMC_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/full/logs/genMC.out" /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.slurm ${dataset} "genMC" ${maxL} "-1" ${largerMaxL})
#submit unsampled acceptance matrix job
acceptance_matrix_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/full/logs/acceptance_matrix.out" /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.slurm ${dataset} ${maxL} "-1" ${largerMaxL})

#submit unsampled acceptance correcting job
jobs_to_wait_for_acceptance_correcting="${largerMaxL_job_id}:${data_job_id}:${accMC_job_id}:${genMC_job_id}:${acceptance_matrix_job_id}"
acceptance_correct_data_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${jobs_to_wait_for_acceptance_correcting} --output="${output_directory}/${dataset}/${maxL}/full/logs/data_acceptance_correcting.out" /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.slurm ${dataset} "data" ${maxL} "-1")
acceptance_correct_accMC_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${jobs_to_wait_for_acceptance_correcting} --output="${output_directory}/${dataset}/${maxL}/full/logs/accMC_acceptance_correcting.out" /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.slurm ${dataset} "accMC" ${maxL} "-1")

jobs_to_wait_for="${acceptance_correct_data_job_id}:${acceptance_correct_accMC_job_id}"
for i in $(seq 0 $((${n_bootstraps}-1))); do
    mkdir -p "${output_directory}/${dataset}/${maxL}/bootstrap_${i}/logs"
    mkdir -p "${output_directory}/${dataset}/${maxL}/bootstrap_${i}/data"
    mkdir -p "${output_directory}/${dataset}/${maxL}/bootstrap_${i}/accMC"
    mkdir -p "${output_directory}/${dataset}/${maxL}/bootstrap_${i}/genMC"
    #submit unsampled data job
    data_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/bootstrap_${i}/logs/data.out" /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.slurm ${dataset} "data" ${maxL} ${i} ${largerMaxL})
    #submit unsampled accMC job
    accMC_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/bootstrap_${i}/logs/accMC.out" /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.slurm ${dataset} "accMC" ${maxL} ${i} ${largerMaxL})
    #submit unsampled genMC job
    genMC_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/bootstrap_${i}/logs/genMC.out" /N/u/rdube/Quartz/work/moments_workflow/calculate_measured_moments_parasite.slurm ${dataset} "genMC" ${maxL} ${i} ${largerMaxL})
    #submit unsampled acceptance matrix job
    acceptance_matrix_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${largerMaxL_job_id} --output="${output_directory}/${dataset}/${maxL}/bootstrap_${i}/logs/acceptance_matrix.out" /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.slurm ${dataset} ${maxL} ${i} ${largerMaxL})
    #submit unsampled acceptance correcting job
    jobs_to_wait_for_acceptance_correcting="${largerMaxL_job_id}:${data_job_id}:${accMC_job_id}:${genMC_job_id}:${acceptance_matrix_job_id}"
    acceptance_correct_data_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${jobs_to_wait_for_acceptance_correcting} --output="${output_directory}/${dataset}/${maxL}/bootstrap_${i}/logs/data_acceptance_correcting.out" /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.slurm ${dataset} "data" ${maxL} ${i})
    acceptance_correct_accMC_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${jobs_to_wait_for_acceptance_correcting} --output="${output_directory}/${dataset}/${maxL}/bootstrap_${i}/logs/accMC_acceptance_correcting.out" /N/u/rdube/Quartz/work/moments_workflow/acceptance_correct.slurm ${dataset} "accMC" ${maxL} ${i})
    jobs_to_wait_for="${jobs_to_wait_for}:${acceptance_correct_data_job_id}:${acceptance_correct_accMC_job_id}"
done
if [ ${n_bootstraps} -gt 0 ]; then
    aggregate_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${jobs_to_wait_for} --output="${output_directory}/${dataset}/${maxL}/full/logs/aggregation.out" /N/u/rdube/Quartz/work/moments_workflow/aggregate_results.slurm "${dataset}" ${maxL} ${n_bootstraps})
    jobs_to_wait_for="${jobs_to_wait_for}:${aggregate_job_id}"
fi
sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${jobs_to_wait_for} --output="${output_directory}/${dataset}/${maxL}/full/logs/comparison_plots.out" /N/u/rdube/Quartz/work/moments_workflow/comparison_plots.slurm "${dataset}" ${maxL}
