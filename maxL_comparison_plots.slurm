#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH -A r00039

echo "Using $SLURM_CPUS_PER_TASK CPUs"
time root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/maxL_comparison_plots.C+(\"${1}\",std::vector<int>(${2}))"
