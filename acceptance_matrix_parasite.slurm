#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH -A r00039

time root -b -q -e ".x /N/u/rdube/Quartz/work/moments_workflow/acceptance_matrix_parasite.C+(\"${1}\",${2}, ${3}, ${4})"