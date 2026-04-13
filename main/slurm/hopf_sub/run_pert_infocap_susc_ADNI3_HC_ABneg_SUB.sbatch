#!/bin/bash
#SBATCH --job-name=HC_ABneg_pert_SUB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=noelia.martinez@upf.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-54

ml MATLAB

matlab -nojvm -nodisplay <<EOF
pert_infocapacity_susc_stage_SUB('HC_ABneg', ${SLURM_ARRAY_TASK_ID}, 1);
exit;
EOF