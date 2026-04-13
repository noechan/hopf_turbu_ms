#!/bin/bash
#SBATCH --job-name=MCI_ABpos_pert_SUB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=noelia.martinez@upf.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-33

ml MATLAB

matlab -nojvm -nodisplay <<EOF
pert_infocapacity_susc_stage_SUB('MCI_ABpos', ${SLURM_ARRAY_TASK_ID}, 1);
exit;
EOF