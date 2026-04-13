#!/bin/bash
#SBATCH --job-name=MCI_ABpos_pert_ADNI3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=noelia.martinez@upf.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-100

# Load MATLAB (adjust to your site’s module name/version)
ml MATLAB

# Be explicit about single-threaded math
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "SLURM_NODELIST=$SLURM_NODELIST"
echo "SLURM_JOB_ID=$SLURM_JOB_ID  SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"

# R2017a: no -batch; use here-doc or -r
matlab -singleCompThread -nojvm -nodisplay << EOF
for cond = 1
  pert_infocapacity_susc_MCI_ABpos_GEC(${SLURM_ARRAY_TASK_ID}, cond);
end
exit
EOF


#matlab -singleCompThread -nojvm -nodisplay < -r "run_hbif(${slurmArrayID},2,'grans');exit;" 
