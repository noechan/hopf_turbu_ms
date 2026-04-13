#!/bin/bash
#SBATCH --job-name=HC_ABpos_Hopf_GEC
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=noelia.martinez@upf.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-301

ml MATLAB

matlab -nojvm -nodisplay <<-EOF
addpath(genpath('/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/hopf_turbu_paper_project/main/code'));
hopf_DTI_Grange_stage('HC_ABpos', ${SLURM_ARRAY_TASK_ID}, 1);
EOF