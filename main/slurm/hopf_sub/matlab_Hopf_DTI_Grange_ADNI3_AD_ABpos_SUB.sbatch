#!/bin/bash
#SBATCH --job-name=AD_ABpos_Hopf_SUB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=noelia.martinez@upf.edu
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-26 #NSUB AD AB neg

# Load Matlab module
ml MATLAB

matlab -nojvm -nodisplay << EOF
    nsub = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    cond = 1;
    G_range = 0:0.01:0.4;
    NG = numel(G_range);

    fprintf('Starting subject %d, NG = %d\n', nsub, NG);

    for s = 1:NG
        fprintf('Subject %d: s = %d / %d\n', nsub, s, NG);
        hopf_DTI_Grange_AD_ABpos_GEC_SUB(nsub, s, cond);
    end
    fprintf('Finished subject %d\n', nsub);
    exit
EOF
