#!/bin/bash 
#SBATCH --job-name=AD_ABpos_Hopf_GEC
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=noelia.martinez@upf.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-300

#Load Matlab 2017a module

ml MATLAB

matlab -nojvm -nodisplay<<-EOF
for cond=1;
hopf_DTI_Grange_AD_ABpos_GEC(${SLURM_ARRAY_TASK_ID},cond);
end
EOF




#matlab -singleCompThread -nojvm -nodisplay < -r "run_hbif(${slurmArrayID},2,'grans');exit;" 
