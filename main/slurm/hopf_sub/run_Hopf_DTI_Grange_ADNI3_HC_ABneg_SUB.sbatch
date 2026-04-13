#!/bin/bash
#SBATCH --job-name=HC_ABneg_Hopf_SUB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=noelia.martinez@upf.edu
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-54

set -euo pipefail

ml MATLAB

matlab -nojvm -nodisplay <<-EOF
try
    addpath(genpath('/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/hopf_turbu_paper_project/main/code'));
    hopf_DTI_Grange_stage_SUB('HC_ABneg', ${SLURM_ARRAY_TASK_ID}, 1);
catch ME
    fprintf(2, '%s\n', getReport(ME, 'extended'));
    exit(1);
end
exit
EOF