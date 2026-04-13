clear
clc

addpath(genpath('/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/hopf_turbu_paper_project/main/code/matlab/hopf_sub/'))
stages = {'HC_ABneg', 'HC_ABpos', 'MCI_ABpos', 'AD_ABpos'};

for i = 1:numel(stages)
    fprintf('\n=============================\n');
    fprintf('Processing stage: %s\n', stages{i});
    fprintf('=============================\n');

    compute_hopf_freq_subject_stage(stages{i});
end