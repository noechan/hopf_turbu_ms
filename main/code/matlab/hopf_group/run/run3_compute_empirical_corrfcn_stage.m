clear
clc

stages = {'HC_ABneg', 'HC_ABpos', 'MCI_ABpos', 'AD_ABpos'};

for i = 1:numel(stages)
    compute_empirical_corrfcn_stage(stages{i});
end