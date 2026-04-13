clear all; close all; clc;

paths = get_paths();

hc_file  = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/participants/cross_sectional/ADNI3_HC_ALL_MRI_sc_analysis_N238rev_idx4MRIselection_v2.xlsx';
mci_file = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/participants/cross_sectional/ADNI3_MCI_ALL_MRI_sc_analysis_N238rev_idx4MRIselection_v2.xlsx';
ad_file  = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/participants/cross_sectional/ADNI3_AD_ALL_MRI_sc_analysis_N238rev_idx4MRIselection_v2.xlsx';

demo_hc  = get_demographics_table(hc_file);
demo_mci = get_demographics_table(mci_file);
demo_ad  = get_demographics_table(ad_file);

demo_tbl = [demo_hc; demo_mci; demo_ad];

stage_keys = {'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'};
stage_info = struct();

for i = 1:numel(stage_keys)

    stage_key = stage_keys{i};
    fprintf('\nProcessing stage: %s\n', stage_key);

    cfg = get_stage_config(stage_key);
    [stage_mask, PTID] = get_stage_mask(cfg, paths);

    stage_ptids = PTID(stage_mask);

    [age_all, sex_all, group_all] = get_covariates_for_ptids(stage_ptids, demo_tbl);

    stage_info.(stage_key).PTID_all  = stage_ptids;
    stage_info.(stage_key).age_all   = age_all;
    stage_info.(stage_key).sex_all   = sex_all;
    stage_info.(stage_key).group_all = group_all;

    fprintf('  n subjects = %d\n', numel(stage_ptids));
end

fprintf('\nAll stages processed successfully.\n');