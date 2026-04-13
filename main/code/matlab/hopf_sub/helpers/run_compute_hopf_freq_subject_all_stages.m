function run_compute_hopf_freq_subject_all_stages()
% RUN_COMPUTE_HOPF_FREQ_SUBJECT_ALL_STAGES
%
% Runs subject-level Hopf frequency computation for all stages.

    stage_keys = { ...
        'HC_ABneg', ...
        'HC_ABpos', ...
        'MCI_ABpos', ...
        'AD_ABpos' ...
    };

    for i = 1:numel(stage_keys)
        fprintf('\n\n########################################\n');
        fprintf('Running stage %s (%d/%d)\n', stage_keys{i}, i, numel(stage_keys));
        fprintf('########################################\n');

        compute_hopf_freq_subject_stage(stage_keys{i});
    end
end