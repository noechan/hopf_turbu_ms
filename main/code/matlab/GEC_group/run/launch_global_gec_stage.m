clear;
clc;

% ------------------------------------------------------------
% Configuration
% ------------------------------------------------------------
cnd = 1;

stage_labels = { ...
    %'HC_ABneg', ...
    'HC_ABpos', ...
    'MCI_ABpos', ...
    'AD_ABpos' ...
    };

% ------------------------------------------------------------
% Run stages
% ------------------------------------------------------------
results = struct( ...
    'stage_label', {}, ...
    'success', {}, ...
    'error_message', {} ...
    );

fprintf('=============================================\n');
fprintf('Launching global GEC stage runs\n');
fprintf('Condition: %d\n', cnd);
fprintf('=============================================\n');

for k = 1:numel(stage_labels)
    stage_label = stage_labels{k};

    fprintf('\n---------------------------------------------\n');
    fprintf('Running stage %d/%d: %s\n', k, numel(stage_labels), stage_label);
    fprintf('---------------------------------------------\n');

    results(k).stage_label = stage_label;
    results(k).success = false;
    results(k).error_message = '';

    try
        run_global_gec_stage(stage_label, cnd);
        results(k).success = true;
        fprintf('Finished successfully: %s\n', stage_label);

    catch ME
        results(k).success = false;
        results(k).error_message = ME.message;
        fprintf('ERROR in %s\n', stage_label);
        fprintf('%s\n', ME.message);
    end
end

% ------------------------------------------------------------
% Summary
% ------------------------------------------------------------
fprintf('\n=============================================\n');
fprintf('Run summary\n');
fprintf('=============================================\n');

for k = 1:numel(results)
    if results(k).success
        fprintf('[OK]   %s\n', results(k).stage_label);
    else
        fprintf('[FAIL] %s --> %s\n', results(k).stage_label, results(k).error_message);
    end
end