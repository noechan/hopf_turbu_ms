function Ceffgroup = run_global_gec_stage(stage_label, cnd)
%RUN_GLOBAL_GEC_STAGE Compute group-level GEC for one diagnostic/Aβ stage.

    if nargin < 1 || isempty(stage_label)
        error('You must provide stage_label, e.g. ''HC_ABneg''.');
    end

    if nargin < 2 || isempty(cnd)
        cnd = 1;
    end

    clc;

    % ---------------------------------------------------------------------
    % Resolve paths relative to this file
    % ---------------------------------------------------------------------
    this_file_dir = fileparts(mfilename('fullpath'));   % .../ADNI3/scripts
    project_dir   = fileparts(this_file_dir);           % .../ADNI3
    helpers_dir   = fullfile(project_dir, 'helpers');   % .../ADNI3/helpers

    addpath(this_file_dir);

    if ~exist(helpers_dir, 'dir')
        error('Helpers folder not found: %s', helpers_dir);
    end
    addpath(helpers_dir);

    % Optional sanity check
    if exist('get_gec_paths', 'file') ~= 2
        error('get_gec_paths.m is not visible on the MATLAB path.');
    end

    % ---------------------------------------------------------------------
    % Load centralised paths and config
    % ---------------------------------------------------------------------
    paths = get_gec_paths();
    cfg   = get_gec_config();

    cfg.cnd = cnd;

    % ---------------------------------------------------------------------
    % Parse stage label -> clinical group file block
    % ---------------------------------------------------------------------
    switch stage_label
        case {'HC_ABneg', 'HC_ABpos'}
            base_group = 'HC';
        case {'MCI_ABpos'}
            base_group = 'MCI';
        case {'AD_ABpos'}
            base_group = 'AD';
        otherwise
            error('Unknown stage_label: %s', stage_label);
    end

    stage_info = paths.stage_files.(base_group);

    % ---------------------------------------------------------------------
    % Load inputs
    % ---------------------------------------------------------------------
    load(stage_info.tseries_file, 'tseries');
    load(stage_info.f_diff_file, 'f_diff');
    load(paths.sc_file, 'sc_schaefer');

    abeta_table = readtable(paths.abeta_table);

    %ptid_file_full = fullfile(paths.gec_data_dir, stage_info.ptid_file);
    load(stage_info.ptid_file, 'PTID');

    % Ensure PTID is cellstr
    PTID = cellstr(PTID);

    % ---------------------------------------------------------------------
    % Select the requested stage
    % ---------------------------------------------------------------------
    [tseries_stage, ptid_stage, meta] = select_stage_timeseries( ...
        tseries, PTID, abeta_table, stage_label);

    NSUB = size(tseries_stage, 1);
    N    = cfg.N;

    fprintf('---------------------------------------------\n');
    fprintf('Running GEC for stage: %s\n', stage_label);
    fprintf('Condition: %d\n', cnd);
    fprintf('Selected subjects: %d\n', NSUB);
    fprintf('---------------------------------------------\n');

    if NSUB == 0
        error('No subjects were selected for stage %s.', stage_label);
    end

    % ---------------------------------------------------------------------
    % Bandpass filter settings
    % ---------------------------------------------------------------------
    fnq = 1 / (2 * cfg.TR);
    Wn  = [cfg.flp / fnq, cfg.fhi / fnq];
    [bfilt, afilt] = butter(cfg.filter_order, Wn);

    % ---------------------------------------------------------------------
    % Structural connectivity initialisation
    % ---------------------------------------------------------------------
    C = sc_schaefer;
    C = C / max(C(:)) * cfg.sc_init_max;

    % ---------------------------------------------------------------------
    % Compute empirical FC and lagged covariance targets
    % ---------------------------------------------------------------------
    FC     = nan(NSUB, N, N);
    COVtau = nan(NSUB, N, N);

    for nsub = 1:NSUB
        ts = tseries_stage{nsub, cnd};

        % Defensive checks
        if isempty(ts)
            error('Empty time series for subject %d in stage %s.', nsub, stage_label);
        end

        if size(ts, 1) ~= N
            error(['Time series for subject %d has %d parcels, but cfg.N=%d. ' ...
                   'Please check the parcellation consistency.'], ...
                  nsub, size(ts,1), N);
        end

        % -------------------------------------------------------------
        % Preprocess signal
        % -------------------------------------------------------------
        signal_filt = zeros(size(ts));

        % Clean once before seed loop
        ts(isnan(ts)) = 0;
        ts(isinf(ts)) = 0;

        for seed = 1:N
            x = ts(seed, :);
            x = detrend(x - nanmean(x));
            signal_filt(seed, :) = filtfilt(bfilt, afilt, x);
        end

        ts2 = signal_filt;
        tst = ts2';

        % -------------------------------------------------------------
        % Empirical FC
        % -------------------------------------------------------------
        FCemp_sub = corrcoef(ts2');

        % -------------------------------------------------------------
        % Empirical covariance and lagged covariance
        % -------------------------------------------------------------
        COVemp_sub    = cov(ts2');
        COVtauemp_sub = zeros(N, N);
        sigratio      = zeros(N, N);

        for i = 1:N
            for j = 1:N
                sigratio(i, j) = 1 / sqrt(COVemp_sub(i, i)) / sqrt(COVemp_sub(j, j));

                [clag, lags] = xcov(tst(:, i), tst(:, j), cfg.Tau);
                indx = find(lags == cfg.Tau, 1, 'first');

                COVtauemp_sub(i, j) = clag(indx) / size(tst, 1);
            end
        end

        COVtauemp_sub = COVtauemp_sub .* sigratio;

        FC(nsub, :, :)     = FCemp_sub;
        COVtau(nsub, :, :) = COVtauemp_sub;

        fprintf('Processed subject %d / %d\n', nsub, NSUB);
    end

    FCemp     = squeeze(nanmean(FC, 1));
    COVtauemp = squeeze(nanmean(COVtau, 1));

    % ---------------------------------------------------------------------
    % Fit GEC using helper
    % ---------------------------------------------------------------------
    fit = fit_gec_from_targets(C, FCemp, COVtauemp, f_diff(:, :), cfg);

    Ceffgroup = fit.Ceffgroup;

    % ---------------------------------------------------------------------
    % Save output using original naming convention
    % ---------------------------------------------------------------------
    if ~exist(paths.output_dir, 'dir')
        mkdir(paths.output_dir);
    end

    out_name = sprintf( ...
        'ADNI3_%s_GECglo_cnd_%03d_sch1000_mod_eps_1_2_ratio_SC_max02.mat', ...
        stage_label, cnd);

    out_file = fullfile(paths.output_dir, out_name);

    save(out_file, ...
        'Ceffgroup', ...
        'FCemp', ...
        'COVtauemp', ...
        'fit', ...
        'stage_label', ...
        'cnd', ...
        'ptid_stage', ...
        'meta', ...
        '-v7.3');

    fprintf('Saved: %s\n', out_file);
end