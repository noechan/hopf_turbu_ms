function compute_hopf_freq_stage(stage_key)
% COMPUTE_HOPF_FREQ_STAGE
%
% Stage-driven version of the original frequency / empirical FC script.
% Processes one analysis stage at a time:
%   - HC_ABneg
%   - HC_ABpos
%   - MCI_ABpos
%   - AD_ABpos
%
% For each stage, the function:
%   - loads the parent-group tseries
%   - selects only subjects belonging to that stage
%   - computes empirical FC from filtered time series
%   - computes dominant frequency per parcel
%
% Outputs are saved to:
%   <project_root>/results/hopf_group/freq/<STAGE>/

    paths = get_paths();
    cfg   = get_stage_config(stage_key);

    addpath(genpath(paths.code_root));

    % -----------------------------------------
    % Create results directory
    % -----------------------------------------
    results_dir = fullfile(paths.freq_root, cfg.results_subdir);
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fprintf('\n========================================\n');
    fprintf('compute_hopf_freq_stage | Stage: %s\n', cfg.stage_key);
    fprintf('Loading tseries from: %s\n', fullfile(paths.data_root, cfg.tseries_file));
    fprintf('Saving results to:   %s\n', results_dir);
    fprintf('========================================\n');

    % -----------------------------------------
    % Load data
    % -----------------------------------------
    S = load(fullfile(paths.data_root, cfg.tseries_file));
    if ~isfield(S, 'tseries')
        error('Variable "tseries" not found in file: %s', fullfile(paths.data_root, cfg.tseries_file));
    end
    tseries = S.tseries;

    % -----------------------------------------
    % Parameters
    % -----------------------------------------
    TR        = cfg.TR;
    NPARCELLS = cfg.NPARCELLS;
    NCOND     = cfg.NCOND;

    % -----------------------------------------
    % Bandpass filter settings
    % -----------------------------------------
    fnq = 1 / (2 * TR);
    flp = 0.008;
    fhi = 0.08;
    Wn  = [flp / fnq, fhi / fnq];
    k   = 2;
    [bfilt, afilt] = butter(k, Wn);

    % -----------------------------------------
    % Load stage mask
    % -----------------------------------------
    [stage_mask, PTID] = get_stage_mask(cfg, paths);

    if ~any(stage_mask)
        error('No subjects found for stage %s.', cfg.stage_key);
    end

    % -----------------------------------------
    % Loop over conditions
    % -----------------------------------------
    for cond = 1:NCOND

        xs_all = tseries(:, cond);

        % Keep only non-empty subjects first
        non_empty_idx = find(~cellfun(@isempty, xs_all));
        xs_nonempty   = xs_all(non_empty_idx);

        % Align stage mask with non-empty entries
        stage_mask_nonempty = stage_mask(non_empty_idx);

        % Select only subjects belonging to this stage
        xs = xs_nonempty(stage_mask_nonempty);
        NSUB = numel(xs);

        if NSUB == 0
            error('No non-empty subjects found for stage %s, condition %d.', cfg.stage_key, cond);
        end

        % Infer common time length from first selected subject
        first_ts = xs{1};
        [Ns0, TT0] = size(first_ts);

        if Ns0 ~= NPARCELLS
            error('First subject in stage %s has %d parcels; expected %d.', ...
                cfg.stage_key, Ns0, NPARCELLS);
        end

        Ts0     = TT0 * TR;
        freq    = (0:TT0/2 - 1) / Ts0;
        nfreqs0 = length(freq);

        fce_stage  = zeros(NSUB, NPARCELLS, NPARCELLS);
        PowSpect   = zeros(nfreqs0, NPARCELLS, NSUB);

        % -----------------------------------------
        % Subject loop
        % -----------------------------------------
        for sub = 1:NSUB
            fprintf('%s | cond %d | sub %d/%d\n', cfg.stage_key, cond, sub, NSUB);

            ts = xs{sub};
            [Ns, TT] = size(ts);

            if Ns ~= NPARCELLS
                error('Subject %d in stage %s has %d parcels; expected %d.', ...
                    sub, cfg.stage_key, Ns, NPARCELLS);
            end

            if TT ~= TT0
                error(['Inconsistent time length detected in stage %s, condition %d. ', ...
                       'Expected %d time points, found %d for subject %d.'], ...
                       cfg.stage_key, cond, TT0, TT, sub);
            end

            tss = zeros(NPARCELLS, TT);

            for seed = 1:NPARCELLS
                x = detrend(ts(seed,:) - mean(ts(seed,:)));
                tss(seed,:) = filtfilt(bfilt, afilt, x);

                pw = abs(fft(tss(seed,:)));
                PowSpect(:, seed, sub) = pw(1:floor(TT/2)).^2 / (TT / TR);
            end

            fce_stage(sub,:,:) = corrcoef(tss', 'rows', 'pairwise');
        end

        % -----------------------------------------
        % Stage-level empirical FC
        % -----------------------------------------
        fce = squeeze(mean(fce_stage, 1));

        % -----------------------------------------
        % Stage-level dominant frequencies
        % -----------------------------------------
        Power_Areas = squeeze(mean(PowSpect, 3));

        for seed = 1:NPARCELLS
            Power_Areas(:, seed) = gaussfilt(freq, Power_Areas(:, seed)', 0.01);
        end

        [~, index] = max(Power_Areas);
        f_diff = freq(index);

        if any(f_diff ~= 0)
            f_diff(f_diff == 0) = mean(f_diff(f_diff ~= 0));
        end

        % -----------------------------------------
        % Save outputs
        % -----------------------------------------
        outfile = fullfile(results_dir, ...
            sprintf('results_f_diff_fce_cond%d_%s_sch1000.mat', cond, cfg.save_prefix));

        save(outfile, 'f_diff', 'fce');

        fprintf('Saved:\n  %s\n', outfile);
    end

    fprintf('Done: %s\n', cfg.stage_key);
end