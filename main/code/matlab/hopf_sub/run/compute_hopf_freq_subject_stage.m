function compute_hopf_freq_subject_stage(stage_key)
% COMPUTE_HOPF_FREQ_SUBJECT_STAGE
%
% Stage-specific subject-level Hopf frequency computation.
%
% Valid stage_key values:
%   - 'HC_ABneg'
%   - 'HC_ABpos'
%   - 'MCI_ABpos'
%   - 'AD_ABpos'

    paths = get_subject_paths();
    cfg   = get_subject_stage_config(stage_key);

    addpath(genpath(paths.code_root));

    results_dir = fullfile(paths.freq_root, cfg.results_subdir);
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fprintf('\n========================================\n');
    fprintf('compute_hopf_freq_subject_stage | Stage: %s\n', cfg.stage_key);
    fprintf('Loading tseries: %s\n', fullfile(paths.data_root, cfg.tseries_file));
    fprintf('Loading geometry: %s\n', paths.precomp_file);
    fprintf('Saving to:       %s\n', results_dir);
    fprintf('========================================\n');

    % --------------------------------------------------
    % Load tseries
    % --------------------------------------------------
    S = load(fullfile(paths.data_root, cfg.tseries_file));
    if ~isfield(S, 'tseries')
        error('Variable "tseries" not found in %s', ...
            fullfile(paths.data_root, cfg.tseries_file));
    end
    tseries = S.tseries;

    % --------------------------------------------------
    % Load precomputed geometry
    % --------------------------------------------------
    if ~isfile(paths.precomp_file)
        error(['Precomputed geometry file not found: %s\n' ...
               'Please run precompute_hopf_geometry first.'], ...
               paths.precomp_file);
    end

    G = load(paths.precomp_file, 'rr', 'delta', 'xrange');

    required_geom = {'rr', 'delta', 'xrange'};
    for iVar = 1:numel(required_geom)
        if ~isfield(G, required_geom{iVar})
            error('Variable "%s" not found in %s', ...
                required_geom{iVar}, paths.precomp_file);
        end
    end

    rr     = G.rr;
    delta  = G.delta;
    xrange = G.xrange;

    % --------------------------------------------------
    % Parameters
    % --------------------------------------------------
    TR         = cfg.TR;
    NPARCELLS  = cfg.NPARCELLS;
    NCOND      = cfg.NCOND;
    TT_target  = cfg.TT_target;
    NR         = cfg.NR;

    if NCOND ~= 1
        error('This script currently expects NCOND = 1.');
    end

    if size(rr,1) ~= NPARCELLS || size(rr,2) ~= NPARCELLS
        error('Geometry matrix rr has size [%d x %d], expected [%d x %d].', ...
            size(rr,1), size(rr,2), NPARCELLS, NPARCELLS);
    end

    if numel(xrange) ~= NR
        error('xrange has length %d, expected %d.', numel(xrange), NR);
    end

    % --------------------------------------------------
    % Filter settings
    % --------------------------------------------------
    fnq = 1 / (2 * TR);
    Wn  = [cfg.flp / fnq, cfg.fhi / fnq];
    [bfilt, afilt] = butter(cfg.butter_order, Wn);

    % --------------------------------------------------
    % Frequency grid
    % --------------------------------------------------
    Ts     = TT_target * TR;
    freq   = (0:TT_target/2 - 1) ./ Ts;
    nfreqs = numel(freq);

    % --------------------------------------------------
    % Extract condition and remove empty subjects
    % --------------------------------------------------
    xs_all = tseries(:, 1);
    valid_mask = ~cellfun(@isempty, xs_all);
    xs_valid   = xs_all(valid_mask);

    % --------------------------------------------------
    % Stage-specific mask
    % --------------------------------------------------
    [stage_mask, PTID_valid] = get_subject_stage_mask(cfg, paths, valid_mask);

    if numel(stage_mask) ~= numel(PTID_valid)
        error('Stage mask and PTID_valid have inconsistent lengths.');
    end

    if numel(xs_valid) ~= numel(PTID_valid)
        error('Mismatch between xs_valid and PTID_valid lengths.');
    end

    xs_stage   = xs_valid(stage_mask);
    PTID_stage = PTID_valid(stage_mask);
    NSUB       = numel(xs_stage);

    fprintf('Valid subjects: %d\n', numel(xs_valid));
    fprintf('Selected for stage %s: %d\n', cfg.stage_key, NSUB);

    if NSUB == 0
        error('No subjects found for stage %s.', cfg.stage_key);
    end

    % --------------------------------------------------
    % Preallocate
    % --------------------------------------------------
    FCemp_SUB    = zeros(NSUB, NPARCELLS, NPARCELLS);
    PowSpect_SUB = zeros(nfreqs, NPARCELLS, NSUB);
    f_diff_SUB   = zeros(NSUB, NPARCELLS);
    corrfcn_SUB  = zeros(NSUB, NPARCELLS, NR);

    % --------------------------------------------------
    % Subject loop
    % --------------------------------------------------
    parfor ii = 1:NSUB

        ts = xs_stage{ii};
        [Ns, Tsub] = size(ts);

        if Ns ~= NPARCELLS
            error('Subject %d in stage %s has %d parcels; expected %d.', ...
                ii, cfg.stage_key, Ns, NPARCELLS);
        end

        % Detrend + filter + enforce common TT_target
        tss = zeros(NPARCELLS, TT_target);

        for p = 1:NPARCELLS
            x = ts(p, :);
            x = detrend(x - mean(x, 'omitnan'));
            x(~isfinite(x)) = 0;

            xf = filtfilt(bfilt, afilt, x);

            if Tsub >= TT_target
                tss(p, :) = xf(1:TT_target);
            else
                pad = zeros(1, TT_target);
                pad(1:Tsub) = xf;
                tss(p, :) = pad;
            end
        end

        % ---------------- FC ----------------
        FC_emp = corrcoef(tss', 'Rows', 'pairwise');
        FCemp_SUB(ii,:,:) = FC_emp;

        % ---------------- Spectra + peak frequency ----------------
        pw = abs(fft(tss, [], 2)).^2 / (TT_target / TR);
        pw = pw(:, 1:floor(TT_target/2));

        local_PowSpect = zeros(nfreqs, NPARCELLS);
        local_f_diff   = zeros(1, NPARCELLS);

        for p = 1:NPARCELLS
            sm = gaussfilt(freq, squeeze(pw(p,:))', cfg.gauss_sigma);
            local_PowSpect(:, p) = sm;

            [~, idx_peak] = max(sm);
            local_f_diff(p) = freq(idx_peak);
        end

        PowSpect_SUB(:,:,ii) = local_PowSpect;
        f_diff_SUB(ii,:)     = local_f_diff;

        % ---------------- Spatial correlation function ----------------
        local_corrfcn = zeros(NPARCELLS, NR);

        for i = 1:NPARCELLS
            numind    = zeros(1, NR);
            corrfcn_1 = zeros(1, NR);

            for j = 1:NPARCELLS
                r     = rr(i,j);
                index = floor(r / delta) + 1;

                if index == NR + 1
                    index = NR;
                end

                mcc = FC_emp(i,j);

                if ~isnan(mcc)
                    corrfcn_1(index) = corrfcn_1(index) + mcc;
                    numind(index)    = numind(index) + 1;
                end
            end

            corrfcn_1(numind > 0) = corrfcn_1(numind > 0) ./ numind(numind > 0);
            local_corrfcn(i,:)    = corrfcn_1;
        end

        corrfcn_SUB(ii,:,:) = local_corrfcn;
    end

    % --------------------------------------------------
    % Defensive correction for zero-frequency peaks
    % --------------------------------------------------
    mask_zero = (f_diff_SUB == 0);
    if any(mask_zero(:))
        nonzero_vals = f_diff_SUB(~mask_zero);
        if ~isempty(nonzero_vals)
            f_diff_SUB(mask_zero) = mean(nonzero_vals, 'omitnan');
        end
    end

    % --------------------------------------------------
    % Stage-level summaries
    % --------------------------------------------------
    FCemp_mean    = squeeze(mean(FCemp_SUB, 1, 'omitnan'));
    PowSpect_mean = mean(PowSpect_SUB, 3, 'omitnan');

    [~, idx_stage] = max(PowSpect_mean, [], 1);
    f_diff_stage = freq(idx_stage);

    if any(f_diff_stage == 0)
        nz = f_diff_stage(f_diff_stage ~= 0);
        if ~isempty(nz)
            f_diff_stage(f_diff_stage == 0) = mean(nz, 'omitnan');
        end
    end

    corrfcn_mean = squeeze(mean(corrfcn_SUB, 1, 'omitnan'));

    % --------------------------------------------------
    % Save subject-level file
    % --------------------------------------------------
    outfile_sub = fullfile(results_dir, ...
        sprintf('results_f_diff_fce_cond1_%s_sch1000_SUB.mat', cfg.save_prefix));

    save(outfile_sub, ...
        'freq', 'xrange', ...
        'PTID_stage', ...
        'FCemp_SUB', 'PowSpect_SUB', 'f_diff_SUB', 'corrfcn_SUB', ...
        'FCemp_mean', 'PowSpect_mean', 'f_diff_stage', 'corrfcn_mean', ...
        '-v7.3');

    fprintf('Saved subject-level file:\n  %s\n', outfile_sub);

    % --------------------------------------------------
    % Group-level file
    % --------------------------------------------------
    if cfg.save_group_compatibility
        outfile_grp = fullfile(results_dir, ...
            sprintf('results_f_diff_fce_cond1_%s_sch1000.mat', cfg.save_prefix));

        f_diff = f_diff_stage;
        fce    = FCemp_mean;

        save(outfile_grp, 'f_diff', 'fce');

        fprintf('Saved compatibility group-level file:\n  %s\n', outfile_grp);
    end

    fprintf('Done: %s\n', cfg.stage_key);
end