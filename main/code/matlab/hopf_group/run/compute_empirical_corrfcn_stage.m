function compute_empirical_corrfcn_stage(stage_key)
% COMPUTE_EMPIRICAL_CORRFCN_STAGE
%
% Stage-driven version of empirical spatial correlation function.
% Computes distance-dependent FC (corrfcn) for one stage:
%   - HC_ABneg
%   - HC_ABpos
%   - MCI_ABpos
%   - AD_ABpos
%
% Uses precomputed geometry (rr, delta, xrange).
%
% Outputs saved to:
%   results/hopf_group/empirical_spacorr/<STAGE>/

    paths = get_paths();
    cfg   = get_stage_config(stage_key);

    addpath(genpath(paths.code_root));

    % -----------------------------------------
    % Create results directory
    % -----------------------------------------
    results_dir = fullfile(paths.results_root, 'empirical_spacorr', cfg.results_subdir);

    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fprintf('\n========================================\n');
    fprintf('compute_empirical_corrfcn_stage | Stage: %s\n', cfg.stage_key);
    fprintf('========================================\n');

    % -----------------------------------------
    % Load tseries
    % -----------------------------------------
    S = load(fullfile(paths.data_root, cfg.tseries_file));
    tseries = S.tseries;

    % -----------------------------------------
    % Load precomputed geometry
    % -----------------------------------------
    P = load(paths.precomp_file, 'rr', 'delta', 'xrange');
    rr     = P.rr;
    delta  = P.delta;
    xrange = P.xrange;

    % -----------------------------------------
    % Parameters
    % -----------------------------------------
    TR        = cfg.TR;
    NPARCELLS = cfg.NPARCELLS;
    NR        = cfg.NR;
    NCOND     = cfg.NCOND;

    % -----------------------------------------
    % Bandpass filter
    % -----------------------------------------
    fnq = 1/(2*TR);
    flp = 0.008;
    fhi = 0.08;
    Wn  = [flp/fnq fhi/fnq];
    k   = 2;
    [bfilt, afilt] = butter(k, Wn);

    % -----------------------------------------
    % Load stage mask
    % -----------------------------------------
    [stage_mask, ~] = get_stage_mask(cfg, paths);

    if ~any(stage_mask)
        error('No subjects found for stage %s.', cfg.stage_key);
    end

    % -----------------------------------------
    % Loop over conditions
    % -----------------------------------------
    for cond = 1:NCOND

        xs_all = tseries(:, cond);

        % Remove empty subjects
        non_empty_idx = find(~cellfun(@isempty, xs_all));
        xs_nonempty   = xs_all(non_empty_idx);

        % Align mask
        stage_mask_nonempty = stage_mask(non_empty_idx);

        % Select stage subjects
        xs = xs_nonempty(stage_mask_nonempty);
        NSUB = numel(xs);

        if NSUB == 0
            error('No valid subjects for stage %s.', cfg.stage_key);
        end

        % Determine time length
        [~, TT] = size(xs{1});

        corrfcnsub = zeros(NSUB, NPARCELLS, NR);

        % -----------------------------------------
        % Subject loop
        % -----------------------------------------
        for sub = 1:NSUB
            fprintf('%s | cond %d | sub %d/%d\n', cfg.stage_key, cond, sub, NSUB);

            ts = xs{sub};

            signal_filt = zeros(NPARCELLS, TT);

            for seed = 1:NPARCELLS
                ts(seed,:) = detrend(ts(seed,:) - mean(ts(seed,:)));
                signal_filt(seed,:) = filtfilt(bfilt, afilt, ts(seed,:));
            end

            fce = corrcoef(signal_filt');

            % -----------------------------------------
            % Spatial correlation function
            % -----------------------------------------
            for i = 1:NPARCELLS
                numind = zeros(1, NR);
                corrfcn_1 = zeros(1, NR);

                for j = 1:NPARCELLS
                    r = rr(i,j);

                    index = floor(r/delta) + 1;
                    if index == NR + 1
                        index = NR;
                    end

                    mcc = fce(i,j);

                    if ~isnan(mcc)
                        corrfcn_1(index) = corrfcn_1(index) + mcc;
                        numind(index) = numind(index) + 1;
                    end
                end

                corrfcnsub(sub,i,:) = corrfcn_1 ./ numind;
            end
        end

        % -----------------------------------------
        % Average across subjects
        % -----------------------------------------
        corrfcn = squeeze(nanmean(corrfcnsub, 1));

        % -----------------------------------------
        % Save outputs
        % -----------------------------------------
        outfile = fullfile(results_dir, ...
            sprintf('empirical_spacorr_rest_cond_%d_%s_sch1000.mat', cond, cfg.save_prefix));

        save(outfile, 'corrfcn', 'corrfcnsub');

        fprintf('Saved:\n  %s\n', outfile);
    end

    fprintf('Done: %s\n', cfg.stage_key);
end