function compute_empirical_spacorr_subject_stage(stage_key)
% COMPUTE_EMPIRICAL_SPACORR_SUBJECT_STAGE
%
% Stage-specific subject-level empirical spatial correlation computation.
%
% Valid stage_key values:
%   - 'HC_ABneg'
%   - 'HC_ABpos'
%   - 'MCI_ABpos'
%   - 'AD_ABpos'
%
% This function:
%   1. Loads stage-specific FCemp_SUB from the frequency step
%   2. Loads precomputed geometry (rr, xrange, delta)
%   3. Computes subject-level spatial correlation functions
%   4. Saves:
%        - subject-level CorrFcn_SUB
%        - group-level mean corrfcn
%
% Output files are saved under:
%   results/hopf_subject/empirical_spacorr/<STAGE>/

    paths = get_subject_paths();
    cfg   = get_subject_stage_config(stage_key);

    addpath(genpath(paths.code_root));

    results_dir = fullfile(paths.empcorr_root, cfg.results_subdir);
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    % --------------------------------------------------
    % Input file from frequency step
    % --------------------------------------------------
    sub_level_file = fullfile(paths.freq_root, cfg.results_subdir, ...
        sprintf('results_f_diff_fce_cond1_%s_sch1000_SUB.mat', cfg.save_prefix));

    fprintf('\n========================================\n');
    fprintf('compute_empirical_spacorr_subject_stage | Stage: %s\n', cfg.stage_key);
    fprintf('Loading FC file:  %s\n', sub_level_file);
    fprintf('Loading geometry: %s\n', paths.precomp_file);
    fprintf('Saving to:        %s\n', results_dir);
    fprintf('========================================\n');

    % --------------------------------------------------
    % Load stage-specific subject-level FC
    % --------------------------------------------------
    if ~isfile(sub_level_file)
        error(['Stage-specific frequency output not found: %s\n' ...
               'Please run compute_hopf_freq_subject_stage first.'], ...
               sub_level_file);
    end

    S = load(sub_level_file, 'FCemp_SUB', 'PTID_stage');
    if ~isfield(S, 'FCemp_SUB')
        error('Variable "FCemp_SUB" not found in file: %s', sub_level_file);
    end
    if ~isfield(S, 'PTID_stage')
        error('Variable "PTID_stage" not found in file: %s', sub_level_file);
    end

    FCemp_SUB  = S.FCemp_SUB;
    PTID_stage = S.PTID_stage;

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
    % Basic parameters
    % --------------------------------------------------
    NPARCELLS = cfg.NPARCELLS;
    NR        = cfg.NR;
    NCOND     = cfg.NCOND;

    if NCOND ~= 1
        error('This script currently expects NCOND = 1.');
    end

    NSUB = size(FCemp_SUB, 1);

    if size(FCemp_SUB,2) ~= NPARCELLS || size(FCemp_SUB,3) ~= NPARCELLS
        error(['FCemp_SUB has unexpected size [%d x %d x %d]; expected [NSUB x %d x %d].'], ...
            size(FCemp_SUB,1), size(FCemp_SUB,2), size(FCemp_SUB,3), ...
            NPARCELLS, NPARCELLS);
    end

    if size(rr,1) ~= NPARCELLS || size(rr,2) ~= NPARCELLS
        error('rr has unexpected size [%d x %d]; expected [%d x %d].', ...
            size(rr,1), size(rr,2), NPARCELLS, NPARCELLS);
    end

    if numel(xrange) ~= NR
        error('xrange has length %d; expected %d.', numel(xrange), NR);
    end

    if numel(PTID_stage) ~= NSUB
        error('PTID_stage length (%d) does not match NSUB (%d).', ...
            numel(PTID_stage), NSUB);
    end

    fprintf('Number of subjects in stage %s: %d\n', cfg.stage_key, NSUB);

    % --------------------------------------------------
    % Allocate subject-level container
    % CorrFcn_SUB: [NSUB x NPARCELLS x NR]
    % --------------------------------------------------
    CorrFcn_SUB = zeros(NSUB, NPARCELLS, NR);

    % --------------------------------------------------
    % Loop over conditions (kept for naming compatibility)
    % --------------------------------------------------
    for cond = 1:NCOND

        fprintf('Computing corrfcn for condition %d...\n', cond);

        % ---------------- SUBJECT-LEVEL LOOP ----------------
        parfor s = 1:NSUB
            FCs = squeeze(FCemp_SUB(s,:,:));  % [NPARCELLS x NPARCELLS]

            corrfcn_sub_i = nan(NPARCELLS, NR);

            for ii = 1:NPARCELLS
                numind    = zeros(1, NR);
                corrfcn_1 = zeros(1, NR);

                for jj = 1:NPARCELLS
                    r = rr(ii,jj);

                    index = floor(r / delta) + 1;
                    if index > NR
                        index = NR;
                    end

                    mcc = FCs(ii,jj);
                    if ~isnan(mcc)
                        corrfcn_1(index) = corrfcn_1(index) + mcc;
                        numind(index)    = numind(index) + 1;
                    end
                end

                valid_bins = numind > 0;
                corrfcn_tmp = nan(1, NR);
                corrfcn_tmp(valid_bins) = corrfcn_1(valid_bins) ./ numind(valid_bins);

                corrfcn_sub_i(ii,:) = corrfcn_tmp;
            end

            CorrFcn_SUB(s,:,:) = corrfcn_sub_i;
        end

        % --------------------------------------------------
        % Group-level mean over subjects
        % --------------------------------------------------
        corrfcn_mean = squeeze(mean(CorrFcn_SUB, 1, 'omitnan'));

        % Backward-compatible variable name
        corrfcn = corrfcn_mean; %#ok<NASGU>

        % --------------------------------------------------
        % Save subject-level output
        % --------------------------------------------------
        outfile_sub = fullfile(results_dir, ...
            sprintf('empirical_spacorr_rest_cond_%d_%s_sch1000_SUB.mat', ...
            cond, cfg.save_prefix));

        save(outfile_sub, ...
            'xrange', 'rr', ...
            'PTID_stage', ...
            'CorrFcn_SUB', ...
            '-v7.3');

        fprintf('Saved subject-level spatial correlations:\n  %s\n', outfile_sub);

        % --------------------------------------------------
        % Save group-level mean output
        % --------------------------------------------------
        outfile_grp = fullfile(results_dir, ...
            sprintf('empirical_spacorr_rest_cond_%d_%s_sch1000.mat', ...
            cond, cfg.save_prefix));

        save(outfile_grp, 'xrange', 'rr', 'corrfcn');

        fprintf('Saved group-level mean spatial correlation:\n  %s\n', outfile_grp);
    end

    fprintf('Done: %s\n', cfg.stage_key);
end