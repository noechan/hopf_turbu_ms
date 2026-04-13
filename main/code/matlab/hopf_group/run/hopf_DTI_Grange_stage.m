function hopf_DTI_Grange_stage(stage_key, s, condition)
% HOPF_DTI_GRANGE_STAGE
%
% Stage-driven refactor of hopf_DTI_Grange_* for Slurm/HPC execution.
%
% Keeps the original numerical simulation core as close as possible to the
% stage-specific version, but removes hard-coded stage names and redundant loading.
%
% Inputs:
%   stage_key : 'HC_ABneg', 'HC_ABpos', 'MCI_ABpos', 'AD_ABpos'
%   s         : index into G_range
%   condition : usually 1
%
% Output:
%   Saves one file per G index:
%   WG_%03d_%d_<save_prefix>_GEC_sch1000.mat

    paths = get_paths();
    cfg   = get_stage_config(stage_key);

    addpath(genpath(paths.code_root));

    % -------------------------------------------------
    % Create results directory
    % -------------------------------------------------
    results_dir = fullfile(paths.hopf_grange_root, cfg.results_subdir);
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fprintf('\n========================================\n');
    fprintf('hopf_DTI_Grange_stage | Stage: %s | s=%d | condition=%d\n', ...
        cfg.stage_key, s, condition);
    fprintf('Saving results to: %s\n', results_dir);
    fprintf('========================================\n');

    % -------------------------------------------------
    % Load stage-specific inputs
    % -------------------------------------------------
    freq_file = fullfile(paths.freq_root, cfg.results_subdir, ...
        sprintf('results_f_diff_fce_cond%d_%s_sch1000.mat', condition, cfg.save_prefix));

    empcorr_file = fullfile(paths.empcorr_root, cfg.results_subdir, ...
        sprintf('empirical_spacorr_rest_cond_%d_%s_sch1000.mat', condition, cfg.save_prefix));

    gec_file = fullfile(paths.data_root, cfg.gec_file);

    if ~isfile(freq_file)
        error('Frequency file not found:\n%s', freq_file);
    end

    if ~isfile(empcorr_file)
        error('Empirical spatial correlation file not found:\n%s', empcorr_file);
    end

    if ~isfile(gec_file)
        error('GEC file not found:\n%s', gec_file);
    end

    % Frequency / empirical FC
    Sfreq = load(freq_file);
    if ~isfield(Sfreq, 'f_diff')
        error('Variable "f_diff" not found in:\n%s', freq_file);
    end
    f_diff = Sfreq.f_diff;

    % Empirical spatial correlation
    Semp = load(empcorr_file);
    if ~isfield(Semp, 'corrfcn')
        error('Variable "corrfcn" not found in:\n%s', empcorr_file);
    end
    empcorrfcn = Semp.corrfcn;

    % GEC
    Sgec = load(gec_file);
    if isfield(Sgec, 'Ceffgroup')
        Ceffgroup = Sgec.Ceffgroup;
    else
        error('Variable "Ceffgroup" not found in:\n%s', gec_file);
    end

    % -------------------------------------------------
    % Load precomputed geometry
    % -------------------------------------------------
    Spre = load(paths.precomp_file, ...
        'rr', 'range', 'delta', 'xrange', 'LAMBDA', 'NLAMBDA', 'C1');

    rr      = Spre.rr;
    range   = Spre.range; %#ok<NASGU>
    delta   = Spre.delta;
    xrange  = Spre.xrange; %#ok<NASGU>
    LAMBDA  = Spre.LAMBDA;
    NLAMBDA = Spre.NLAMBDA;
    C1      = Spre.C1;

    % -------------------------------------------------
    % Parameters and G range
    % -------------------------------------------------
    NPARCELLS = cfg.NPARCELLS;
    NR        = cfg.NR;
    NRini     = cfg.NRini;
    NRfin     = cfg.NRfin;
    NSUBSIM   = cfg.NSUBSIM;

    lambda = cfg.lambda;
    [~, indsca] = min(abs(LAMBDA - lambda));

    G_range =0:0.01:3;

    if s < 1 || s > numel(G_range)
        error('s=%d is out of bounds for G_range of length %d.', s, numel(G_range));
    end

    G = G_range(s);

    % -------------------------------------------------
    % Parameters of the data
    % -------------------------------------------------
    TR = cfg.TR;

    % Bandpass filter settings
    fnq = 1/(2*TR);
    flp = 0.008;
    fhi = 0.08;
    Wn  = [flp/fnq fhi/fnq];
    k   = 2;
    [bfilt, afilt] = butter(k, Wn);

    % -------------------------------------------------
    % Parameters HOPF
    % -------------------------------------------------
    Tmax  = cfg.Tmax;
    omega = repmat(2*pi*f_diff', 1, 2);
    omega(:,1) = -omega(:,1);

    dt   = 0.1*TR/2;
    sig  = 0.01;
    dsig = sqrt(dt)*sig;

    % -------------------------------------------------
    % Preallocation
    % -------------------------------------------------
    corrfcn = zeros(NPARCELLS, NR);
    lam_mean_spatime_enstrophy = zeros(NLAMBDA, NPARCELLS, Tmax);
    err_hete    = zeros(1, NSUBSIM);
    InfoFlow    = zeros(NSUBSIM, NLAMBDA-1);
    InfoCascade = zeros(1, NSUBSIM);
    Turbulence  = zeros(1, NSUBSIM);
    mutinfo     = zeros(1, NLAMBDA);

    % Connectivity scaling
    factor = max(max(Ceffgroup));
    C = Ceffgroup / factor * 0.2;

    fprintf('G = %.4f\n', G);

    % =================================================
    % FROM HERE ON: ORIGINAL MODEL SIMULATIONS
    % =================================================
    for sub = 1:NSUBSIM
        fprintf('%s | G index %d/%d | sim %d/%d\n', ...
            cfg.stage_key, s, numel(G_range), sub, NSUBSIM);

        wC = G * C;
        sumC = repmat(sum(wC,2), 1, 2);

        %% Hopf Simulation
        a  = -0.02 * ones(NPARCELLS, 2);
        xs = zeros(Tmax, NPARCELLS);

        z  = 0.1 * ones(NPARCELLS, 2); % x = z(:,1), y = z(:,2)
        nn = 0;

        % discard first 2000 time steps
        for t = 0:dt:2000
            suma = wC*z - sumC.*z;
            zz   = z(:,end:-1:1);
            z    = z + dt*(a.*z + zz.*omega - z.*(z.*z + zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        end

        % actual modeling
        for t = 0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z;
            zz   = z(:,end:-1:1);
            z    = z + dt*(a.*z + zz.*omega - z.*(z.*z + zz.*zz) + suma) + dsig*randn(NPARCELLS,2);

            if abs(mod(t,TR)) < 0.01
                nn = nn + 1;
                xs(nn,:) = z(:,1)';
            end
        end

        ts = xs';

        %% Computing simulated FC
        signal_filt = zeros(NPARCELLS, Tmax);
        Phases      = zeros(NPARCELLS, Tmax);

        for seed = 1:NPARCELLS
            ts(seed,:) = detrend(ts(seed,:) - mean(ts(seed,:)));
            signal_filt(seed,:) = filtfilt(bfilt, afilt, ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end

        fcsimul = corrcoef(signal_filt');

        %% Computing simulated observables
        for i = 1:NPARCELLS
            numind   = zeros(1, NR);
            corrfcn_1 = zeros(1, NR);

            for j = 1:NPARCELLS
                r = rr(i,j);
                index = floor(r/delta) + 1;

                if index == NR + 1
                    index = NR;
                end

                mcc = fcsimul(i,j);

                if ~isnan(mcc)
                    corrfcn_1(index) = corrfcn_1(index) + mcc;
                    numind(index) = numind(index) + 1;
                end
            end

            corrfcn(i,:) = corrfcn_1 ./ numind;

            % enstrophy
            ilam = 1;
            for lam = LAMBDA %#ok<NASGU>
                enstrophy = nansum(repmat(squeeze(C1(ilam,i,:)), 1, Tmax) .* ...
                    complex(cos(Phases), sin(Phases))) / sum(C1(ilam,i,:));
                lam_mean_spatime_enstrophy(ilam,i,:) = abs(enstrophy);
                ilam = ilam + 1;
            end
        end

        Rspatime = squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
        Rsub = nanstd(Rspatime(:));

        ilam = 1;
        for lam = LAMBDA %#ok<NASGU>
            if ilam > 1
                [cc, pp] = corr( ...
                    (squeeze(lam_mean_spatime_enstrophy(ilam,:,2:end)))', ...
                    (squeeze(lam_mean_spatime_enstrophy(ilam-1,:,1:end-1)))' );
                mutinfo(ilam) = nanmean(cc(pp(:) < 0.05));
            end
            ilam = ilam + 1;
        end

        %% Computing difference between simulated and empirical FC
        err1  = zeros(1, NPARCELLS);
        err11 = zeros(1, NR);

        for i = 1:NPARCELLS
            for k = NRini:NRfin
                err11(k) = (corrfcn(i,k) - empcorrfcn(i,k))^2;
            end
            err1(i) = nanmean(err11(NRini:NRfin));
        end

        err_hete(sub) = sqrt(nanmean(err1));

        %% Observables
        Inflam2 = mutinfo(2:end);
        InfoFlow(sub,:) = Inflam2;
        InfoCascade(sub) = nanmean(Inflam2);
        Turbulence(sub) = Rsub;
    end

    % -------------------------------------------------
    % Save
    % -------------------------------------------------
    outfile = fullfile(results_dir, ...
        sprintf('WG_%03d_%d_%s_GEC_sch1000.mat', s, condition, cfg.save_prefix));

    save(outfile, 'InfoCascade', 'InfoFlow', 'Turbulence', 'err_hete');

    fprintf('Saved:\n  %s\n', outfile);
end