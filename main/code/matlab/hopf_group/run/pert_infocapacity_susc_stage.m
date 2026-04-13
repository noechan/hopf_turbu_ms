function pert_infocapacity_susc_stage(stage_key, s, condition)
% PERT_INFOCAPACITY_SUSC_STAGE
%
% Stage-driven refactor of pert_infocapacity_susc_*_GEC for Slurm/HPC execution.
%
% Keeps the original numerical simulation core as close as possible to the
% stage-specific version, but removes hard-coded stage names and redundant loading.
%
% Inputs:
%   stage_key : 'HC_ABneg', 'HC_ABpos', 'MCI_ABpos', 'AD_ABpos'
%   s         : random seed / Slurm array index
%   condition : usually 1
%
% Output:
%   Saves one file per seed:
%   Wtrials_%03d_%d_err_hete_<save_prefix>_GEC_sch1000.mat

    paths = get_paths();
    cfg   = get_stage_config(stage_key);

    addpath(genpath(paths.code_root));

    % -------------------------------------------------
    % Create results directory
    % -------------------------------------------------
    results_dir = fullfile(paths.results_root, 'W_trials', cfg.results_subdir);
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fprintf('\n========================================\n');
    fprintf('pert_infocapacity_susc_stage | Stage: %s | s=%d | condition=%d\n', ...
        cfg.stage_key, s, condition);
    fprintf('Saving results to: %s\n', results_dir);
    fprintf('========================================\n');

    rng(s);
    measure = 'err_hete';

    % -------------------------------------------------
    % Build stage-specific input files
    % -------------------------------------------------
    optG_file = fullfile(paths.optG_root, cfg.stage_key, sprintf('optG_%s.mat', cfg.stage_key));

    freq_file = fullfile(paths.freq_root, cfg.results_subdir, ...
        sprintf('results_f_diff_fce_cond%d_%s_sch1000.mat', condition, cfg.save_prefix));

    empcorr_file = fullfile(paths.empcorr_root, cfg.results_subdir, ...
        sprintf('empirical_spacorr_rest_cond_%d_%s_sch1000.mat', condition, cfg.save_prefix));

    gec_file = fullfile(paths.data_root, cfg.gec_file);
    tseries_file = fullfile(paths.data_root, cfg.tseries_file);
    ptid_file = fullfile(paths.data_root, cfg.ptid_file);
    schaefer_file = fullfile(paths.data_root, 'SchaeferCOG.mat');
    abeta_file = paths.abeta_table;

    % -------------------------------------------------
    % Validate files
    % -------------------------------------------------
    required_files = {optG_file, freq_file, empcorr_file, gec_file, ...
                      tseries_file, ptid_file, schaefer_file, abeta_file};

    for i = 1:numel(required_files)
        if ~isfile(required_files{i})
            error('Required file not found:\n%s', required_files{i});
        end
    end

    % -------------------------------------------------
    % Load stage-specific inputs
    % -------------------------------------------------
    if condition == 1
        Sopt = load(optG_file);
        opt_field = sprintf('optG_%s', cfg.stage_key);

        if isfield(Sopt, opt_field)
            G = Sopt.(opt_field);
        else
            fn = fieldnames(Sopt);
            if numel(fn) == 1
                G = Sopt.(fn{1});
            else
                error('Could not identify optimal G variable in:\n%s', optG_file);
            end
        end
    else
        error('Unsupported condition=%d. Only condition==1 is currently implemented.', condition);
    end

    Sfreq = load(freq_file);
    if ~isfield(Sfreq, 'f_diff')
        error('Variable "f_diff" not found in:\n%s', freq_file);
    end
    f_diff = Sfreq.f_diff;

    Semp = load(empcorr_file);
    if ~isfield(Semp, 'corrfcn')
        error('Variable "corrfcn" not found in:\n%s', empcorr_file);
    end
    empcorrfcn = Semp.corrfcn; %#ok<NASGU>

    Sgec = load(gec_file);
    if isfield(Sgec, 'Ceffgroup')
        Ceffgroup = Sgec.Ceffgroup;
    else
        error('Variable "Ceffgroup" not found in:\n%s', gec_file);
    end

    Sschaefer = load(schaefer_file);
    if isfield(Sschaefer, 'SchaeferCOG')
        SchaeferCOG = Sschaefer.SchaeferCOG;
    else
        error('Variable "SchaeferCOG" not found in:\n%s', schaefer_file);
    end

    Sts = load(tseries_file);
    if isfield(Sts, 'tseries')
        tseries = Sts.tseries;
    else
        fn = fieldnames(Sts);
        if numel(fn) == 1
            tseries = Sts.(fn{1});
        else
            error('Could not identify tseries variable in:\n%s', tseries_file);
        end
    end

    Sptid = load(ptid_file);
    if isfield(Sptid, 'PTID')
        PTID = cellstr(Sptid.PTID);
    else
        error('Variable "PTID" not found in:\n%s', ptid_file);
    end

    abeta_table = readtable(abeta_file);

    % -------------------------------------------------
    % Filter tseries by Aβ status using stage config
    % -------------------------------------------------
    is_group = strcmp(abeta_table.GROUP, cfg.group_label);
    is_abeta = abeta_table.ABeta_pvc == cfg.abeta_value;
    stage_ptids = abeta_table.PTID(is_group & is_abeta);

    stage_mask = ismember(PTID, stage_ptids);
    tseries = tseries(stage_mask,:); %#ok<NASGU>

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
    % Parameters and lambda
    % -------------------------------------------------
    NPARCELLS = cfg.NPARCELLS;
    NR        = cfg.NR;
    NRini     = cfg.NRini; %#ok<NASGU>
    NRfin     = cfg.NRfin; %#ok<NASGU>
    NSUBSIM   = cfg.NSUBSIM;

    lambda = cfg.lambda;
    [~, indsca] = min(abs(LAMBDA - lambda));

    % -------------------------------------------------
    % Parameters of the data
    % -------------------------------------------------
    TR = cfg.TR;

    fnq = 1/(2*TR);
    flp = 0.008;
    fhi = 0.08;
    Wn  = [flp/fnq fhi/fnq];
    k   = 2;
    [bfilt, afilt] = butter(k, Wn);
    Isubdiag = find(tril(ones(NPARCELLS),-1)); %#ok<NASGU>

    % -------------------------------------------------
    % Parameters HOPF
    % -------------------------------------------------
    Tmax  = cfg.Tmax;
    omega = repmat(2*pi*f_diff',1,2);
    omega(:,1) = -omega(:,1);

    dt   = 0.1*TR/2;
    sig  = 0.01;
    dsig = sqrt(dt)*sig;

    % -------------------------------------------------
    % Preallocation
    % -------------------------------------------------
    lam_mean_spatime_enstrophy = zeros(NLAMBDA, NPARCELLS, Tmax);
    ensspasub  = zeros(NSUBSIM, NPARCELLS);
    ensspasub1 = zeros(NSUBSIM, NPARCELLS);

    factor = max(max(Ceffgroup));
    C = Ceffgroup / factor * 0.2;

    fprintf('G = %.4f\n', G);

    % =================================================
    % FROM HERE ON: ORIGINAL MODEL SIMULATIONS
    % =================================================
    for sub = 1:NSUBSIM
        fprintf('%s | seed %d | sim %d/%d\n', ...
            cfg.stage_key, s, sub, NSUBSIM);

        wC   = G * C;
        sumC = repmat(sum(wC,2),1,2);

        %% HOPF SIMULATION
        a  = -0.02 * ones(NPARCELLS,2);
        xs = zeros(Tmax,NPARCELLS);
        z  = 0.1 * ones(NPARCELLS,2);
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

        %% Compute the Kuramoto local parameter of the unperturbed model
        signal_filt = zeros(NPARCELLS, Tmax);
        Phases      = zeros(NPARCELLS, Tmax);

        for seed = 1:NPARCELLS
            ts(seed,:) = detrend(ts(seed,:) - mean(ts(seed,:)));
            signal_filt(seed,:) = filtfilt(bfilt, afilt, ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end

        for i = 1:NPARCELLS
            ilam = 1;
            for lam = LAMBDA %#ok<NASGU>
                enstrophy = nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax) .* ...
                    complex(cos(Phases),sin(Phases))) / sum(C1(ilam,i,:));
                lam_mean_spatime_enstrophy(ilam,i,:) = abs(enstrophy);
                ilam = ilam + 1;
            end
        end

        Rspatime = squeeze(lam_mean_spatime_enstrophy(indsca,:,:));
        ensspasub(sub,:) = (nanmean(Rspatime,2))';

        %% HERE WE INTRODUCE THE PERTURBATION
        a  = -0.02 + 0.02 * repmat(rand(NPARCELLS,1),1,2) .* ones(NPARCELLS,2);
        nn = 0;

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
        Rspatime1 = zeros(NPARCELLS,Tmax);

        %% Compute the Kuramoto local order parameter of the perturbed model
        signal_filt = zeros(NPARCELLS, Tmax);
        Phases      = zeros(NPARCELLS, Tmax);

        for seed = 1:NPARCELLS
            ts(seed,:) = detrend(ts(seed,:) - mean(ts(seed,:)));
            signal_filt(seed,:) = filtfilt(bfilt, afilt, ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end

        for i = 1:NPARCELLS
            enstrophy = nansum(repmat(squeeze(C1(indsca,i,:)),1,Tmax) .* ...
                complex(cos(Phases),sin(Phases))) / sum(C1(indsca,i,:));
            Rspatime1(i,:) = abs(enstrophy);
        end

        ensspasub1(sub,:) = (nanmean(Rspatime1,2))';
    end

    %% Compute susceptibility and information capacity
    infocapacity   = nanmean(nanstd(ensspasub1 - ones(NSUBSIM,1)*nanmean(ensspasub)));
    susceptibility = nanmean(nanmean(ensspasub1 - ones(NSUBSIM,1)*nanmean(ensspasub)));

    % -------------------------------------------------
    % Save
    % -------------------------------------------------
    save_file = fullfile(results_dir, ...
        sprintf('Wtrials_%03d_%d_%s_%s_GEC_sch1000.mat', ...
        s, condition, measure, cfg.save_prefix));

    save(save_file, 'infocapacity', 'susceptibility');

    fprintf('Saved: %s\n', save_file);
end