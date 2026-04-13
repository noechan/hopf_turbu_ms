function cfg = get_subject_stage_config(stage_key)
% GET_SUBJECT_STAGE_CONFIG
%
% Returns configuration for one subject-level stage.

    switch upper(stage_key)

        case 'HC_ABNEG'
            cfg.stage_key      = 'HC_ABneg';
            cfg.parent_group   = 'HC';
            cfg.group_label    = 'HC';
            cfg.abeta_value    = 0;
            cfg.abeta_label    = 'ABneg';

            cfg.tseries_file   = 'tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file      = 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat';

            cfg.save_prefix    = 'ADNI3_HC_ABneg';
            cfg.results_subdir = 'HC_ABneg';

        case 'HC_ABPOS'
            cfg.stage_key      = 'HC_ABpos';
            cfg.parent_group   = 'HC';
            cfg.group_label    = 'HC';
            cfg.abeta_value    = 1;
            cfg.abeta_label    = 'ABpos';

            cfg.tseries_file   = 'tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file      = 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat';

            cfg.save_prefix    = 'ADNI3_HC_ABpos';
            cfg.results_subdir = 'HC_ABpos';

        case 'MCI_ABPOS'
            cfg.stage_key      = 'MCI_ABpos';
            cfg.parent_group   = 'MCI';
            cfg.group_label    = 'MCI';
            cfg.abeta_value    = 1;
            cfg.abeta_label    = 'ABpos';

            cfg.tseries_file   = 'tseries_ADNI3_MCI_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file      = 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat';

            cfg.save_prefix    = 'ADNI3_MCI_ABpos';
            cfg.results_subdir = 'MCI_ABpos';

        case 'AD_ABPOS'
            cfg.stage_key      = 'AD_ABpos';
            cfg.parent_group   = 'AD';
            cfg.group_label    = 'AD';
            cfg.abeta_value    = 1;
            cfg.abeta_label    = 'ABpos';

            cfg.tseries_file   = 'tseries_ADNI3_AD_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file      = 'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat';

            cfg.save_prefix    = 'ADNI3_AD_ABpos';
            cfg.results_subdir = 'AD_ABpos';

        otherwise
            error('Unknown stage_key: %s', stage_key);
    end

    % ----------------------------
    % Stage filtering meta
    % ----------------------------
    cfg.requires_abeta = true;

    % ----------------------------
    % Shared parameters
    % ----------------------------
    cfg.TR         = 3;
    cfg.NPARCELLS  = 1000;
    cfg.NCOND      = 1;
    cfg.TT_target  = 197;
    cfg.NR         = 400;
    cfg.gauss_sigma = 0.01;

    % Band-pass
    cfg.flp = 0.008;
    cfg.fhi = 0.08;
    cfg.butter_order = 2;

    % Compatibility
    cfg.save_group_compatibility = true;
end