function cfg = get_gec_config()
%GET_GEC_CONFIG Centralised configuration for the ADNI3 GEC pipeline.

cfg = struct();

% Basic dimensions / condition
cfg.cnd = 1;
cfg.N = 1000;
cfg.NPARCELLS = cfg.N;

% Linear Hopf / optimisation parameters
cfg.Tau = 1;              % fit time-lagged statistics at 1 TR
cfg.sigma = 0.01;         % additive noise scale in the linear Hopf model
cfg.epsFC = 0.00018;      % learning rate for FC error
cfg.epsFCtau = 0.00036;   % learning rate for lagged-covariance error
cfg.maxC = 0.1;           % renormalise connectivity after each update
cfg.maxIter = 5000;
cfg.stop_check_every = 100;
cfg.stop_rel_improvement = 0.001;

% MRI acquisition / filtering
cfg.TR = 3;
cfg.flp = 0.008;
cfg.fhi = 0.08;
cfg.filter_order = 2;

% SC initialisation scaling
cfg.sc_init_max = 0.2;

% Preserve your original connectivity constraint exactly as written.
% I recommend checking whether this offset is the intended one for Sch1000.
cfg.homologue_offset = 200;

% Stage naming conventions
cfg.valid_stage_labels = {'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'};
end