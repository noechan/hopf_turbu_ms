function compute_perturbation_stagewise(trials, condition)
% SUMMARIZE_PERTURBATION_STAGEWISE
%
% Stage-driven refactor of the perturbation summary/statistics script.
%
% Inputs:
%   trials    : number of perturbation trials, usually 100
%   condition : usually 1
%
% Outputs:
%   1) Saves per-stage perturbation summary files
%   2) Saves boxplots for Information Capacity and Susceptibility
%   3) Saves raw permutation-test p-values and FDR-corrected results

    if nargin < 1
        trials = 100;
    end
    if nargin < 2
        condition = 1;
    end

    clearvars -except trials condition
    clc

    paths = get_paths();
    addpath(genpath(paths.code_root));

    stage_keys = {'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'};
    groups     = stage_keys;

    results_dir = fullfile(paths.results_root,'plots');
    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    fprintf('========================================\n');
    fprintf('summarize_perturbation_stagewise | trials=%d | condition=%d\n', ...
        trials, condition);
    fprintf('Saving results to: %s\n', results_dir);
    fprintf('========================================\n');

    % -------------------------------------------------
    % Load all stages using stage-specific loaders
    % -------------------------------------------------
    stage_data = struct();

    for i = 1:numel(stage_keys)
        stage_key = stage_keys{i};
        cfg = get_stage_config(stage_key); %#ok<NASGU>

        loader_name = sprintf('load_condition_InfoCap_ADNI3_%s', stage_key);
        if exist(loader_name, 'file') ~= 2
            error('Loader function not found: %s', loader_name);
        end

        fprintf('Loading perturbation summary for stage: %s\n', stage_key);
        tmp = feval(loader_name, trials, condition);

        if ~isfield(tmp, 'infocap_all')
            error('Field "infocap_all" not found for stage %s', stage_key);
        end
        if ~isfield(tmp, 'suscep_all')
            error('Field "suscep_all" not found for stage %s', stage_key);
        end

        stage_data.(stage_key) = tmp;
    end

    % -------------------------------------------------
    % Export per-stage perturbation summaries
    % -------------------------------------------------
    for i = 1:numel(stage_keys)
        stage_key = stage_keys{i};

        infocap_all = stage_data.(stage_key).infocap_all;
        suscep_all  = stage_data.(stage_key).suscep_all;

        save_file = fullfile(results_dir, sprintf('perturbation_%s.mat', stage_key));
        save(save_file, 'infocap_all', 'suscep_all');

        fprintf('Saved stage summary: %s\n', save_file);
    end

    % -------------------------------------------------
    % Build matrices for plotting/stats
    % -------------------------------------------------
    infocap_mat = [
        stage_data.HC_ABneg.infocap_all, ...
        stage_data.HC_ABpos.infocap_all, ...
        stage_data.MCI_ABpos.infocap_all, ...
        stage_data.AD_ABpos.infocap_all ...
    ];

    suscep_mat = [
        stage_data.HC_ABneg.suscep_all, ...
        stage_data.HC_ABpos.suscep_all, ...
        stage_data.MCI_ABpos.suscep_all, ...
        stage_data.AD_ABpos.suscep_all ...
    ];

    comparison_pairs = {
        [1,2], ...
        [2,3], ...
        [3,4], ...
        [1,3], ...
        [1,4]
    };

    % -------------------------------------------------
    % Information Capacity
    % -------------------------------------------------
    figure
    boxplot(infocap_mat, 'labels', groups)
    title('Information capacity');

    p1P = zeros(1,5);

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABneg.infocap_all, ...
        stage_data.HC_ABpos.infocap_all, 1000);
    p1P(1,1) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABpos.infocap_all, ...
        stage_data.MCI_ABpos.infocap_all, 1000);
    p1P(1,2) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.MCI_ABpos.infocap_all, ...
        stage_data.AD_ABpos.infocap_all, 1000);
    p1P(1,3) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABneg.infocap_all, ...
        stage_data.MCI_ABpos.infocap_all, 1000);
    p1P(1,4) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABneg.infocap_all, ...
        stage_data.AD_ABpos.infocap_all, 1000);
    p1P(1,5) = pval;

    sigstar(comparison_pairs, p1P);

    png_file = fullfile(results_dir, 'Info Capacity ADNI3 N=238 Permutation 1000 both.png');
    fig_file = fullfile(results_dir, 'Info Capacity ADNI3 N=238 Permutation 1000 both.fig');
    print(gcf, png_file, '-dpng');
    saveas(gcf, fig_file);
    close()

    [h_t, ~, ~, adj_p_t] = fdr_bh(p1P(1,1:5), 0.05, 'pdep', 'yes');
    save(fullfile(results_dir, 'InfoCap_rejectedH0.mat'), 'h_t', 'adj_p_t', 'p1P');

    % -------------------------------------------------
    % Susceptibility
    % -------------------------------------------------
    figure
    boxplot(suscep_mat, 'labels', groups)
    title('Susceptibility');

    p2P = zeros(1,5);

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABneg.suscep_all, ...
        stage_data.HC_ABpos.suscep_all, 1000);
    p2P(1,1) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABpos.suscep_all, ...
        stage_data.MCI_ABpos.suscep_all, 1000);
    p2P(1,2) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.MCI_ABpos.suscep_all, ...
        stage_data.AD_ABpos.suscep_all, 1000);
    p2P(1,3) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABneg.suscep_all, ...
        stage_data.MCI_ABpos.suscep_all, 1000);
    p2P(1,4) = pval;

    [pval, ~, ~] = permutationTest( ...
        stage_data.HC_ABneg.suscep_all, ...
        stage_data.AD_ABpos.suscep_all, 1000);
    p2P(1,5) = pval;

    sigstar(comparison_pairs, p2P);

    png_file = fullfile(results_dir, 'Susceptibility ADNI3 N=238 Permutation 1000 both.png');
    fig_file = fullfile(results_dir, 'Susceptibility ADNI3 N=238 Permutation 1000 both.fig');
    print(gcf, png_file, '-dpng');
    saveas(gcf, fig_file);
    close()

    [h_t, ~, ~, adj_p_t] = fdr_bh(p2P(1,1:5), 0.05, 'pdep', 'yes');
    save(fullfile(results_dir, 'InfoSus_rejectedH0.mat'), 'h_t', 'adj_p_t', 'p2P');

    fprintf('All perturbation summaries, figures, and stats saved to:\n%s\n', results_dir);
end