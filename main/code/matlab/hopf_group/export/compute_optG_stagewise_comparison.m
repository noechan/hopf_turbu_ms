function compute_optG_stagewise_comparison(condition)
% COMPUTE_OPTG_STAGEWISE_COMPARISON
%
% Stage-driven refactor to compute optimal G for all ADNI3 staging groups.
%
% Inputs:
%   condition : usually 1
%
% Outputs:
%   Saves one optG file per stage in:
%   results/hopf_group/optG/<stage_key>/optG_<stage_key>.mat
%
%   Also saves a comparative figure across stages.

    if nargin < 1
        condition = 1;
    end

    clearvars -except condition
    clc

    paths = get_paths();
    addpath(genpath(paths.code_root));

    G = 100; %#ok<NASGU> % kept for compatibility with original loader functions
    G_range = 0:0.01:3;

    stage_keys = {'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'};

    stage_colors = containers.Map( ...
        {'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'}, ...
        {[0.7 0.7 0.7], [0.4940 0.1840 0.5560], [0.8500 0.3250 0.0980], [0 0.4470 0.7410]} ...
    );

    results = struct();

    for i = 1:numel(stage_keys)
        stage_key = stage_keys{i};
        cfg = get_stage_config(stage_key);

        fprintf('========================================\n');
        fprintf('Computing optG | Stage: %s | condition=%d\n', stage_key, condition);
        fprintf('========================================\n');

        % ---------------------------------------------
        % Stage-specific G_range directory
        % ---------------------------------------------
        grange_dir = fullfile(paths.hopf_grange_root, stage_key);
        if ~isfolder(grange_dir)
            error('G_range directory not found:\n%s', grange_dir);
        end
        addpath(grange_dir);

        % ---------------------------------------------
        % Load stage condition using stage-specific loader
        % ---------------------------------------------
        loader_name = sprintf('load_condition_%s', stage_key);
        if exist(loader_name, 'file') ~= 2
            error('Loader function not found: %s', loader_name);
        end

        cond_struct = feval(loader_name, G, condition);

        if ~isfield(cond_struct, 'err_hete_range')
            error('Field "err_hete_range" not found for stage %s', stage_key);
        end

        [minerrhete, ioptG] = min(mean(cond_struct.err_hete_range, 2));
        optG = G_range(ioptG);

        results.(stage_key).cond_struct = cond_struct;
        results.(stage_key).minerrhete  = minerrhete;
        results.(stage_key).ioptG       = ioptG;
        results.(stage_key).optG        = optG;
        results.(stage_key).color       = stage_colors(stage_key);

        % ---------------------------------------------
        % Save optG
        % ---------------------------------------------
        optG_dir = fullfile(paths.optG_root, stage_key);
        if ~isfolder(optG_dir)
            mkdir(optG_dir);
        end

        save_file = fullfile(optG_dir, sprintf('optG_%s.mat', stage_key));
        save(save_file, 'optG');

        fprintf('Saved optG for %s: %.4f\n', stage_key, optG);
        fprintf('File: %s\n', save_file);
    end

    % ---------------------------------------------
    % Comparative plot
    % ---------------------------------------------
    x = G_range(1:50)';

    figure;
    hold on

    for i = 1:numel(stage_keys)
        stage_key = stage_keys{i};
        c = results.(stage_key).color;
        err_hete_range = results.(stage_key).cond_struct.err_hete_range;

        shadedErrorBar( ...
            x, err_hete_range(1:50,:)', {@mean,@std}, ...
            'lineprops', {'k-','color',c,'markerfacecolor',c} ...
        );

        xline(results.(stage_key).optG, 'Color', c, 'LineWidth', 1.5);
    end

    ylabel('Fitting');
    xlabel('Global coupling G');
    legend(stage_keys, 'Location', 'northeastoutside');
    xticks([0 0.1 0.2 0.3 0.4 0.5]);

    % ---------------------------------------------
    % Save figure
    % ---------------------------------------------
    fig_dir = paths.hopf_grange_root;
    if ~isfolder(fig_dir)
        mkdir(fig_dir);
    end

    png_file = fullfile(fig_dir, 'Global coupling Fitting ADNI3 Staging 4(N=238).png');
    fig_file = fullfile(fig_dir, 'Global coupling Fitting ADNI3 Staging 4 (N=238).fig');

    print(gcf, png_file, '-dpng');
    saveas(gcf, fig_file);

    fprintf('Comparative figure saved to:\n%s\n%s\n', png_file, fig_file);
end