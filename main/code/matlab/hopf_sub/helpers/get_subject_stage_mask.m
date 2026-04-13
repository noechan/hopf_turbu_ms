function [stage_mask, PTID_valid] = get_stage_mask(cfg, paths, valid_mask)
% GET_STAGE_MASK
%
% Returns a logical mask selecting subjects belonging to the requested stage.
% The mask is returned in the order of PTIDs after applying valid_mask
% (e.g., non-empty subject data).

    abeta_table = readtable(paths.abeta_table);

    S = load(fullfile(paths.data_root, cfg.ptid_file));
    if ~isfield(S, 'PTID')
        error('Variable "PTID" not found in file: %s', ...
            fullfile(paths.data_root, cfg.ptid_file));
    end

    PTID_all = cellstr(S.PTID);

    if nargin < 3 || isempty(valid_mask)
        valid_mask = true(size(PTID_all));
    end

    PTID_valid = PTID_all(valid_mask);

    stage_ptids = abeta_table.PTID( ...
        strcmp(abeta_table.GROUP, cfg.parent_group) & ...
        abeta_table.ABeta_pvc == cfg.abeta_value ...
    );

    stage_mask = ismember(PTID_valid, stage_ptids);
end