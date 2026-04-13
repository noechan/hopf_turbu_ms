function [tseries_stage, ptid_stage, meta] = select_stage_timeseries(tseries, PTID, abeta_table, stage_label)
%SELECT_STAGE_TIMESERIES Select subject time series for one clinical/Aβ stage.
%
% Inputs
%   tseries      : cell array of subject time series
%   PTID         : PTID list aligned with tseries
%   abeta_table  : table with at least PTID, GROUP, ABeta_pvc
%   stage_label  : one of
%                  'HC_ABneg', 'HC_ABpos', 'MCI_ABneg', 'MCI_ABpos', 'AD_ABpos'
%
% Outputs
%   tseries_stage : subset of tseries for the requested stage
%   ptid_stage    : PTIDs of the selected subjects
%   meta          : struct with group/status/mask bookkeeping

arguments
    tseries
    PTID
    abeta_table table
    stage_label char
end

% Ensure PTID is a cellstr for robust matching
if isstring(PTID)
    PTID = cellstr(PTID);
elseif ischar(PTID)
    PTID = cellstr(PTID);
elseif ~iscell(PTID)
    PTID = cellstr(PTID);
end

switch stage_label
    case 'HC_ABneg'
        group_name = 'HC';
        abeta_value = 0;
    case 'HC_ABpos'
        group_name = 'HC';
        abeta_value = 1;
    case 'MCI_ABneg'
        group_name = 'MCI';
        abeta_value = 0;
    case 'MCI_ABpos'
        group_name = 'MCI';
        abeta_value = 1;
    case 'AD_ABpos'
        group_name = 'AD';
        abeta_value = 1;
    otherwise
        error('Unknown stage_label: %s', stage_label);
end

required_vars = {'PTID','GROUP','ABeta_pvc'};
missing_vars = required_vars(~ismember(required_vars, abeta_table.Properties.VariableNames));
if ~isempty(missing_vars)
    error('ABeta table is missing required variable(s): %s', strjoin(missing_vars, ', '));
end

stage_ptids = abeta_table.PTID( ...
    strcmp(abeta_table.GROUP, group_name) & ...
    abeta_table.ABeta_pvc == abeta_value ...
    );

mask = ismember(PTID, stage_ptids);

tseries_stage = tseries(mask, :);
ptid_stage = PTID(mask);

meta = struct();
meta.stage_label = stage_label;
meta.group_name = group_name;
meta.abeta_value = abeta_value;
meta.mask = mask;
meta.n_selected = sum(mask);
end