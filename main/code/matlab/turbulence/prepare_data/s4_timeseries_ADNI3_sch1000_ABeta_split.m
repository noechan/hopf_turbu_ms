clear all; close all
addpath(genpath('/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/code/turbulence_fulldenoising'))
output_path = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/code/turbulence_fulldenoising/sch1000_N238rev/figures_N152/sch1000/Abeta_Status/raw/';

%% Load ABeta status table
abeta_table = readtable('/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/participants/cross_sectional/ADNI3_N238rev_with_ABETA_Status_CL24.xlsx');
abeta_table.Properties.VariableNames;

%% Load PTIDs
ptid_path='/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/code/turbulence_fulldenoising/sch1000_N238rev/prepare_data/';
load(fullfile(ptid_path, 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat'));    % PTID
PTID_HC = cellstr(PTID);
load(fullfile(ptid_path, 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat'));   % PTID
PTID_MCI = cellstr(PTID);
load(fullfile(ptid_path, 'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat'));    % PTID
PTID_AD = cellstr(PTID);

%% Helper: masks by Aβ status aligned to group PTIDs
get_abeta_masks = @(group_ids, gname) deal( ...
    ismember(group_ids, abeta_table.PTID(strcmpi(abeta_table.GROUP, gname) & abeta_table.ABeta_pvc==1)), ... % Aβ+
    ismember(group_ids, abeta_table.PTID(strcmpi(abeta_table.GROUP, gname) & abeta_table.ABeta_pvc==0)) ...  % Aβ−
    );

% Keep HC Aβ− and HC Aβ+; keep only MCI Aβ+
[isAbetaPos_HC, isAbetaNeg_HC] = get_abeta_masks(PTID_HC,  'HC');
[isAbetaPos_MCI, ~]            = get_abeta_masks(PTID_MCI, 'MCI');
[isAbetaPos_AD, ~]            = get_abeta_masks(PTID_AD, 'AD');


%% Load & save results
data_path='/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/timeseries/inputs_fulldenoising/';
cd(fullfile(data_path,'ADNI3_HC_MPRAGE_IRFSPGR_N238rev'))
tseries_HC=  load('tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat');

tseries_HC_ABneg=tseries_HC.tseries(isAbetaNeg_HC,:);
tseries_HC_ABpos=tseries_HC.tseries(isAbetaPos_HC,:);
save('tseries_ADNI3_HC_ABetaNeg_MPRAGE_IRFSPGR_sch1000_N238rev.mat','tseries_HC_ABneg')
save('tseries_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_sch1000_N238rev.mat','tseries_HC_ABpos')


cd(fullfile(data_path,'ADNI3_MCI_MPRAGE_IRFSPGR_N238rev'))
tseries_MCI= load('tseries_ADNI3_MCI_MPRAGE_IRFSPGR_sch1000_N238rev.mat');

tseries_MCI_ABpos=tseries_MCI.tseries(isAbetaPos_MCI,:);
save('tseries_ADNI3_MCI_ABetaPos_MPRAGE_IRFSPGR_sch1000_N238rev.mat','tseries_MCI_ABpos')


cd(fullfile(data_path,'ADNI3_AD_MPRAGE_IRFSPGR_N238rev'))
tseries_AD=  load('tseries_ADNI3_AD_MPRAGE_IRFSPGR_sch1000_N238rev');

tseries_AD_ABpos=tseries_AD.tseries(isAbetaPos_AD,:);
save('tseries_ADNI3_AD_ABetaPos_MPRAGE_IRFSPGR_sch1000_N238rev.mat','tseries_AD_ABpos')
