%% Compute HC–AD turbulence difference (lam = 0.01, Schaefer-1000)
% Input: Turbu_ComBat_ADNI3_HC_AD_ABeta_lam1_sch1000_N145.xlsx
% Output: diff_hc_ad_lam1_sch1000_ComBat_ADNI3.mat
%         variable: diff_hc_ad_lam1 (1 x 1000 double)

clear; clc;

% --- Paths ---
in_xlsx  = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/code/turbulence_fulldenoising/sch1000_N238rev/harmonization_allfeat/Turbu_ComBat_ADNI3_HC_AD_ABeta_lam1_sch1000_N145.xlsx';
out_mat  = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/neuromaps_analysis/neuromaps-data/annotations/turbu/MNI152/brain_nodes_diff_hc_ad_lam1_sch1000_ComBat_N145.mat';

% --- Read table ---
T = readtable(in_xlsx);

% Expect columns: PTID, Group, Schaefer_1 ... Schaefer_1000
% Check column names (optional sanity check)
disp(T.Properties.VariableNames(1:10));

% Identify Schaefer columns
isSchaefer = startsWith(T.Properties.VariableNames, 'Schaefer_');
schaeferCols = T.Properties.VariableNames(isSchaefer);

% Convert to numeric array (subjects x parcels)
X = table2array(T(:, schaeferCols));   % size: [Nsubj x 1000]

% --- Define groups ---
% HC group:   Group == 'HC_ABneg'
% AD group:   Group == 'AD_ABpos'
hc_idx = strcmp(T.Group, 'HC_ABneg');
ad_idx = strcmp(T.Group, 'AD_ABpos');

X_HC = X(hc_idx, :);   % HC subjects x 1000
X_AD = X(ad_idx, :);   % AD subjects x 1000

fprintf('N(HC) = %d\n', sum(hc_idx));
fprintf('N(AD) = %d\n', sum(ad_idx));

% --- Compute parcel-wise group means ---
HC_mean = mean(X_HC, 1, 'omitnan');   % 1 x 1000
AD_mean = mean(X_AD, 1, 'omitnan');   % 1 x 1000

% --- Difference: HC - AD ---
diff_hc_ad_lam1 = HC_mean - AD_mean;  % 1 x 1000

% Optional: check a few values
disp('First 5 values of diff_hc_ad_lam1:');
disp(diff_hc_ad_lam1(1:5));

% --- Save to .mat ---
save(out_mat, 'diff_hc_ad_lam1');

fprintf('Saved %s with variable diff_hc_ad_lam1 (size: %d x %d)\n', ...
    out_mat, size(diff_hc_ad_lam1,1), size(diff_hc_ad_lam1,2));
