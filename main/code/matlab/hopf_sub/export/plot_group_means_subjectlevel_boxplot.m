%% plot_subjectlevel_boxplots_replication.m
% Subject-level replication plots using BOXPLOTS
% - For each subject: mean across perturbation trials
% - For each group: show distribution of subject means
%
% Expected folder structure:
%   <baseDir>/<group>/Wtrials_<trial>_<condition>_<measure>_ADNI3_<group>_GEC_sch1000_SUB<subject>.mat
%
% Representative files contain:
%   nsub, s, condition, G, infocapacity, susceptibility

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
baseDir   = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/results/Wtrials';
trials    = 100;
condition = 1;
measure   = 'err_hete';

groups = {'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'};
nSub   = [54, 39, 33, 26];

nPermutations = 1000;
doSave = false;

%% ---------------- OUTPUT DIRECTORY ----------------
outFigDir = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/plots/subject-level';

if doSave
    if ~exist(outFigDir,'dir')
        mkdir(outFigDir);
    end
end

%% ---------------- BASIC CHECKS ----------------
if numel(groups) ~= numel(nSub)
    error('groups and nSub must have the same length.');
end

if any(isnan(nSub)) || any(nSub <= 0)
    error('Please set valid subject counts in nSub.');
end

if ~exist(baseDir, 'dir')
    error('baseDir does not exist:\n%s', baseDir);
end

for k = 1:numel(groups)
    thisGroupDir = fullfile(baseDir, groups{k});
    if ~exist(thisGroupDir, 'dir')
        warning('Group folder not found: %s', thisGroupDir);
    end
end

%% ---------------- AGGREGATE SUBJECT-LEVEL DATA ----------------
G = cell(1, numel(groups));

for k = 1:numel(groups)
    G{k} = aggregate_subject_level_group( ...
        groups{k}, baseDir, nSub(k), trials, condition, measure);
end

%% ---------------- EXTRACT SUBJECT-LEVEL VARIABLES ----------------
% Information capacity
HC_ABneg_infocap  = G{1}.IC_subject_means(~isnan(G{1}.IC_subject_means));
HC_ABpos_infocap  = G{2}.IC_subject_means(~isnan(G{2}.IC_subject_means));
MCI_ABpos_infocap = G{3}.IC_subject_means(~isnan(G{3}.IC_subject_means));
AD_ABpos_infocap  = G{4}.IC_subject_means(~isnan(G{4}.IC_subject_means));

% Susceptibility
HC_ABneg_susc  = G{1}.S_subject_means(~isnan(G{1}.S_subject_means));
HC_ABpos_susc  = G{2}.S_subject_means(~isnan(G{2}.S_subject_means));
MCI_ABpos_susc = G{3}.S_subject_means(~isnan(G{3}.S_subject_means));
AD_ABpos_susc  = G{4}.S_subject_means(~isnan(G{4}.S_subject_means));

%% ---------------- PRINT SUMMARY ----------------
fprintf('\n=========== SUBJECT-LEVEL SUMMARY ===========\n');
fprintf('Trials per subject: %d | condition: %d | measure: %s\n\n', ...
    trials, condition, measure);

fprintf('--- Information Capacity ---\n');
fprintf('Group\t\tExpectedN\tValidN\tMean\t\tSD\t\tSEM\n');
fprintf('---------------------------------------------------------------\n');
allInfocap = {HC_ABneg_infocap, HC_ABpos_infocap, MCI_ABpos_infocap, AD_ABpos_infocap};
for k = 1:numel(groups)
    x = allInfocap{k};
    fprintf('%s\t%d\t\t%d\t%.6f\t%.6f\t%.6f\n', ...
        groups{k}, nSub(k), numel(x), mean(x), std(x), std(x)/sqrt(numel(x)));
end

fprintf('\n--- Susceptibility ---\n');
fprintf('Group\t\tExpectedN\tValidN\tMean\t\tSD\t\tSEM\n');
fprintf('---------------------------------------------------------------\n');
allSusc = {HC_ABneg_susc, HC_ABpos_susc, MCI_ABpos_susc, AD_ABpos_susc};
for k = 1:numel(groups)
    x = allSusc{k};
    fprintf('%s\t%d\t\t%d\t%.6f\t%.6f\t%.6f\n', ...
        groups{k}, nSub(k), numel(x), mean(x), std(x), std(x)/sqrt(numel(x)));
end
fprintf('=============================================\n\n');

%% ---------------- COMMON PAIRWISE CONTRASTS ----------------
pairLabels = { ...
    'HC_ABneg vs HC_ABpos', ...
    'HC_ABpos vs MCI_ABpos', ...
    'MCI_ABpos vs AD_ABpos', ...
    'HC_ABneg vs MCI_ABpos', ...
    'HC_ABneg vs AD_ABpos'};

sigPairs = {[1,2],[2,3],[3,4],[1,3],[1,4]};

%% ============================================================
%% INFORMATION CAPACITY
%% ============================================================
figure('Name', 'Information Capacity (Subject-level boxplot)');
boxplot([HC_ABneg_infocap; HC_ABpos_infocap; MCI_ABpos_infocap; AD_ABpos_infocap], ...
    [ones(size(HC_ABneg_infocap)); ...
     2*ones(size(HC_ABpos_infocap)); ...
     3*ones(size(MCI_ABpos_infocap)); ...
     4*ones(size(AD_ABpos_infocap))], ...
    'Labels', groups);

title('Information capacity');
ylabel('Subject-level mean information capacity');
box off;
grid on;

p_IC       = nan(1,5);
obsDiff_IC = nan(1,5);
effSize_IC = nan(1,5);

[pval, observeddifference, effectsize] = permutationTest(HC_ABneg_infocap, HC_ABpos_infocap, nPermutations);
p_IC(1) = pval; obsDiff_IC(1) = observeddifference; effSize_IC(1) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(HC_ABpos_infocap, MCI_ABpos_infocap, nPermutations);
p_IC(2) = pval; obsDiff_IC(2) = observeddifference; effSize_IC(2) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(MCI_ABpos_infocap, AD_ABpos_infocap, nPermutations);
p_IC(3) = pval; obsDiff_IC(3) = observeddifference; effSize_IC(3) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(HC_ABneg_infocap, MCI_ABpos_infocap, nPermutations);
p_IC(4) = pval; obsDiff_IC(4) = observeddifference; effSize_IC(4) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(HC_ABneg_infocap, AD_ABpos_infocap, nPermutations);
p_IC(5) = pval; obsDiff_IC(5) = observeddifference; effSize_IC(5) = effectsize;

sigstar(sigPairs, p_IC);

[h_IC, ~, ~, adj_p_IC] = fdr_bh(p_IC, 0.05, 'pdep', 'yes');

fprintf('=========== PAIRWISE PERMUTATION TESTS: INFORMATION CAPACITY ===========\n');
fprintf('Comparison\t\t\tRaw p\t\tAdj p\t\tObservedDiff\tEffectSize\tReject H0\n');
fprintf('----------------------------------------------------------------------------------------------\n');
for i = 1:5
    fprintf('%s\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n', ...
        pairLabels{i}, p_IC(i), adj_p_IC(i), obsDiff_IC(i), effSize_IC(i), h_IC(i));
end
fprintf('==============================================================================================\n\n');

if doSave
    print(gcf, fullfile(outFigDir, 'InfoCapacity_ADNI3_subject_level_boxplot_perm1000'), '-dpng');
    saveas(gcf, fullfile(outFigDir, 'InfoCapacity_ADNI3_subject_level_boxplot_perm1000.fig'));
end

%% ============================================================
%% SUSCEPTIBILITY
%% ============================================================
figure('Name', 'Susceptibility (Subject-level boxplot)');
boxplot([HC_ABneg_susc; HC_ABpos_susc; MCI_ABpos_susc; AD_ABpos_susc], ...
    [ones(size(HC_ABneg_susc)); ...
     2*ones(size(HC_ABpos_susc)); ...
     3*ones(size(MCI_ABpos_susc)); ...
     4*ones(size(AD_ABpos_susc))], ...
    'Labels', groups);

title('Susceptibility');
ylabel('Subject-level mean susceptibility');
box off;
grid on;

p_S       = nan(1,5);
obsDiff_S = nan(1,5);
effSize_S = nan(1,5);

[pval, observeddifference, effectsize] = permutationTest(HC_ABneg_susc, HC_ABpos_susc, nPermutations);
p_S(1) = pval; obsDiff_S(1) = observeddifference; effSize_S(1) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(HC_ABpos_susc, MCI_ABpos_susc, nPermutations);
p_S(2) = pval; obsDiff_S(2) = observeddifference; effSize_S(2) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(MCI_ABpos_susc, AD_ABpos_susc, nPermutations);
p_S(3) = pval; obsDiff_S(3) = observeddifference; effSize_S(3) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(HC_ABneg_susc, MCI_ABpos_susc, nPermutations);
p_S(4) = pval; obsDiff_S(4) = observeddifference; effSize_S(4) = effectsize;

[pval, observeddifference, effectsize] = permutationTest(HC_ABneg_susc, AD_ABpos_susc, nPermutations);
p_S(5) = pval; obsDiff_S(5) = observeddifference; effSize_S(5) = effectsize;

sigstar(sigPairs, p_S);

[h_S, ~, ~, adj_p_S] = fdr_bh(p_S, 0.05, 'pdep', 'yes');

fprintf('=========== PAIRWISE PERMUTATION TESTS: SUSCEPTIBILITY ===========\n');
fprintf('Comparison\t\t\tRaw p\t\tAdj p\t\tObservedDiff\tEffectSize\tReject H0\n');
fprintf('------------------------------------------------------------------------------------------\n');
for i = 1:5
    fprintf('%s\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n', ...
        pairLabels{i}, p_S(i), adj_p_S(i), obsDiff_S(i), effSize_S(i), h_S(i));
end
fprintf('==========================================================================================\n\n');

if doSave
    print(gcf, fullfile(outFigDir, 'Susceptibility_ADNI3_subject_level_boxplot_perm1000'), '-dpng');
    saveas(gcf, fullfile(outFigDir, 'Susceptibility_ADNI3_subject_level_boxplot_perm1000.fig'));
end

%% ---------------- SAVE STATISTICAL OUTPUTS ----------------
if doSave
    save(fullfile(outFigDir, 'SubjectLevel_BoxplotStats_ADNI3.mat'), ...
        'groups', 'nSub', 'nPermutations', ...
        'HC_ABneg_infocap', 'HC_ABpos_infocap', 'MCI_ABpos_infocap', 'AD_ABpos_infocap', ...
        'HC_ABneg_susc', 'HC_ABpos_susc', 'MCI_ABpos_susc', 'AD_ABpos_susc', ...
        'p_IC', 'obsDiff_IC', 'effSize_IC', 'h_IC', 'adj_p_IC', ...
        'p_S', 'obsDiff_S', 'effSize_S', 'h_S', 'adj_p_S');
end