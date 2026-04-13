%% plot_groupmeans_subjectlevel_replication.m
% Replication plot: subject-level perturbations aggregated to GROUP MEANS
% - For each subject: mean across 100 perturbation trials
% - For each group: mean across subjects (and SEM across subjects)
%
% Expected folder structure:
%   <baseDir>/<group>/Wtrials_<trial>_<condition>_<measure>_ADNI3_<group>_GEC_sch1000_SUB<subject>.mat
%
% Example:
%   /.../results/Wtrials/HC_ABneg/Wtrials_001_1_err_hete_ADNI3_HC_ABneg_GEC_sch1000_SUB001.mat
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

% Subject counts per group (same order as "groups")
nSub = [54, 39, 33, 26];

% Plot options
doSave = false;
outFigDir = fullfile(baseDir, 'plots', 'subject_level');
if doSave && ~exist(outFigDir, 'dir')
    mkdir(outFigDir);
end

%% ---------------- BASIC CHECKS ----------------
if numel(groups) ~= numel(nSub)
    error('groups and nSub must have the same length.');
end

if any(isnan(nSub)) || any(nSub <= 0)
    error('Please set valid subject counts in nSub = [N_HCneg, N_HCpos, N_MCIpos, N_ADpos].');
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

%% ---------------- AGGREGATE ----------------
G = cell(1, numel(groups));

for k = 1:numel(groups)
    G{k} = aggregate_subject_level_group( ...
        groups{k}, baseDir, nSub(k), trials, condition, measure);
end

%% ---------------- EXTRACT MEANS + SEMS ----------------
IC_means = nan(numel(groups), 1);
IC_sems  = nan(numel(groups), 1);
S_means  = nan(numel(groups), 1);
S_sems   = nan(numel(groups), 1);

for k = 1:numel(groups)
    IC_means(k) = G{k}.IC_group_mean;
    IC_sems(k)  = G{k}.IC_group_sem;
    S_means(k)  = G{k}.S_group_mean;
    S_sems(k)   = G{k}.S_group_sem;
end

%% ---------------- PRINT REPLICATION SUMMARY ----------------
fprintf('\n=========== SUBJECT-LEVEL -> GROUP MEANS (Replication Check) ===========\n');
fprintf('Trials per subject: %d | condition: %d | measure: %s\n\n', ...
    trials, condition, measure);
fprintf('Group\t\tNsub\tIC_mean\t\tIC_SEM\t\tS_mean\t\tS_SEM\n');
fprintf('-----------------------------------------------------------------------\n');

for k = 1:numel(groups)
    fprintf('%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
        groups{k}, nSub(k), IC_means(k), IC_sems(k), S_means(k), S_sems(k));
end

fprintf('=======================================================================\n\n');

%% ---------------- PLOT 1: INFORMATION CAPACITY ----------------
figure('Name', 'Info Capacity (Group Means)');
hold on;

x = 1:numel(groups);

bar(x, IC_means);
errorbar(x, IC_means, IC_sems, 'LineStyle', 'none');

set(gca, 'XTick', x, 'XTickLabel', groups);
xtickangle(20);
ylabel('Information capacity (group mean of subject means)');
title('Replication plot (subject-level aggregation): Information Capacity');
box off;
grid on;

if doSave
    saveas(gcf, fullfile(outFigDir, 'InfoCapacity_GROUPMEANS_subject_level.png'));
    savefig(gcf, fullfile(outFigDir, 'InfoCapacity_GROUPMEANS_subject_level.fig'));
end

%% ---------------- PLOT 2: SUSCEPTIBILITY ----------------
figure('Name', 'Susceptibility (Group Means)');
hold on;

bar(x, S_means);
errorbar(x, S_means, S_sems, 'LineStyle', 'none');

set(gca, 'XTick', x, 'XTickLabel', groups);
xtickangle(20);
ylabel('Susceptibility (group mean of subject means)');
title('Replication plot (subject-level aggregation): Susceptibility');
box off;
grid on;

if doSave
    saveas(gcf, fullfile(outFigDir, 'Susceptibility_GROUPMEANS_subject_level.png'));
    savefig(gcf, fullfile(outFigDir, 'Susceptibility_GROUPMEANS_subject_level.fig'));
end