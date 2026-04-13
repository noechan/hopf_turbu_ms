clear all, close all
% S1_precompute_geometry.m
code_path = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/hofp_turbu_paper_project/main/code/matlab/hopf_subject/helpers/';
addpath(genpath(code_path))
data_path = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/hopf_turbu_paper_project/main/data/hopf_subject/';
output_path ='/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/hopf_turbu_paper_project/main/data/hopf_subject/precomputed_subject_geometry'; 
lambda_target = 0.18;

precompute_subject_hopf_geometry(data_path,output_path, lambda_target);