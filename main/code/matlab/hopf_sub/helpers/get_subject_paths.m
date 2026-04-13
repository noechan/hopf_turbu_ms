function paths = get_subject_paths()
% GET_SUBJECT_PATHS
%
% Centralised path definition for the subject-level Hopf pipeline.

    paths.project_root = '/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/hopf_turbu_paper_project/main/';

    % Code
    paths.code_root    = fullfile(paths.project_root, 'matlab', 'code');

    % Subject-level data
    paths.data_root    = fullfile(paths.project_root, 'data', 'hopf_subject');
    paths.precomp_root = fullfile(paths.data_root, 'precomputed_subject_geometry');

    % Results
    paths.results_root = fullfile(paths.project_root, 'results', 'hopf_subject');
    paths.freq_root    = fullfile(paths.results_root, 'freq');
    paths.empcorr_root    = fullfile(paths.results_root, 'empirical_spacorr');

    % Shared metadata
    paths.abeta_table  = fullfile(paths.data_root, 'ADNI3_N238rev_with_ABETA_Status_CL24.xlsx');
    paths.cog_file     = fullfile(paths.data_root, 'SchaeferCOG.mat');

    %Precomputed geometry file for later steps
    paths.precomp_file = fullfile(paths.precomp_root, 'precomputed_subject_geometry.mat');
end