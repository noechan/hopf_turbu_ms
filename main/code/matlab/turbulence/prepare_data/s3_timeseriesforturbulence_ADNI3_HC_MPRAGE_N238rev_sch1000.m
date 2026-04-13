%% Create timeseries for turbulence

clear all
addpath(genpath('/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/code/turbulence_fulldenoising/')) %code path

% Get timeseries from batch 1
data_path='/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/timeseries/prepro_fulldenoising/'; 
input_path='/Volumes/ADNI/Projects/2024/ADNI/LONI_IDA/ADNI3/timeseries/inputs_fulldenoising/';
batch_path={'ADNI3_HC_MPRAGE_batch_1','ADNI3_HC_MPRAGE_batch_2','ADNI3_HC_MPRAGE_batch_3'};
tseries_path={'ADNI3_HC_MPRAGE_batch_1_N238rev','ADNI3_HC_MPRAGE_batch_2_N238rev','ADNI3_HC_MPRAGE_batch_3_N238rev'};

for batch=1:3

data_path_batch=fullfile(data_path,batch_path{1,batch});
fnames=dir(data_path_batch);
fnames_filt = fnames(~startsWith({fnames.name}, '._') & ~startsWith({fnames.name}, '.') & ~startsWith({fnames.name}, '..'));
connIDs_filename = sprintf('connIDs_BIDS_MPRAGE_60_89_batch_%d_HC', batch);
load(connIDs_filename, sprintf('connIDs_BIDS_MPRAGE_60_89_batch_%d_HC', batch));
connIDs_varname = sprintf('connIDs_BIDS_MPRAGE_60_89_batch_%d_HC', batch);
connIDs = eval(connIDs_varname); % Get connIDs for the current batch

    % Filter file names based on connIDs
    fnames_matching = fnames_filt(connIDs);
    
    % Initialize a container for time series
    tseries = cell(length(fnames_matching), 1);
    
    % Loop through each subject
    for s = 1:length(fnames_matching)
        % Load the data for the current subject
        load(fullfile(data_path_batch, fnames_matching(s).name));
        
        % Process the data for 1000 Schaefer nodes
        data_sch1000 = zeros(1000, 197); % Initialize array for efficiency
        for i = 1:1000
            if batch==1
                data_sch1000(i, :) = data{1, 945 + i}(1:197)'; % Extract and transpose
            else
                data_sch1000(i, :) = data{1, 567 + i}(1:197)'; % Extract and transpose
            end
            tseries{s} = data_sch1000; % Store the result
        end
    end
    
    % Save the time series data for the current batch
    input_path_batch=fullfile(input_path,tseries_path{1,batch});
    cd(input_path_batch);
    save_filename = sprintf('tseries_ADNI3_HC_MPRAGE_batch%d_sch1000_matching_QC.mat', batch);
    save(save_filename, 'tseries');
    
    % Clear variables except for those needed in the next iteration
    clearvars -except data_path input_path batch_path tseries_path batch
end


