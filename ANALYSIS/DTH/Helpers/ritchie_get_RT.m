function ritchie_get_RT(subjects)

%RITCHIE_GET_RT Calculate the average RT for each object for the
%replication of the DTH analysis of Ritchie (2015, CompBiol) results.
%
%
%Input: subject ID(s) (e.g., 1515,1605)
%
%Returns a Nx1 vector of RTs, N being the number of conditions.
%

%% Add paths 
data_dir = '/scratch/agnek95/PDM/ritchie2015_data/';
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS/';
addpath(data_dir);

%% Load the structure containing RTs
for s = 1:numel(subjects)
    subname = num2str(subjects(s));
    load(fullfile(data_dir,sprintf('s%s_PCA_S1_5RTanalysis_ClassifyResults.mat',subname)));
    meanRT = squeeze(Results.ActivePassive(1).RTdistanceMatrix(1,:,1)); %take the mean across trials

    %% Save
    save(fullfile(results_dir,subname,'RT_active_mean'),'meanRT');
end
end