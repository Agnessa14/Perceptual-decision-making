function glm_eeg_rnn_rt(conditions)
%GLM_EEG_RNN_RT Perform the GLM analysis to determine how much variance in
%human RT is explained by EEG data and rCNN RTs. 
%
%Input: conditions ('artificial','natural' or 'both')
%
%Output: results from the GLM analysis.
%
%Performs the GLM analysis on the selected scenes, using EEG data and rCNN RTs as regressors.
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';

%% Load the EEG data (distances to the hyperplane for all subjects), rCNN RTs and human RTs 
load(fullfile(results_avg_dir,'distances_all_subjects.mat'),'distances_all_subjects');
load(fullfile(results_avg_dir,'02.11_2_rnn/Model_RDM_redone','correlation_RT_human_RNN_cross-validated.mat'),'data');
rCNN_RTs = data;
load(fullfile(results_avg_dir,'RT_all_subjects_5_35_categorization.mat'),'RTs');
human_RTs = RTs;

%Average over subjects
mean_EEG_distances = squeeze(mean(distances_all_subjects,1));
mean_human_RTs = squeeze(mean(human_RTs,1));

%Select the appropriate scenes
if strcmp(conditions,'artificial')
    conds = 1:30;
elseif strcmp(conditions,'natural')
    conds = 31:60;
elseif strcmp(conditions,'both')
    conds = 1:60;
end
selected_EEG_distances = mean_EEG_distances(conds,:);
selected_human_RTs = mean_human_RTs(conds);
selected_rCNN_RTs = rCNN_RTs(conds);

%% 
for roi = 1:3
    
    mean_distances_fmri = mean(distances_fmri(:,:,roi),2); 
    
    for sub = 1:30
    for t = 1:size(distances_all_subjects,3)
        
        mean_distances_eeg = distances_all_subjects(sub,:,t); 
        [~,~,~,~,r2_behav_eeg] = regress(mean_distances_eeg', [ones(size(mean_distances_fmri,1),1) mean_RTs']); 
        [~,~,~,~,r2_fmri_eeg] = regress(mean_distances_eeg', [ones(size(mean_distances_fmri,1),1) mean_distances_fmri]); 
        [~,~,~,~,r2_full] = regress(mean_distances_eeg', [ones(size(mean_distances_fmri,1),1), mean_distances_fmri, mean_RTs']);
        
        dth_shared(t,sub,roi) = r2_fmri_eeg(1)+r2_behav_eeg(1)-r2_full(1); 
        dth_behav_eeg(t,sub,roi) = r2_behav_eeg(1); 
    end 
    end
end