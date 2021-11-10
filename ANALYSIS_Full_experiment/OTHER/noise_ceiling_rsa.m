function noise_ceiling_rsa(subjects,conditions) 
%NOISE_CEILING Calculate the noise ceiling for the representational similarity analysis
%of object decoding (categorization task). 
%
%Input: subject IDs (eg., 1:13), conditions ('all','artificial' or
%'natural')
%
%Output: 
% 1)Three 1xP noise ceiling vectors: all scenes, natural scenes and artificial scenes
% (P = # timepoints)
%
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load existing RDMs or construct them from individual ones
if strcmp(conditions,'all')
    conds = 1:60;
elseif strcmp(conditions,'artificial')
    conds = 1:30;
elseif strcmp(conditions,'natural')
    conds = 31:60;
else
    error('Wrong conditions specified');
end
numConditions = numel(conds);
numTimepoints = 200;
numSubjects = numel(subjects);
rdm_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints); 

%% Loop: collect results from all subjects 
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'split_within_task_rdm_pearson_categorization.mat'),'rdm_avg');           
    rdm_all_subjects(subject,:,:,:) = squeeze(mean(rdm_avg(:,conds,conds,:),1)); %average over the RDMs of each half of trials      
end   

%% Remove all subject nans 
rdm_subjects = rdm_all_subjects(subjects,:,:,:);
noise_ceiling = NaN(numSubjects,numTimepoints);

%% Average over N-1 subjects and correlate with the left-out subject 
for s = 1:numSubjects
    leftout_subject = squeeze(rdm_subjects(s,:,:,:));
    other_subjects = squeeze(nanmean(rdm_subjects(1:end~=s,:,:,:),1)); %avg over remaining subjects
   
    %take upper triangular & flatten
    leftout_subject(isnan(leftout_subject)) = 0;
    other_subjects(isnan(other_subjects)) = 0;
    leftout_flattened_cell = arrayfun(@(x) squareform(leftout_subject(:,:,x)+(leftout_subject(:,:,x))'),...
        1:numTimepoints,'UniformOutput',false);
    leftout_flattened = reshape(cell2mat(leftout_flattened_cell),[],numTimepoints);
    other_flattened_cell = arrayfun(@(x) squareform(other_subjects(:,:,x)+(other_subjects(:,:,x))'),...
        1:numTimepoints,'UniformOutput',false);
    other_flattened = reshape(cell2mat(other_flattened_cell),[],numTimepoints);   
    for tp = 1:numTimepoints
        noise_ceiling(s,tp) = corr(leftout_flattened(:,tp),other_flattened(:,tp));
    end
end

%% Average over iterations
noise_ceiling_lower_bound = mean(noise_ceiling,1); %avg over subjects
plot(noise_ceiling_lower_bound);

% Save matrices
save(fullfile(results_avg_dir,sprintf('noise_ceiling_rsa_subjects_%d_%d_%s_scenes',subjects(1),subjects(end),conditions)),'noise_ceiling_lower_bound');

end
    
