function within_task_pearson_object_decoding_removed_scene(subject,task) 
%WITHIN_TASK_PEARSON_OBJECT_DECODING_REMOVED_SCENE Perform object decoding (average of pairwise object decoding) 
%using correlation distance (1-Pearson's correlation). Meant for
%within-task RSA, therefore the trials are split in half, resulting in 2
%RDMs per timepoint. For the subjects with a removed scene.
%
%Input: subject ID, task (1=categorization, 2=fixation)
%
%Output: 2xNxNxP vector of correlations, where N is the number of conditions and
%P is the number of timepoints. 
%
% Author: Agnessa Karapetian, 2021

%% Set-up prereqs
%add paths
addpath(genpath('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT'));
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
addpath('/home/agnek95/OR/TOOLBOX/MVNN/MEG_SVM_decoding_MVNN'); %MVNN toolbox
addpath(genpath('/home/agnek95/OR/ANALYSIS/DECODING/libsvm')); %libsvm toolbox
addpath('/home/agnek95/OR/ANALYSIS/DECODING'); %MVNN function
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

%subject and task name strings
subname = get_subject_name(subject);
task_name = get_task_name(task); 

%check if there's a directory for that subject, otherwise create one
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
if ~isfolder(fullfile(results_dir,subname))
    mkdir(results_dir,subname);
end

%% Prepare data
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/%s/',subname);
load(fullfile(data_dir,sprintf('timelock_%s',task_name)),'timelock'); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)),'behav');

%only keep the trials with a positive RT & correct response
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

%% Define the required variables
numConditionsAll = 60;
[~, trials_per_condition] = min_number_trials(timelock_triggers, numConditionsAll); %minimum number of trials per scene
removed_condition = find(trials_per_condition==min(trials_per_condition));
low_minnumtrials = min(trials_per_condition);
numTrials = min(trials_per_condition(trials_per_condition>low_minnumtrials));
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=1; 

%exclude trials from removed scene 
included_conditions = find(trials_per_condition>=numTrials);
numConditionsIncluded = numel(included_conditions);

%Preallocate 
numPseudotrials = 6;
numTrialsPerBin = round(numTrials/numPseudotrials);
rdm=NaN(numPermutations,2,numPseudotrials/2,numConditionsIncluded,numConditionsIncluded,numTimepoints);

%Split data in two halves
data = create_data_matrix(numConditionsAll,timelock_triggers,numTrials,timelock_data);
data = data(included_conditions,:,:,:);
data = multivariate_noise_normalization(data);
half_data_1 = data(:,1:round(numTrials/2),:,:);
half_data_2 = data(:,round(numTrials/2)+1:end,:,:);

%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Permuting the trials');
    half_data_1_perm = half_data_1(:,randperm(size(half_data_1,2)),:,:);
    half_data_2_perm = half_data_2(:,randperm(size(half_data_2,2)),:,:);
    
    disp('Binning data into pseudotrials');
    [pseudoTrials_1,numPseudotrials_1] = create_pseudotrials(numTrialsPerBin,half_data_1_perm);    
    [pseudoTrials_2,numPseudotrials_2] = create_pseudotrials(numTrialsPerBin,half_data_2_perm);    
    min_num_PTs = min([numPseudotrials_1 numPseudotrials_2]);
    
    %only get the upper diagonal
    for condA=1:numConditionsIncluded-1 %1:59
        for condB = condA+1:numConditionsIncluded %2:60
            for tp = 1:numTimepoints 
                for pt = 1:min_num_PTs
                    rdm(perm,1,pt,condA,condB,tp) = 1-corr(squeeze(pseudoTrials_1(condA,pt,:,tp)),squeeze(pseudoTrials_1(condB,pt,:,tp)),'type','Pearson');
                    rdm(perm,2,pt,condA,condB,tp) = 1-corr(squeeze(pseudoTrials_2(condA,pt,:,tp)),squeeze(pseudoTrials_2(condB,pt,:,tp)),'type','Pearson');
                end
            end 
        end 
    end 
    toc
end


%% Add NaN to the removed scene
rdm_avg = squeeze(nanmean(nanmean(rdm,1),3)); %average over permutations and pseudotrials
rdm_1 = [rdm_avg(:,1:removed_condition-1,:,:),NaN(2,1,numConditionsAll-1,numTimepoints),rdm_avg(:,removed_condition:end,:,:)];
rdm_2 = cat(3,rdm_1(:,:,1:removed_condition-1,:),NaN(2,numConditionsAll,1,200),rdm_1(:,:,removed_condition:end,:));
rdm_avg = rdm_2;

%% Save the representational dissimilarity matrix
save(fullfile(results_dir,subname,sprintf('split_within_task_rdm_pearson_%s.mat',task_name)),'rdm_avg');

