function pearson_object_decoding_removed_scene_sl(subject,task)
%PEARSON_OBJECT_DECODING_REMOVED_SCENE Perform object decoding (average of pairwise object decoding) using 1-Pearson's correlation.
%For subjects with a removed scene due to not enough trials.
%
%Input: subject ID, task (1=categorization, 2=fixation)
%
%Output: NxNxP vector of correlations, where N is the number of conditions and
%P is the number of timepoints.
%
%Author: Agnessa Karapetian, 2021
%

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

load(fullfile(data_dir,sprintf('timelock_%s',task_name))); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)));

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
numChannels = 63;
chanIdx = 1:numChannels;

%exclude trials from removed scene
included_conditions = find(trials_per_condition>=numTrials);
numConditionsIncluded = numel(included_conditions);

%Preallocate
numTrialsPerBin = 5;
numPseudotrials = round(numTrials/numTrialsPerBin);
rdm=NaN(numPermutations,numPseudotrials,numConditionsIncluded,numConditionsIncluded,numTimepoints,numChannels);

%% Decoding
for perm = 1:numPermutations
    tic
    disp('Creating the data matrix');
    data = create_data_matrix(numConditionsAll,timelock_triggers,numTrials,timelock_data);
    data = data(included_conditions,:,:,:);
    
    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);
    
    disp('Binning data into pseudotrials');
    [pseudoTrials,numPseudotrials] = create_pseudotrials(numTrialsPerBin,data);
    
     % Modifing neighbourhood and pseudotrial matrix due to missing channels
    [neighbourhoods, missing_channel_ids,pseudoTrials] = correct_neighbourhood_map(timelock.label,pseudoTrials);
    
    %only get the upper diagonal
    for condA=1:numConditionsIncluded-1
        for condB = condA+1:numConditionsIncluded
            for tp = 1:numTimepoints
                for pt = 1:numPseudotrials
                    for iChan = chanIdx
                        if  ~ismember(iChan,missing_channel_ids) || isempty(missing_channel_ids)
                            neighbours = neighbourhoods(iChan , :);
                            neighbours = neighbours(~isnan(neighbours));
                            rdm(perm,pt,condA,condB,tp,iChan) = 1-corr(squeeze(pseudoTrials(condA,pt,neighbours,tp)),squeeze(pseudoTrials(condB,pt,neighbours,tp)),'type','Pearson');
                        end
                    end
                end
            end
        end
    end
    toc
end

%% Add NaN to the removed scene
rdm_avg = squeeze(mean(mean(rdm,1),2)); %average over permutations and pseudotrials
rdm_1 = [rdm_avg(1:removed_condition-1,:,:);NaN(1,numConditionsAll-1,200,numChannels);rdm_avg(removed_condition:end,:,:)];
rdm_2 = [rdm_1(:,1:removed_condition-1,:,:),NaN(numConditionsAll,1,200,numChannels),rdm_1(:,removed_condition:end,:,:)];
rdm_avg = rdm_2;

%% Save the representational dissimilarity matrix
save(fullfile(results_dir,subname,sprintf('rdm_pearson_%s.mat',task_name)),'rdm_avg');

