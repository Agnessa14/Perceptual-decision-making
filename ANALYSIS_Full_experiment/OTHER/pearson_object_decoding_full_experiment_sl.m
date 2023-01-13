function pearson_object_decoding_full_experiment_sl(subject,task)
%PEARSON_OBJECT_DECODING_FULL_EXPERIMENT Perform object decoding (average of pairwise object decoding) using 1-Pearson's correlation.
%
%Input: subject ID, task (1=categorization, 2=fixation)
%
%Output: NxNxP vector of correlations, where N is the number of conditions and
%P is the number of timepoints.
%
% Author: Agnessa Karapetian, 2021, modified by Muthukumar Pandaram (2022)

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
times = -195:100:705;
time_2_idx = (times/5)+40;
numConditions = 60;
[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); %minimum number of trials per scene
numPermutations=100;
numChannels = 63;
chanIdx = 1:numChannels;

%Preallocate
numTrialsPerBin = 5;
numPseudotrials = round(numTrials/numTrialsPerBin);
% rdm=NaN(numPermutations,numPseudotrials,numConditions,numConditions,numel(times),numChannels);
rdm_1=NaN(round(numPermutations/2),numPseudotrials,numConditions,numConditions,numel(times),numChannels);
rdm_2=NaN(round(numPermutations/2),numPseudotrials,numConditions,numConditions,numel(times),numChannels);

%% Decoding
for perm = 1:numPermutations
    tic
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_triggers,numTrials,timelock_data);
    
    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);
    
    disp('Binning data into pseudotrials');
    [pseudoTrials,numPseudotrials] = create_pseudotrials(numTrialsPerBin,data);
    
     % Modifing neighbourhood and pseudotrial matrix due to missing channels
    [neighbourhoods, missing_channel_ids,pseudoTrials] = correct_neighbourhood_map(timelock.label,pseudoTrials);
    
    %only get the upper diagonal
    for condA=1:numConditions-1 %1:59
        for condB = condA+1:numConditions %2:60
            for tp = 1:numel(time_2_idx)
                t = time_2_idx(tp);
                for pt = 1:numPseudotrials
                    for iChan = chanIdx
                        if ~ismember(iChan,missing_channel_ids) || isempty(missing_channel_ids)
                            neighbours = neighbourhoods(iChan , :);
                            neighbours = neighbours(~isnan(neighbours));
%                             rdm(perm,pt,condA,condB,tp,iChan) = 1-corr(squeeze(pseudoTrials(condA,pt,neighbours,t)),squeeze(pseudoTrials(condB,pt,neighbours,t)),'type','Pearson');
                            if perm <= 50
                                rdm_1(perm,pt,condA,condB,tp,iChan) = 1-corr(squeeze(pseudoTrials(condA,pt,neighbours,t)),squeeze(pseudoTrials(condB,pt,neighbours,t)),'type','Pearson');
                            else
                                rdm_2(perm-50,pt,condA,condB,tp,iChan) = 1-corr(squeeze(pseudoTrials(condA,pt,neighbours,t)),squeeze(pseudoTrials(condB,pt,neighbours,t)),'type','Pearson');
                        
                            end
                        end    
                    end
                end
            end
        end
    end
    toc
end

%% Save the representational dissimilarity matrix
% rdm_avg = squeeze(mean(mean(rdm,1),2)); %average over permutations and pseudotrials
rdm_avg_both = NaN(2,numConditions,numConditions,numel(times),numChannels);
rdm_avg_both(1,:,:,:,:) = squeeze(mean(rdm_1,1:2));%average over permutations and pseudotrials
rdm_avg_both(2,:,:,:,:) = squeeze(mean(rdm_2,1:2));
rdm_avg = squeeze(mean(rdm_avg_both,1));
save(fullfile(results_dir,subname,sprintf('rdm_pearson_searchlight_%s.mat',task_name)),'rdm_avg');

