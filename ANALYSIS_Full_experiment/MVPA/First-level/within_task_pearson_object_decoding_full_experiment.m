function within_task_pearson_object_decoding_full_experiment(subject,task) 
%PEARSON_OBJECT_DECODING_FULL_EXPERIMENT Perform object decoding (average of pairwise object decoding) 
%using correlation distance (1-Pearson's correlation). Meant for
%within-task RSA, therefore the trials are split in half, resulting in 2
%RDMs per timepoint.
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
numConditions = 60;
[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); %minimum number of trials per scene
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=100; 

%Preallocate 
numPseudotrials = 6;
numTrialsPerBin = round(numTrials/numPseudotrials);
rdm=NaN(2,numPermutations,numPseudotrials,numConditions,numConditions,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); 

    disp('Binning data into pseudotrials');
    [pseudoTrials,numPseudotrials] = create_pseudotrials(numTrialsPerBin,data);    
    
    %only get the upper diagonal
    for h = 1:2 %each half of pseudotrials
        if h == 1
            half_pseudotr = pseudoTrials(:,1:numPseudotrials/2,:,:);
        else
            half_pseudotr = pseudoTrials(:,(numPseudotrials/2)+1:end,:,:);
        end
        for condA=1:numConditions-1 %1:59
            for condB = condA+1:numConditions %2:60
                for tp = 1:numTimepoints 
                    for pt = 1:numPseudotrials/2
                        rdm(perm,h,pt,condA,condB,tp) = 1-corr(squeeze(half_pseudotr(condA,pt,:,tp)),squeeze(half_pseudotr(condB,pt,:,tp)),'type','Pearson');
                        disp(rdm(perm,h,pt,condA,condB,tp));
                    end
                end 
            end 
        end 
    end
    toc
end

%% Save the representational dissimilarity matrix
rdm_avg = squeeze(mean(rdm,[1,3])); %average over permutations and pseudotrials
save(fullfile(results_dir,subname,sprintf('within_task_rdm_pearson_%s.mat',task_name)),'rdm_avg');

