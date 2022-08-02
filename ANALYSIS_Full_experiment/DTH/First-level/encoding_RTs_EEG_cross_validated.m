function encoding_RTs_EEG_cross_validated(subject,task)
%encoding_RTs_EEG_cross_validated
%Performs the regression RT analysis using ridge regression to map eeg
%features on to RTs 
%
%Input: subject ID (integer), task (1=categorization, 2=fixation)
%
%Author: Agnessa Karapetian, Johannes Singer 2022

%% Add paths
%toolboxes and helper functions
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
addpath('/home/agnek95/OR/TOOLBOX/MVNN/MEG_SVM_decoding_MVNN'); %MVNN toolbox
addpath(genpath('/home/agnek95/OR/ANALYSIS/DECODING/libsvm')); %libsvm toolbox
addpath('/home/agnek95/OR/ANALYSIS/DECODING'); %MVNN function
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

subname = get_subject_name(subject);
task_name = get_task_name(task);

%check if there's a directory for that subject, otherwise create one
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
if ~isfolder(fullfile(results_dir,subname))
    mkdir(results_dir,subname);
end

%data and results
data_dir = fullfile('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/',subname);
results_dir = fullfile(results_dir,subname);
addpath(genpath(data_dir));
addpath(results_dir);

%% Prepare data
%load data
load(fullfile(data_dir,sprintf('timelock_%s',task_name)),'timelock'); %eeg

% load the group level behavioral data here 
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
filename_RT = sprintf('RT_all_subjects_5_35_%s.mat',task_name);
load(fullfile(results_avg_dir,filename_RT),'RTs');

% I am not sure how the variable that the group-level RTs contains is called but I refer to it later as "RTs"
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)),'behav'); 

%only keep the trials with a positive RT & correct response
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); 
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); 

%% Define the required variables
numConditions = 60;

[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); 
numTimepoints = size(timelock_data,3);
numPermutations=100; 

%Preallocate 
encodingAccuracy=NaN(numPermutations,numTimepoints);

%% Running the MVPA
rng('shuffle');
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);
    
    disp('Split into pseudotrials'); 
    numTrialsPerBin = 5; %try different combinations of bins/numTrialsPerBin
    [bins,~] = create_pseudotrials(numTrialsPerBin,data);
    
    % get the indices to split the behavioral RTs into two random
    % halfs (15 subjects in train, 15 subjects in test) 
    itrain_RTs = randperm(size(RTs,1),size(RTs,1)/2);
    itest_RTs  = setdiff(1:size(RTs,1),itrain_RTs); 
    
    % split the behavioral RTs into train and test
    train_RTs = squeeze(nanmean(RTs(itrain_RTs,:))); %here using nanmean bcs some participants dont have RTs for some scenes 
    test_RTs = squeeze(nanmean(RTs(itest_RTs,:)));

    for t = 1:numTimepoints
        % need to average training data across pseudotrials such that the
        % output shape is 60xn_channels
        training_data = squeeze(mean(bins(:,1:end-1,:,t),2));  %train on all pseudotrials
        testing_data = squeeze(bins(:,end,:,t));  %test on all pseudotrials
             
        % regress the eeg patterns onto the RTs and obtain weights
        weights = ridge(train_RTs',training_data,0.01,0);
        % get predicted RTs
        y_hat = weights(1) + testing_data*weights(2:end);
        % get prediction accuracy -> correlation between test RTs and
        % predicted RTs 
        encodingAccuracy(perm,t) = corr(y_hat,test_RTs');
        
    end
    toc
end

%% Save the decision values + decoding accuracy
encodingAccuracy_avg = squeeze(mean(encodingAccuracy,1)); 
filename = 'cross_validated_regression_RTs';

save(fullfile(results_dir,sprintf('%s_encodingAccuracy_%s.mat',filename,task_name)),'encodingAccuracy_avg');

end
   