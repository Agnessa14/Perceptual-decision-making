function SVM_artificial_vs_natural_decoding_PDM_removed_scene(subject,task) 
%SVM_ARTIFICIAL_VS_NATURAL_DECODING_PDM_REMOVED_SCENE Perform category decoding (artificial vs natural) using the SVM classifier. 
%For the subjects who didn't have enough correct trials for a scene and had
%to have it removed. 
%
%Input: subject ID, task (1=categorization, 2=fixation)
%
%Output: NxNxP vector of accuracies in %, where N is the number of conditions and
%P is the number of timepoints. 
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

%% Prepare data
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/%s/',subname);
load(fullfile(data_dir,sprintf('timelock_%s',task_name))); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)));

%only keep the trials with a positive RT & correct response
% timelock_data = timelock;
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

%% Define the required variables
numConditionsAll = 60;
[~, trials_per_condition] = min_number_trials(timelock_triggers, numConditionsAll); %minimum number of trials per scene
removed_condition = find(trials_per_condition==min(trials_per_condition));
low_minnumtrials = min(trials_per_condition);
numTrials = min(trials_per_condition(trials_per_condition>low_minnumtrials));
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=100; 

%exclude trials from removed scene 
included_conditions = find(trials_per_condition>=numTrials);
numConditionsIncluded = numel(included_conditions);

%Preallocate 
decodingAccuracy=NaN(numPermutations,numConditionsIncluded,numConditionsIncluded,numTimepoints);
if removed_condition<=30
    num_conditions_artificial = 29;
    num_conditions_natural = 30; 
else
    num_conditions_artificial = 30;
    num_conditions_natural = 29; 
end

%% Decoding
rng('shuffle');
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditionsAll,timelock_triggers,numTrials,timelock_data);
    data = data(included_conditions,:,:,:); 

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); %returns: numConditions x numTrials x numElectrodes x numTimepoints

    disp('Split into artificial and natural');
    data_artificial = data(1:num_conditions_artificial,:,:,:);
    data_natural = data(num_conditions_natural+1:end,:,:,:);
       
    disp('Average over trials and scenes');
    data_artificial_avg = squeeze(mean(data_artificial,2));
    data_natural_avg = squeeze(mean(data_natural,2));
    
    disp('Permute the conditions (scenes)');
    conditions_order_artificial = randperm(num_conditions_artificial)';
    conditions_order_natural = randperm(num_conditions_natural)';
    data_artificial_avg = data_artificial_avg(conditions_order_artificial,:,:);
    data_natural_avg = data_natural_avg(conditions_order_natural,:,:);
    
    disp('Insert NaN instead of the removed scene');
    size_data_2 = size(data_artificial_avg);
    if removed_condition<=30
        data_artificial_avg = [data_artificial_avg(1:removed_condition-1,:,:);
            NaN([1, size_data_2(2:3)]);data_artificial_avg(removed_condition:end,:,:)];
    else
        removed_condition_nat = removed_condition-30;
        data_natural_avg = [data_natural_avg(1:removed_condition_nat-1,:,:);
            NaN([1, size_data_2(2:3)]);data_natural_avg(removed_condition_nat:end,:,:)];
    end
    disp('Put both categories into one matrix');
    data_both_categories = NaN([2,size(data_artificial_avg)]);
    data_both_categories(1,:,:,:) = data_artificial_avg;
    data_both_categories(2,:,:,:) = data_natural_avg;
    
    disp('Split into bins of scenes');
    numScenesPerBin = 5;
    [bins,numBins] = create_pseudotrials(numScenesPerBin,data_both_categories);
    
    num_bins_testing = 3;  
    if removed_condition<=30
        testing_conditions_artificial = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins-1;
        testing_conditions_natural = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins;
    else
        testing_conditions_artificial = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins;
        testing_conditions_natural = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins-1;
    end
    
    for t = 1:numTimepoints 
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,1:end-num_bins_testing,:,t)); ...
            squeeze(bins(2,1:end-num_bins_testing,:,t))]; %train on half of the bins
        testing_data  = [squeeze(data_both_categories(1,testing_conditions_artificial,:,t));...
            squeeze(data_both_categories(2,testing_conditions_natural,:,t))]; %test on the other half
                 
        labels_train  = [ones(numBins-num_bins_testing,1);2*ones(numBins-num_bins_testing,1)]; %one label for each pseudotrial
        labels_test   = [ones(numel(testing_conditions_artificial),1);2*ones(numel(testing_conditions_natural),1)]; 
        
        disp('Train the SVM');
        train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
        model=svmtrain_01(labels_train,training_data,train_param_str); 

        disp('Test the SVM');
        [~, accuracy, ~] = svmpredict(labels_test,testing_data,model);
        decodingAccuracy(perm,t)=accuracy(1); 
        
    end   

    toc
end

%% Save the decoding accuracy
decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
save(sprintf('/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/%s/ssvm_artificial_vs_natural_decoding_accuracy_%s.mat',subname,task_name),'decodingAccuracy_avg');

