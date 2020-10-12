function ritchie_trial_level_dth_SVM_pilot_2(subject)
%RITCHIE_TRIAL_LEVEL_DTH_SVM_PILOT_2 Performs the distance-to-hyperplane analysis using SVM on
%a balanced dataset. 
%
%Input: subject ID (integer)
%
%Output: 
%   - NxP matrix of decision values, where N is the number of conditions
%   and P is the number of timepoints.
%   - Nx1 vector of RTs. 
%   - Nx1 vector of minimum trial #s. 
 
%% Add paths
%toolboxes and helper functions
addpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS/OTHER');
addpath(genpath('/home/agnek95/SMST/PDM/ANALYSIS/'));
addpath('/home/agnek95/OR/TOOLBOX/MVNN/MEG_SVM_decoding_MVNN'); %MVNN toolbox
addpath(genpath('/home/agnek95/OR/ANALYSIS/DECODING/libsvm')); %libsvm toolbox
addpath('/home/agnek95/OR/ANALYSIS/DECODING'); %MVNN function
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

subname = get_subject_name(subject);

%% Prepare data
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/ritchie_subject_%s',subname);
addpath(genpath(data_dir));
load(fullfile(data_dir,sprintf('s%s_PCA_S1_50hz.mat',subname))); %eeg
 
%take only the needed data 
triggers = data.TrialList(:,1);
category = data.TrialList(:,3);
task = data.TrialList(:,4);
response = data.TrialList(:,6);
RT = data.TrialList(:,7);
trials_final = [];
for t = 1:size(data.TrialList,1)
    if  category(t) == response(t) && task(t) == 1 %only take the active task and the correct trials 
        trials_final = [trials_final;t];
    end
end
triggers_final = triggers(trials_final);
data_final = permute(data.class_dat(trials_final,:,:),[1 3 2]);
RT_final = RT(trials_final);

%% Define the required variables
numConditions = 24;
[numTrials, ~] = min_number_trials(triggers_final, numConditions); 
numTimepoints = size(data_final,3);
numPermutations=1; 
last_artificial_sample = numConditions/2;%double check
numConditionsPerCategory = numConditions/2; 

%split trials into training and testing
trials_training = 1:numTrials*numConditions/4;
trials_testing = numTrials*(numConditions/4)+1:numTrials*numConditions/2;

%Preallocate 
decisionValues_artificial=NaN(numPermutations,numel(trials_training)*2,numTimepoints);
decisionValues_natural=NaN(numPermutations,numel(trials_training)*2,numTimepoints);

%% Running the MVPA
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,triggers_final,numTrials,data_final);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); %numConditions x numTrials x numElectrodes x numTimepoints

    disp('Split into artificial and natural');
    size_data = size(data);
    
    %permute the rows so that during reshaping we have the trials of all
    %conditions in the right order
    data_artificial = permute(data(1:last_artificial_sample,:,:,:),[2,1,3,4]);
    data_natural = permute(data(last_artificial_sample+1:end,:,:,:),[2,1,3,4]);

    %reshape to have the trials and conditions dimensions in one
    data_artificial_reshaped = reshape(data_artificial,[numTrials*numConditionsPerCategory,size_data(3:4)]);
    data_natural_reshaped = reshape(data_natural,[numTrials*numConditionsPerCategory,size_data(3:4)]);
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [data_artificial_reshaped(trials_training,:,t); data_natural_reshaped(trials_training,:,t)]; 
        testing_data  = [data_artificial_reshaped(trials_testing,:,t); data_natural_reshaped(trials_testing,:,t)];
        labels_train  = [ones(numel(trials_training),1); 2*ones(numel(trials_training),1)]; %one label for each pseudotrial
        labels_test   = labels_train; % we have the same size of training and testing data bcs we don't care about accuracy
        
        disp('Train the SVM: run 1');
        train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
        model=svmtrain_01(labels_train,training_data,train_param_str); 
        
        disp('Test the SVM: run 1');
        [~, ~, decision_values_run1] = svmpredict(labels_test,testing_data,model);  %for 60 conditions, you get a 30x1 vector of decision_values
        
        disp('Train the SVM: run 2');
        train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
        model=svmtrain_01(labels_test,testing_data,train_param_str); 
        
        disp('Test the SVM: run 2');
        [~, ~, decision_values_run2] = svmpredict(labels_train,training_data,model);  %for 60 conditions, you get a 30x1 vector of decision_values
        
        disp('Putting the decision values into the big matrix');
        decisionValues_artificial(perm,trials_training,t) = abs(decision_values_run2(1:numel(decision_values_run2)/2));
        decisionValues_artificial(perm,trials_testing,t) = abs(decision_values_run1(1:numel(decision_values_run1)/2));
        decisionValues_natural(perm,trials_training,t) = abs(decision_values_run2(1:numel(decision_values_run2)/2));
        decisionValues_natural(perm,trials_testing,t) = abs(decision_values_run1(1:numel(decision_values_run1)/2));              
    end
    toc
end

%% Rework the decision values
decisionValuesAll = cat(2,decisionValues_artificial,decisionValues_natural);
decisionValuesAvg = squeeze(mean(decisionValuesAll,1)); %avg over permutations

%% Get the average (over trials) reaction time for each condition
RT_per_condition = NaN(numConditions,1);
dv_per_condition = NaN(numConditions,numTimepoints);

for c = 1:numConditions
    RT_per_condition(c) = mean(RT_final(triggers_final==c));
    if c > 1
        trials = (numTrials*c-1)+1:numTrials*c;
    else 
        trials = 1:numTrials;
    end
    dv_per_condition(c,:) = mean(decisionValuesAvg(trials,:),1);
end

%% Save decision values and RTs
results_dir = fullfile('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/',subname);
save(fullfile(results_dir,'ritchie_trial_level_decisionValues'),'dv_per_condition');
save(fullfile(results_dir,'ritchie_trial_level_RTs_correct_trials'),'RT_per_condition');
end
   
    
