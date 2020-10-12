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
numPermutations=100; 
last_artificial_sample = numConditions/2;%double check
num_conditions_batch = numConditions/4; %each training and testing set (for both artificial and natural) has the quarter of all conditions

%Preallocate 
decisionValues=NaN(numPermutations,numConditions,numTimepoints);

%define the conditions batches that go into the training and testing sets
conditions_art_batch_1 = 1:numConditions/4;
conditions_art_batch_2 = numConditions/4+1:numConditions/2;
conditions_nat_batch_1 = 1:numConditions/4;
conditions_nat_batch_2 = numConditions/4+1:numConditions/2;

%% Running the MVPA
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,triggers_final,numTrials,data_final);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);

    disp('Split into artificial and natural');
    data_artificial = data(1:last_artificial_sample,:,:,:);
    data_natural = data(last_artificial_sample+1:end,:,:,:);
       
    disp('Average over trials');
    data_artificial_avg = squeeze(mean(data_artificial,2));
    data_natural_avg = squeeze(mean(data_natural,2));
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [data_artificial_avg(conditions_art_batch_1,:,t); data_natural_avg(conditions_nat_batch_1,:,t)]; 
        testing_data  = [data_artificial_avg(conditions_art_batch_2,:,t); data_natural_avg(conditions_nat_batch_2,:,t)];
        labels_train  = [ones(num_conditions_batch,1); 2*ones(num_conditions_batch,1)]; %one label for each pseudotrial
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
        decisionValues(perm,conditions_art_batch_1,t) = abs(decision_values_run2(1:numel(decision_values_run2)/2));
        decisionValues(perm,conditions_art_batch_2,t) = abs(decision_values_run1(1:numel(decision_values_run1)/2));
        decisionValues(perm,conditions_nat_batch_1+last_artificial_sample,t) = abs(decision_values_run2(1:numel(decision_values_run2)/2));
        decisionValues(perm,conditions_nat_batch_2+last_artificial_sample,t) = abs(decision_values_run1(1:numel(decision_values_run1)/2));              
    end
    toc
end

%% Save the decision values
decisionValuesAvg = squeeze(mean(decisionValues,1)); %avg over permutations
results_dir = fullfile('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/',subname);
save(fullfile(results_dir,'all_trials_ritchie_decisionValues'),'decisionValuesAvg');

%% Get the average (over trials) reaction time for each condition
RT_per_condition = NaN(numConditions,1);

for c = 1:numConditions
    RT_per_condition(c) = mean(RT_final(triggers_final==c));
end

save(fullfile(results_dir,'all_trials_ritchie_RTs_correct_trials'),'RT_per_condition');
end
   
    
