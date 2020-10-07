function dth_SVM_pilot_2(subject)
%DTH_SVM_PILOT_2 Performs the distance-to-hyperplane analysis using SVM on
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

%data and results
data_dir = fullfile('/scratch/agnek95/PDM/DATA/DATA_PILOT_2/',subname);
results_dir = fullfile('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/',subname);
addpath(genpath(data_dir));
addpath(results_dir);
%% Prepare data
%load data
load(fullfile(data_dir,'timelock')); %eeg
load(fullfile(data_dir,'preprocessed_behavioural_data'));

%only keep the trials with a positive RT & correct response
triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); 
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); 

%% Define the required variables
numConditions = 60;
[numTrials, ~] = min_number_trials(triggers, numConditions); 
numTimepoints = size(timelock_data,3);
numPermutations=1; 
last_artificial_sample = 30;
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
    data = create_data_matrix(numConditions,triggers,numTrials,timelock_data);

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
        labels_train = [ones(num_conditions_batch,1); 2*ones(num_conditions_batch,1)]; %one label for each pseudotrial
        labels_test = labels_train; % we have the same size of training and testing data bcs we don't care about accuracy
        
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
        decisionValues(perm,conditions_art_batch_1,t) = decision_values_run2(1:numel(decision_values_run2)/2);
        decisionValues(perm,conditions_art_batch_2,t) = decision_values_run1(1:numel(decision_values_run1)/2);
        decisionValues(perm,conditions_nat_batch_1,t) = decision_values_run2(1:numel(decision_values_run2)/2);
        decisionValues(perm,conditions_nat_batch_1,t) = decision_values_run1(1:numel(decision_values_run1)/2);              
    end
    toc
end


%% Save the decision values
decisionValuesAvg = squeeze(mean(decisionValues,1)); %avg over permutations
save(fullfile(results_dir,'decisionValues'),'decisionValuesAvg');

%% Get the average (over trials) reaction time for each condition
RT_per_condition = NaN(numConditions,1);
RT_correct = behav.RT(behav.RT > 0 & behav.points == 1);

for c = 1:numConditions
    RT_per_condition(c) = mean(RT_correct(triggers==c));
end

save(fullfile(results_dir,'RTs_correct_trials'),'RT_per_condition');

end
   
    
