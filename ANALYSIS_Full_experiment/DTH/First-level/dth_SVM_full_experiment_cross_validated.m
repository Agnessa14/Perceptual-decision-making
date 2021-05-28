function dth_SVM_full_experiment_cross_validated(subject,task)
%DTH_SVM_PILOT_FULL_EXPERIMENT_CROSS_VALIDATED Performs the distance-to-hyperplane analysis using SVM on
%a balanced dataset. 
%
%Input: subject ID (integer)
%
%Output: 
%   - NxP matrix of decision values, where N is the number of conditions
%   and P is the number of timepoints.
%   - Nx1 vector of RTs. 
%   - Nx1 vector of minimum trial #s. 
%
%Author: Agnessa Karapetian, 2021

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

%data and results
data_dir = fullfile('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/',subname);
results_dir = fullfile('/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/',subname);
addpath(genpath(data_dir));
addpath(results_dir);

%% Prepare data
%load data
load(fullfile(data_dir,sprintf('timelock_%s',task_name))); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)));

%only keep the trials with a positive RT & correct response
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); 
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); 

%% Define the required variables
numConditions = 60;
num_categories = 2; %categories to decode
num_conditions_per_category = numConditions/num_categories;

[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); 
numTimepoints = size(timelock_data,3);
numPermutations=100; 

%Preallocate 
decisionValues_Artificial=NaN(numPermutations,num_conditions_per_category,numTimepoints);
decisionValues_Natural = NaN(numPermutations,num_conditions_per_category,numTimepoints);

%% Running the MVPA
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);

    disp('Split into artificial and natural');
    data_artificial = data(1:num_conditions_per_category,:,:,:);
    data_natural = data(num_conditions_per_category+1:end,:,:,:);
       
    disp('Split trials into two groups (training and testing) and average over trials');
    numTrials_training = round(size(data_artificial,2)/2);
    data_artificial_avg_training = squeeze(mean(data_artificial(:,1:numTrials_training,:,:),2));
    data_natural_avg_training = squeeze(mean(data_natural(:,1:numTrials_training,:,:),2));
    data_artificial_avg_testing = squeeze(mean(data_artificial(:,numTrials_training+1:end,:,:),2));
    data_natural_avg_testing = squeeze(mean(data_natural(:,numTrials_training+1:end,:,:),2));
   
    disp('Permute the conditions (scenes)');
    conditions_order = randperm(num_conditions_per_category)';
    data_artificial_avg_training = data_artificial_avg_training(conditions_order,:,:);
    data_natural_avg_training = data_natural_avg_training(conditions_order,:,:);
    
    disp('Put both categories into one matrix');
    data_both_categories_training = NaN([num_categories,size(data_artificial_avg_training)]);
    data_both_categories_training(1,:,:,:) = data_artificial_avg_training;
    data_both_categories_training(2,:,:,:) = data_natural_avg_training;
    
    data_both_categories_testing = NaN([num_categories,size(data_artificial_avg_testing)]);
    data_both_categories_testing(1,:,:,:) = data_artificial_avg_testing;
    data_both_categories_testing(2,:,:,:) = data_natural_avg_testing;
    
    disp('Split into bins of scenes');
    numScenesPerBin = 4;
    [bins,numBins] = create_pseudotrials(numScenesPerBin,data_both_categories_training);
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,:,:,t)); squeeze(bins(2,:,:,t))];   
        testing_data  = [squeeze(data_both_categories_testing(1,:,:,t)); squeeze(data_both_categories_testing(2,:,:,t))];
       
        labels_train  = [ones(numBins,1);2*ones(numBins,1)]; %one label for each pseudotrial
        labels_test   = [ones(num_conditions_per_category,1);2*ones(num_conditions_per_category,1)]; % we have the same size of training and testing data bcs we don't care about accuracy
        
        disp('Train the SVM');
        train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
        model=svmtrain_01(labels_train,training_data,train_param_str); 
        
        disp('Test the SVM');
        [~, ~, decision_values] = svmpredict(labels_test,testing_data,model);  %for 60 conditions, you get a 30x1 vector of decision_values
        
        disp('Putting the decision values into the big matrix');
        for c = 1:num_conditions_per_category
            condition = conditions_order(c);
            decisionValues_Artificial(perm,condition,t) = abs(decision_values(c));
            decisionValues_Natural(perm,condition,t) = abs(decision_values(c+num_conditions_per_category));
        end
    end
    toc
end

%% Save the decision values
decisionValues_Artificial_Avg = squeeze(mean(decisionValues_Artificial,1)); %avg over permutations
decisionValues_Natural_Avg = squeeze(mean(decisionValues_Natural,1)); %avg over permutations
decisionValues_Avg = [decisionValues_Artificial_Avg;decisionValues_Natural_Avg];
filename = sprintf('cross_validated_decisionValues_%s.mat',task_name);
save(fullfile(results_dir,filename),'decisionValues_Avg');

end
   
    
