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
num_categories = 2; %categories to decode
num_conditions_per_category = numConditions/num_categories;

[numTrials, ~] = min_number_trials(triggers, numConditions); 
numTimepoints = size(timelock_data,3);
numPermutations=100; 

%Preallocate 
decisionValues_Artificial=NaN(numPermutations,num_conditions_per_category,numTimepoints);
decisionValues_Natural = NaN(numPermutations,num_conditions_per_category,numTimepoints);

%% Running the MVPA
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);

    disp('Split into artificial and natural');
    data_artificial = data(1:num_conditions_per_category,:,:,:);
    data_natural = data(num_conditions_per_category+1:end,:,:,:);
       
    disp('Average over trials');
    data_artificial_avg = squeeze(mean(data_artificial,2));
    data_natural_avg = squeeze(mean(data_natural,2));
   
    disp('Permute the conditions (scenes)');
    conditions_order = randperm(num_conditions_per_category)';
    data_artificial_avg = data_artificial_avg(conditions_order,:,:);
    data_natural_avg = data_natural_avg(conditions_order,:,:);
    
    disp('Put both categories into one matrix');
    data_both_categories = NaN([num_categories,size(data_artificial_avg)]);
    data_both_categories(1,:,:,:) = data_artificial_avg;
    data_both_categories(2,:,:,:) = data_natural_avg;
    
    disp('Split into bins of scenes');
    numScenesPerBin = 6;
    [bins,numBins] = create_pseudotrials(numScenesPerBin,data_both_categories);
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,:,:,t)); squeeze(bins(2,:,:,t))];   
        testing_data  = [squeeze(data_both_categories(1,:,:,t)); squeeze(data_both_categories(2,:,:,t))];
       
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
save(fullfile(results_dir,'pseudotrials_decisionValues'),'decisionValues_Avg');

%% Get the average (over trials) reaction time for each condition
RT_per_condition = NaN(numConditions,1);
RT_correct = behav.RT(behav.RT > 0 & behav.points == 1);

for c = 1:numConditions
    RT_per_condition(c) = mean(RT_correct(triggers==c));
end

save(fullfile(results_dir,'RTs_correct_trials'),'RT_per_condition');
end
   
    
