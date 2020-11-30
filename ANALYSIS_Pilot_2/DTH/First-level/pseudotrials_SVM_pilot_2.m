function pseudotrials_SVM_pilot_2(subject)
%PSEUDOTRIALS_DTH_SVM_PILOT_2 Performs the distance-to-hyperplane analysis using SVM on
%a balanced dataset. Instead of creating pseudoconditions out of scenes,
%the trials from across conditions are lumped into pseudotrials.
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
decisionValues = NaN(numPermutations,numConditions,numTimepoints);

%% Running the MVPA
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);

    disp('Split into artificial and natural');
    data_artificial = data(1:num_conditions_per_category,:,:,:); %in Ritchie's dataset, the natural conditions come first
    data_natural = data(num_conditions_per_category+1:end,:,:,:);
       
    disp('Training set: Reshape by taking the trials from all conditions for each category');
    size_data_category = size(data_artificial);
    data_artificial_reshaped = reshape(data_artificial,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    data_natural_reshaped = reshape(data_natural,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    
    disp('Permute the trials')
    data_artificial_permuted = data_artificial_reshaped(randperm(size(data_artificial_reshaped,1)),:,:);
    data_natural_permuted = data_natural_reshaped(randperm(size(data_natural_reshaped,1)),:,:);
       
    disp('Put both categories into one matrix');
    data_both_categories = NaN([num_categories,size(data_artificial_permuted)]);
    data_both_categories(1,:,:,:) = data_artificial_permuted;
    data_both_categories(2,:,:,:) = data_natural_permuted;
    
    disp('Split into pseudotrials');
    numTrialsPerBin = 20; %try different combinations of bins/numTrialsPerBin
    [bins,numBins] = create_pseudotrials(numTrialsPerBin,data_both_categories);
    
    disp('Testing set: Average over trials');
    data_artificial_avg = squeeze(mean(data_artificial,2));
    data_natural_avg = squeeze(mean(data_natural,2));
    data_testing_both = NaN([num_categories,size(data_artificial_avg)]);
    data_testing_both(1,:,:,:) = data_artificial_avg;
    data_testing_both(2,:,:,:) = data_natural_avg;
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,:,:,t)); squeeze(bins(2,:,:,t))];  %train on all pseudotrials
        testing_data  = [squeeze(data_testing_both(1,:,:,t)); squeeze(data_testing_both(2,:,:,t))]; %test on all conditions 
       
        labels_train  = [ones(numBins,1);2*ones(numBins,1)]; 
        labels_test   = [ones(num_conditions_per_category,1);2*ones(num_conditions_per_category,1)]; 
        
        disp('Train the SVM');
        train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
        model=svmtrain_01(labels_train,training_data,train_param_str); 
        
        disp('Test the SVM');
        [~, ~, decision_values] = svmpredict(labels_test,testing_data,model);  
        
        disp('Putting the decision values into the big matrix');
        decisionValues(perm,:,t) = abs(decision_values);
        
    end
    toc
end

%% Save the decision values
decisionValues_Avg = squeeze(mean(decisionValues,1));
results_dir = fullfile('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/',subname);
save(fullfile(results_dir,'true_pseudotrials_svm_decisionValues'),'decisionValues_Avg');

end
   
    
