function  pseudotrials_SVM_artificial_vs_natural_decoding_PDM(subject,task) 
%PSEUDOTRIALS_SVM_ARTIFICIAL_VS_NATURAL_DECODING_PDM_FULL_EXPERIMENT Perform category decoding (artificial vs natural) using the SVM classifier. 
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

%subject and task name string
subname = get_subject_name(subject);
task_name = get_task_name(task); 

%check if there's a directory for that subject, otherwise create one
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
if ~isfolder(fullfile(results_dir,subname))
    mkdir(results_dir,subname);
end
 
%% Prepare data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/%s/',subname);
load(fullfile(data_dir,sprintf('timelock_%s',task_name))); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)));

%only keep the trials with a positive RT & correct response
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

%% Define the required variables
numConditions = 60;
num_categories = 2; %categories to decode
num_conditions_per_category = numConditions/num_categories;
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=100; 

%minimum number of trials per scene
[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); 

%Preallocate 
decodingAccuracy=NaN(numPermutations,numTimepoints);

%% Decoding
rng('shuffle');
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); %returns: numConditions x numTrials x numElectrodes x numTimepoints
      
    disp('Split into artificial and natural');
    data_artificial = data(1:num_conditions_per_category,:,:,:);
    data_natural = data(num_conditions_per_category+1:end,:,:,:);
    
    %%% Fix: 1) permute the conditions 2) take half of the conditions and
    %%% mix the trials from the different scenes -> training set 
        
    
    
    
    disp('Permute the conditions (scenes)');
    conditions_order = randperm(num_conditions_per_category)';
    data_artificial_avg = data_artificial_avg(conditions_order,:,:);
    data_natural_avg = data_natural_avg(conditions_order,:,:);
    
    disp('Training set: Reshape by taking the trials from half of the conditions for each category');
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
    
    %%% Fix: 3) testing set is a vector of 30 conditions that were not used
    %%% in the training
%     disp('Testing set: Average over trials');
%     data_artificial_avg = squeeze(mean(data_artificial,2));
%     data_natural_avg = squeeze(mean(data_natural,2));
%     data_testing_both = NaN([num_categories,size(data_artificial_avg)]);
%     data_testing_both(1,:,:,:) = data_artificial_avg;
%     data_testing_both(2,:,:,:) = data_natural_avg;

%     testing_conditions = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins;

    for t = 1:numTimepoints 
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,:,:,t)); squeeze(bins(2,:,:,t))];  %train on half of the bins pseudotrials
        testing_data  = [squeeze(data_testing_both(1,:,:,t)); squeeze(data_testing_both(2,:,:,t))]; %test on all conditions 
        testing_data  = [squeeze(data_both_categories(1,testing_conditions,:,t)); squeeze(data_both_categories(2,testing_conditions,:,t))];
                 
        labels_train  = [ones(numBins-num_bins_testing,1);2*ones(numBins-num_bins_testing,1)]; %one label for each pseudotrial
        labels_test   = [ones(numel(testing_conditions),1);2*ones(numel(testing_conditions),1)];   
        
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
save(fullfile(results_dir,subname,sprintf('svm_artificial_vs_natural_decoding_accuracy_%s.mat',task_name)),'decodingAccuracy_avg');

