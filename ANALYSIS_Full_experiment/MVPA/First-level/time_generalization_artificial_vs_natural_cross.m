function time_generalization_artificial_vs_natural_cross(subject) 
%TIME_GENERALIZATION_ARTIFICIAL_VS_NATURAL_CROSS Perform category decoding
%(artificial vs natural) using the SVM classifier,
%in a time-generalized manner (trained & tested on all timepoints),
%across tasks.
%
%Input: subject ID
%
%Output: PxP vector of accuracies in %, where P is the number of timepoints. 
%
%Author: Agnessa Karapetian, 2021

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

%check if there's a directory for that subject, otherwise create one
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
if ~isfolder(fullfile(results_dir,subname))
    mkdir(results_dir,subname);
end
 
%% Prepare data
%Categorization

%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/%s/',subname);
load(fullfile(data_dir,'timelock_categorization')); %eeg
load(fullfile(data_dir,'preprocessed_behavioural_data_categorization'));

%only keep the trials with a positive RT & correct response
timelock_triggers_categorization = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data_categorization = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

clear behav
clear timelock

%Distraction
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/%s/',subname);
load(fullfile(data_dir,'timelock_fixation')); %eeg
load(fullfile(data_dir,'preprocessed_behavioural_data_fixation'));

%only keep the trials with a positive RT & correct response
timelock_triggers_distraction = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data_distraction = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

clear behav
clear timelock

%% Downsample from 5ms/timepoint to 20ms
numResampledTps = 50; 
numTimepointsOld = size(timelock_data_categorization,3);
steps_per_tp = numTimepointsOld/numResampledTps; %number of old timepoints to average over to obtain one new timepoint 
size_data_cat = size(timelock_data_categorization);
size_data_dis = size(timelock_data_distraction);
downsampled_timelock_data_categorization = NaN([size_data_cat(1:2), numResampledTps]);
downsampled_timelock_data_distraction = NaN([size_data_dis(1:2), numResampledTps]);
i = 1;
for tp=1:numResampledTps
    downsampled_timelock_data_categorization(:,:,tp) = squeeze(mean(timelock_data_categorization(:,:,i:i+steps_per_tp-1),3));
    downsampled_timelock_data_distraction(:,:,tp) = squeeze(mean(timelock_data_distraction(:,:,i:i+steps_per_tp-1),3));
    i = i+steps_per_tp;
end  

%% Define the required variables
numConditions = 60;
num_categories = 2; %categories to decode
num_conditions_per_category = numConditions/num_categories;
numTimepoints = numResampledTps; %number of timepoints
numPermutations=100; 

%minimum number of trials per scene
numTrials = min(min_number_trials(timelock_triggers_categorization, numConditions),...
    min_number_trials(timelock_triggers_distraction,numConditions)); %minimum number of trials per scene

%Preallocate 
decodingAccuracy_1=NaN(numPermutations,numTimepoints,numTimepoints);
decodingAccuracy_2=NaN(numPermutations,numTimepoints,numTimepoints);

%% Decoding
rng('shuffle');
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data_categorization = create_data_matrix(numConditions,timelock_triggers_categorization,numTrials,downsampled_timelock_data_categorization);
    data_distraction = create_data_matrix(numConditions,timelock_triggers_distraction,numTrials,downsampled_timelock_data_distraction);
    
    disp('Performing MVNN');
    data_categorization = multivariate_noise_normalization(data_categorization); 
    data_distraction = multivariate_noise_normalization(data_distraction); 
      
    disp('Split into artificial and natural');
    data_artificial_cat = data_categorization(1:num_conditions_per_category,:,:,:);
    data_natural_cat = data_categorization(num_conditions_per_category+1:end,:,:,:);
    data_artificial_dis = data_distraction(1:num_conditions_per_category,:,:,:);
    data_natural_dis = data_distraction(num_conditions_per_category+1:end,:,:,:);
    
    disp('Average over trials');
    data_artificial_cat_avg = squeeze(mean(data_artificial_cat,2));
    data_natural_cat_avg = squeeze(mean(data_natural_cat,2));
    data_artificial_dis_avg = squeeze(mean(data_artificial_dis,2));
    data_natural_dis_avg = squeeze(mean(data_natural_dis,2));
    
    disp('Permute the conditions (scenes)');
    conditions_order = randperm(num_conditions_per_category)';
    data_artificial_cat_avg = data_artificial_cat_avg(conditions_order,:,:);
    data_natural_cat_avg = data_natural_cat_avg(conditions_order,:,:);
    data_artificial_dis_avg = data_artificial_dis_avg(conditions_order,:,:);
    data_natural_dis_avg = data_natural_dis_avg(conditions_order,:,:);
    
    disp('Put both categories into one matrix');
    data_both_categories_cat = NaN([num_categories,size(data_artificial_cat_avg)]);
    data_both_categories_dis = NaN([num_categories,size(data_artificial_dis_avg)]);
    data_both_categories_cat(1,:,:,:) = data_artificial_cat_avg;
    data_both_categories_cat(2,:,:,:) = data_natural_cat_avg;
    data_both_categories_dis(1,:,:,:) = data_artificial_dis_avg;
    data_both_categories_dis(2,:,:,:) = data_natural_dis_avg;

    disp('Split into bins of scenes');
    numScenesPerBin = 5;
    [bins_cat,numBins_cat] = create_pseudotrials(numScenesPerBin,data_both_categories_cat);
    [bins_dis,numBins_dis] = create_pseudotrials(numScenesPerBin,data_both_categories_dis);

    num_bins_testing = 2;  
    testing_conditions = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins_cat;
    
    for tp1 = 1:numTimepoints 
        %train model 1: on categorization data
        training_data_1 = [squeeze(bins_cat(1,1:end-num_bins_testing,:,tp1)); squeeze(bins_cat(2,1:end-num_bins_testing,:,tp1))];                 
        labels_train  = [ones(numBins_cat-num_bins_testing,1);2*ones(numBins_cat-num_bins_testing,1)]; %one label for each pseudotrial
        train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
        model_1=svmtrain_01(labels_train,training_data_1,train_param_str); 
        
        %train model 2: on distraction data
        training_data_2 = [squeeze(bins_dis(1,1:end-num_bins_testing,:,tp1)); squeeze(bins_dis(2,1:end-num_bins_testing,:,tp1))];               
        labels_train  = [ones(numBins_dis-num_bins_testing,1);2*ones(numBins_dis-num_bins_testing,1)]; 
        train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; 
        model_2=svmtrain_01(labels_train,training_data_2,train_param_str); 
        
        for tp2 = 1:numTimepoints
            %test model 1: on distraction data
            testing_data_1  = [squeeze(data_both_categories_dis(1,testing_conditions,:,tp2)); squeeze(data_both_categories_dis(2,testing_conditions,:,tp2))];
            labels_test   = [ones(numel(testing_conditions),1);2*ones(numel(testing_conditions),1)];   
            [~, accuracy_1, ~] = svmpredict(labels_test,testing_data_1,model_1);
            decodingAccuracy_1(perm,tp1,tp2)=accuracy_1(1);         
            
            %test model 2: on categorization data
            testing_data_2  = [squeeze(data_both_categories_cat(1,testing_conditions,:,tp2)); squeeze(data_both_categories_cat(2,testing_conditions,:,tp2))];
            labels_test   = [ones(numel(testing_conditions),1);2*ones(numel(testing_conditions),1)];   
            [~, accuracy_2, ~] = svmpredict(labels_test,testing_data_2,model_2);
            decodingAccuracy_2(perm,tp1,tp2)=accuracy_2(1);  
        end
    end   
    toc
end

%% Average over both models 
da1 = squeeze(mean(decodingAccuracy_1,1)); %Avg over permutations
da2 = squeeze(mean(decodingAccuracy_2,1));
timeg_decodingAccuracy_avg = (da1+da2)/2;

%% Save the decoding accuracy
save(fullfile(results_dir,subname,'2_models_time_gen_svm_artificial_vs_natural_decoding_accuracy_cross_task.mat'),'timeg_decodingAccuracy_avg');

