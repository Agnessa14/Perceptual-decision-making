function time_generalization_artificial_vs_natural_cross_removed_scene(subject) 
%TIME_GENERALIZATION_ARTIFICIAL_VS_NATURAL_CROSS_REMOVED_SCENE Perform category decoding
%(artificial vs natural) using the SVM classifier,
%in a time-generalized manner (trained & tested on all timepoints),
%across tasks. For the participants who had to have a scene removed due to
%not enough trials.
%
%Input: subject ID
%
%Output: PxP vector of accuracies in %, where P is the number of timepoints. 
%
%Author: Agnessa Karapetian, 2021
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
numConditionsAll = 60;
num_categories = 2; %categories to decode
[~, trials_per_condition_cat] = min_number_trials(timelock_triggers_categorization, numConditionsAll); %minimum number of trials per scene
[~, trials_per_condition_dis] = min_number_trials(timelock_triggers_distraction, numConditionsAll); %minimum number of trials per scene

removed_condition = find(trials_per_condition_cat==min(trials_per_condition_cat)); %the removed condition is always in the cat. task
low_minnumtrials = min(trials_per_condition_cat);
numTrials = min([trials_per_condition_cat(trials_per_condition_cat>low_minnumtrials),trials_per_condition_dis]);

numTimepoints = numResampledTps; %number of timepoints
numPermutations=100; 

%exclude trials from removed scene 
included_conditions = trials_per_condition_cat>=numTrials;
if removed_condition<=30
    num_conditions_artificial = 29;
    num_conditions_natural = 30; 
else
    num_conditions_artificial = 30;
    num_conditions_natural = 29; 
end

%Preallocate 
decodingAccuracy=NaN(numPermutations,numTimepoints,numTimepoints);

%% Decoding
rng('shuffle');
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data_categorization = create_data_matrix(numConditionsAll,timelock_triggers_categorization,numTrials,downsampled_timelock_data_categorization);
    data_distraction = create_data_matrix(numConditionsAll,timelock_triggers_distraction,numTrials,downsampled_timelock_data_distraction);
    data_categorization = data_categorization(included_conditions,:,:,:);
    
    disp('Performing MVNN');
    data_categorization = multivariate_noise_normalization(data_categorization); 
    data_distraction = multivariate_noise_normalization(data_distraction); 
      
    disp('Split into artificial and natural');
    data_artificial_cat = data_categorization(1:num_conditions_artificial,:,:,:);
    data_natural_cat = data_categorization(num_conditions_natural+1:end,:,:,:);
    data_artificial_dis = data_distraction(1:num_conditions_artificial,:,:,:);
    data_natural_dis = data_distraction(num_conditions_natural+1:end,:,:,:);
    
    disp('Average over trials');
    data_artificial_cat_avg = squeeze(mean(data_artificial_cat,2));
    data_natural_cat_avg = squeeze(mean(data_natural_cat,2));
    data_artificial_dis_avg = squeeze(mean(data_artificial_dis,2));
    data_natural_dis_avg = squeeze(mean(data_natural_dis,2));
    
    disp('Permute the conditions (scenes)');
    conditions_order_artificial = randperm(num_conditions_artificial)';
    conditions_order_natural = randperm(num_conditions_natural);
    data_artificial_cat_avg = data_artificial_cat_avg(conditions_order_artificial,:,:);
    data_natural_cat_avg = data_natural_cat_avg(conditions_order_natural,:,:);
    data_artificial_dis_avg = data_artificial_dis_avg(conditions_order_artificial,:,:);
    data_natural_dis_avg = data_natural_dis_avg(conditions_order_natural,:,:);
    
    disp('Insert NaN instead of the removed scene'); %remove the scene from both tasks for balance
    size_data_2 = size(data_artificial_cat_avg);
    if removed_condition<=30
        data_artificial_cat_avg = [data_artificial_cat_avg(1:removed_condition-1,:,:);
            NaN([1, size_data_2(2:3)]);data_artificial_cat_avg(removed_condition:end,:,:)];
        data_artificial_dis_avg = [data_artificial_dis_avg(1:removed_condition-1,:,:);
            NaN([1, size_data_2(2:3)]);data_artificial_dis_avg(removed_condition:end,:,:)];
    else
        removed_condition_nat = removed_condition-30;
        data_natural_cat_avg = [data_natural_cat_avg(1:removed_condition_nat-1,:,:);
            NaN([1, size_data_2(2:3)]);data_natural_cat_avg(removed_condition_nat:end,:,:)];
        data_natural_dis_avg = [data_natural_dis_avg(1:removed_condition_nat-1,:,:);
            NaN([1, size_data_2(2:3)]);data_natural_dis_avg(removed_condition_nat:end,:,:)];
    end
    
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
%     [bins_dis,numBins_dis] = create_pseudotrials(numScenesPerBin,data_both_categories_dis);

    num_bins_testing = 3;  
    if removed_condition<=30
        testing_conditions_artificial = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins_cat-1;
        testing_conditions_natural = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins_cat;
    else
        testing_conditions_artificial = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins_cat;
        testing_conditions_natural = (numScenesPerBin*num_bins_testing)+1:numScenesPerBin*numBins_cat-1;
    end
    for tp1 = 1:numTimepoints 
        disp('Split into training and testing');
        training_data = [squeeze(bins_cat(1,1:end-num_bins_testing,:,tp1)); squeeze(bins_cat(2,1:end-num_bins_testing,:,tp1))]; %train on half of the bins                 
        labels_train  = [ones(numBins_cat-num_bins_testing,1);2*ones(numBins_cat-num_bins_testing,1)]; %one label for each pseudotrial
        
        disp('Train the SVM');
        train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
        model=svmtrain_01(labels_train,training_data,train_param_str); 
        
        for tp2 = 1:numTimepoints
            disp('Test the SVM');
            testing_data  = [squeeze(data_both_categories_dis(1,testing_conditions_artificial,:,tp2)); squeeze(data_both_categories_dis(2,testing_conditions_natural,:,tp2))];
            labels_test   = [ones(numel(testing_conditions_artificial),1);2*ones(numel(testing_conditions_natural),1)];   
            [~, accuracy, ~] = svmpredict(labels_test,testing_data,model);
            decodingAccuracy(perm,tp1,tp2)=accuracy(1);         
        end
    end   
    toc
end

%% Save the decoding accuracy
timeg_decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
save(fullfile(results_dir,subname,'3PT_time_gen_svm_artificial_vs_natural_decoding_accuracy_cross_task.mat'),'timeg_decodingAccuracy_avg');

