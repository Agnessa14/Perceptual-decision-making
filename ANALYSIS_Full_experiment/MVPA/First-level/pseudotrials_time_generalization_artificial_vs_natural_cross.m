function pseudotrials_time_generalization_artificial_vs_natural_cross(subject) 
%PSEUDOTRIALS_TIME_GENERALIZATION_ARTIFICIAL_VS_NATURAL_CROSS Perform time-generalized (on each pair of timepoints) 
%artificial vs natural decoding on data from both tasks using the SVM classifier. 
%
%Input: subject ID
%
%Output: PxP vector of accuracies in %, where P is the number of timepoints. 
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
numTrials = min(min_number_trials(timelock_triggers_categorization, numConditions),...
    min_number_trials(timelock_triggers_distraction,numConditions)); %minimum number of trials per scene
num_categories = 2; %categories to decode
num_conditions_per_category = numConditions/num_categories;
numTimepoints = numResampledTps; %number of timepoints
numPermutations=100; 

%Preallocate 
decodingAccuracy=NaN(numPermutations,numConditions,numConditions,numTimepoints,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrices');
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
    
    disp('Reshape by taking the trials from half of the conditions for each category');
    size_data_category = size(data_artificial_cat);
    data_artificial_cat_reshaped = reshape(data_artificial_cat,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    data_natural_cat_reshaped = reshape(data_natural_cat,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    data_artificial_dis_reshaped = reshape(data_artificial_dis,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    data_natural_dis_reshaped = reshape(data_natural_dis,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    
    disp('Permute the trials')
    data_artificial_cat_permuted = data_artificial_cat_reshaped(randperm(size(data_artificial_cat_reshaped,1)),:,:);
    data_natural_cat_permuted = data_natural_cat_reshaped(randperm(size(data_natural_cat_reshaped,1)),:,:);
    data_artificial_dis_permuted = data_artificial_dis_reshaped(randperm(size(data_artificial_dis_reshaped,1)),:,:);
    data_natural_dis_permuted = data_natural_dis_reshaped(randperm(size(data_natural_dis_reshaped,1)),:,:);

    disp('Put both categories into one matrix');
    data_both_categories_cat = NaN([num_categories,size(data_artificial_cat_permuted)]);
    data_both_categories_cat(1,:,:,:) = data_artificial_cat_permuted;
    data_both_categories_cat(2,:,:,:) = data_natural_cat_permuted;
    
    data_both_categories_dis = NaN([num_categories,size(data_artificial_dis_permuted)]);
    data_both_categories_dis(1,:,:,:) = data_artificial_dis_permuted;
    data_both_categories_dis(2,:,:,:) = data_natural_dis_permuted;
    
    disp('Split into pseudotrials');
    numTrialsPerBin = 20; %try different combinations of bins/numTrialsPerBin
    testing_bins = 5;     
    [bins_cat,numBins_cat] = create_pseudotrials(numTrialsPerBin,data_both_categories_cat);
    [bins_dis,numBins_dis] = create_pseudotrials(numTrialsPerBin,data_both_categories_dis);
    
    for tp1 = 1:numTimepoints 
        disp('Split into training and testing');
        training_data = [squeeze(bins_cat(1,1:end-testing_bins,:,tp1)); squeeze(bins_dis(2,1:end-testing_bins,:,tp1))];  %train on all but one bin                  
        labels_train  = [ones(numBins_cat-testing_bins,1);2*ones(numBins_dis-testing_bins,1)]; %one label for each pseudotrial
        
        disp('Train the SVM');
        train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
        model=svmtrain_01(labels_train,training_data,train_param_str); 
        
        for tp2 = 1:numTimepoints
            disp('Test the SVM');
            testing_data  = [squeeze(bins_cat(1,end-testing_bins+1:end,:,tp2)); squeeze(bins_dis(2,end-testing_bins+1:end,:,tp2))];
            labels_test   = [ones(testing_bins,1);2*ones(testing_bins,1)];   
            [~, accuracy, ~] = svmpredict(labels_test,testing_data,model);
            decodingAccuracy(perm,tp1,tp2)=accuracy(1);         
        end
    end   
    toc
end

%% Save the decoding accuracy
timeg_decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
save(fullfile(results_dir,subname,'timegen_pseudotrials_svm_artificial_vs_natural_decoding_accuracy_cross_task.mat'),'timeg_decodingAccuracy_avg');

end
    
    
