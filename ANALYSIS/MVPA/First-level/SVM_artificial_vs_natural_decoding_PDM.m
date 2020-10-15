function SVM_artificial_vs_natural_decoding_PDM(subject) 
%SVM_ARTIFICIAL_VS_NATURAL_DECODING_PDM Perform category decoding (artificial vs natural) using the SVM classifier. 
%
%Input: subject ID
%
%Output: NxNxP vector of accuracies in %, where N is the number of conditions and
%P is the number of timepoints. 
%
%% Set-up prereqs
%add paths
addpath(genpath('/scratch/agnek95/PDM/DATA/DATA_PILOT_2'));
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS/'));
addpath(genpath('/home/agnek95/SMST/PDM/ANALYSIS/'));
addpath('/home/agnek95/OR/TOOLBOX/MVNN/MEG_SVM_decoding_MVNN'); %MVNN toolbox
addpath(genpath('/home/agnek95/OR/ANALYSIS/DECODING/libsvm')); %libsvm toolbox
addpath('/home/agnek95/OR/ANALYSIS/DECODING'); %MVNN function
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

%subject name string
subname = get_subject_name(subject);
 
%% Prepare data
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_PILOT_2/%s/',subname);
load(fullfile(data_dir,'timelock')); %eeg
load(fullfile(data_dir,'preprocessed_behavioural_data'));

%only keep the trials with a positive RT & correct response
% timelock_data = timelock;
triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

%% Define the required variables
numConditions = 60;
num_conditions_per_category = numConditions/2;
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=1; 

%define the conditions batches that go into the training and testing sets
% conditions_batch_1 = 1:num_conditions_per_category/2;
% conditions_batch_2 = (num_conditions_per_category/2)+1:num_conditions_per_category;

%minimum number of trials per scene
[numTrials, ~] = min_number_trials(triggers, numConditions); 

%Preallocate 
decodingAccuracy=NaN(numPermutations,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); %returns: numConditions x numTrials x numElectrodes x numTimepoints
      
    disp('Split into artificial and natural');
    data_artificial = data(1:num_conditions_per_category,:,:,:);
    data_natural = data(num_conditions_per_category+1:end,:,:,:);
       
    disp('Average over trials');
    data_artificial_avg = squeeze(mean(data_artificial,2));
    data_natural_avg = squeeze(mean(data_natural,2));

    for t = 1:numTimepoints 
        disp('Split into training and testing');
        training_data = [data_artificial_avg(1:end-1,:,t); data_natural_avg(1:end-1,:,t)]; 
        testing_data  = [data_artificial_avg(end,:,t); data_natural_avg(end,:,t)];

        labels_train  = [ones(num_conditions_per_category-1,1); 2*ones(num_conditions_per_category-1,1)]; %one label for each pseudotrial
        labels_test   = [1;2];   

        disp('Train the SVM');
        train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
        model=svmtrain_01(labels_train,training_data,train_param_str); 

        disp('Test the SVM');
        [~, accuracy, ~,] = svmpredict(labels_test,testing_data,model);
        decodingAccuracy(perm,t)=accuracy(1);                
    end   

    toc
end

%% Save the decoding accuracy
decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
save(sprintf('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/%s/trained_on_58_svm_artificial_vs_natural_decoding_accuracy',subname),'decodingAccuracy_avg');

