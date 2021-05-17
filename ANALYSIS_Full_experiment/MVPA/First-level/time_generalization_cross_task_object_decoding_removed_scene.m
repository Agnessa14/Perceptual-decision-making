function time_generalization_cross_task_object_decoding_removed_scene(subject) 
%TIME_GENERALIZATION_CROSS_TASK_OBJECT_DECODING_REMOVED_SCENE Perform time-generalized
%(on each pair of timepoints) object decoding (average of pairwise object decoding) 
% using the SVM classifier. For the participants with a missing condition. 
%
%Input: subject ID
%
%Output: NxNxPxP vector of accuracies in %, where N is the number of conditions and
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
numTimepoints = numResampledTps; 
numPermutations=50; 
numConditionsAll = 60;
[min_cat,trials_per_condition_cat] = min_number_trials(timelock_triggers_categorization, numConditionsAll); 
[min_dis,trials_per_condition_dis] = min_number_trials(timelock_triggers_distraction, numConditionsAll);
if min_cat<min_dis
    removed_condition = find(trials_per_condition_cat==min(trials_per_condition_cat));
    trials_per_condition = trials_per_condition_cat;
    low_minnumtrials = min_cat;
elseif min_dis<min_cat
    removed_condition = find(trials_per_condition_dis==min(trials_per_condition_dis));
    trials_per_condition = trials_per_condition_dis;
    low_minnumtrials = min_dis;
end    
numTrials = min(trials_per_condition(trials_per_condition>low_minnumtrials));

%exclude trials from removed scene 
included_conditions = find(trials_per_condition>=numTrials);
numConditionsIncluded = numel(included_conditions);

%Preallocate 
decodingAccuracy=NaN(numPermutations,numConditionsIncluded,numConditionsIncluded,numTimepoints,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrices');
    data_categorization = create_data_matrix(numConditionsAll,timelock_triggers_categorization,numTrials,downsampled_timelock_data_categorization);
    data_categorization = data_categorization(included_conditions,:,:,:);
    data_distraction = create_data_matrix(numConditionsAll,timelock_triggers_distraction,numTrials,downsampled_timelock_data_distraction);
    data_distraction = data_distraction(included_conditions,:,:,:); 

    disp('Performing MVNN');
    data_categorization = multivariate_noise_normalization(data_categorization); 
    data_distraction = multivariate_noise_normalization(data_distraction); 

    disp('Binning data into pseudotrials');
    numTrialsPerBin = round(numTrials/6);
    [pseudoTrials_categorization,numPTs_categorization] = create_pseudotrials(numTrialsPerBin,data_categorization);
    [pseudoTrials_distraction,~] = create_pseudotrials(numTrialsPerBin,data_distraction);
    
    %loop over conditions and timepoints
    for condA=1:numConditionsIncluded-1
        for condB = condA+1:numConditionsIncluded 
            for timePoint1 = 1:numTimepoints                
                %train the model
                training_data = [squeeze(pseudoTrials_categorization(condA,1:end-1,:,timePoint1)) ; squeeze(pseudoTrials_categorization(condB,1:end-1,:,timePoint1))]; %(numbins-1)x63x1 each
                labels_train=[ones(1,numPTs_categorization-1) 2*ones(1,numPTs_categorization-1)]; %one label for each pseudotrial
                train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
                model=svmtrain_01(labels_train',training_data,train_param_str); 

                for timePoint2 = 1:numTimepoints
                    disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                    ', timepoints ->',num2str(timePoint1), ', and ->', num2str(timePoint2)]);
                    %test the model
                    testing_data = [squeeze(pseudoTrials_distraction(condA,end,:,timePoint2))' ; squeeze(pseudoTrials_distraction(condB,end,:,timePoint2))']; %1x63x1 each
                    labels_test=[1 2]; 
                    [~, accuracy, ~,] = svmpredict(labels_test',testing_data,model);
                    decodingAccuracy(perm,condA,condB,timePoint1,timePoint2)=accuracy(1);  
                end
            end 
        end 
    end 
    toc
end

%% Add NaN to the removed scene
timeg_decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
DA_1 = [timeg_decodingAccuracy_avg(1:removed_condition-1,:,:,:);NaN(1,numConditionsAll-1,numTimepoints,numTimepoints);...
    timeg_decodingAccuracy_avg(removed_condition:end,:,:,:)];
DA_2 = [DA_1(:,1:removed_condition-1,:,:),NaN(numConditionsAll,1,numTimepoints,numTimepoints),...
    DA_1(:,removed_condition:end,:,:)];
timeg_decodingAccuracy_avg = DA_2;

%% Save the decoding accuracy
save(fullfile(results_dir,subname,sprintf('time_generalized_svm_object_decoding_crosstask.mat')),'timeg_decodingAccuracy_avg');

end