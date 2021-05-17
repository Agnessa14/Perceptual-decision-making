function time_generalization_within_task_object_decoding_removed_scene(subject,task) 
%TIME_GENERALIZATION_WITHIN_TASK_OBJECT_DECODING_REMOVED_SCENE Perform time-generalized (on each pair of timepoints) object decoding (average of pairwise object decoding) using the SVM classifier. 
%For subjects with a removed scene due to not enough trials.
%
%Input: subject ID, task (1=categorization, 2=fixation)
%
%Output: NxNxPxP vector of accuracies in %, where N is the number of conditions and
%P is the number of timepoints. 
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

%subject and task name strings
subname = get_subject_name(subject);
task_name = get_task_name(task); 

%check if there's a directory for that subject, otherwise create one
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
if ~isfolder(fullfile(results_dir,subname))
    mkdir(results_dir,subname);
end

%% Prepare data
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/%s/',subname);
load(fullfile(data_dir,sprintf('timelock_%s',task_name))); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)));

%only keep the trials with a positive RT & correct response
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

%% Downsample from 5ms/timepoint to 20ms
numResampledTps = 50; 
numTimepointsOld = size(timelock_data,3);
steps_per_tp = numTimepointsOld/numResampledTps; %number of old timepoints to average over to obtain one new timepoint 
size_data = size(timelock_data);
downsampled_timelock_data = NaN([size_data(1:2), numResampledTps]);
i = 1;
for tp=1:numResampledTps
    downsampled_timelock_data(:,:,tp) = squeeze(mean(timelock_data(:,:,i:i+steps_per_tp-1),3));
    i = i+steps_per_tp;
end 

%% Define the required variables
numConditionsAll = 60;
[~,trials_per_condition] = min_number_trials(timelock_triggers, numConditionsAll); 
removed_condition = find(trials_per_condition==min(trials_per_condition));
low_minnumtrials = min(trials_per_condition);
numTrials = min(trials_per_condition(trials_per_condition>low_minnumtrials));

%exclude trials from removed scene 
included_conditions = find(trials_per_condition>=numTrials);
numConditionsIncluded = numel(included_conditions);

numTimepoints = numResampledTps; 
numPermutations=1; 

%Preallocate 
decodingAccuracy=NaN(numPermutations,numConditionsIncluded,numConditionsIncluded,numTimepoints,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditionsAll,timelock_triggers,numTrials,downsampled_timelock_data);
    data = data(included_conditions,:,:,:);
    
    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); 

    disp('Binning data into pseudotrials');
    numTrialsPerBin = round(numTrials/6);
    [pseudoTrials,numPTs] = create_pseudotrials(numTrialsPerBin,data);
   
    %only get the upper diagonal
    for condA=1:numConditionsIncluded-1 
        for condB = condA+1:numConditionsIncluded%1:numConditions 
            for timePoint1 = 1:numTimepoints 
                % L-1 pseudo trials go to training set, the Lth to testing set
                training_data=[squeeze(pseudoTrials(condA,1:end-1,:,timePoint1)) ; squeeze(pseudoTrials(condB,1:end-1,:,timePoint1))]; %(numbins-1)x63x1 each
                labels_train=[ones(1,numPTs-1) 2*ones(1,numPTs-1)]; %one label for each pseudotrial

                %train the model
                train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
                model=svmtrain_01(labels_train',training_data,train_param_str); 
                for timePoint2 = 1:numTimepoints
                    disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                    ', timepoints ->',num2str(timePoint1), ', and ->', num2str(timePoint1)]);
                    testing_data=[squeeze(pseudoTrials(condA,end,:,timePoint2))' ; squeeze(pseudoTrials(condB,end,:,timePoint2))']; %1x63x1 each
                    labels_test= [1 2]; 
                    
                    %test the model
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
save(fullfile(results_dir,subname,sprintf('time_generalized_svm_object_decoding_%s.mat',task_name)),'timeg_decodingAccuracy_avg');

end