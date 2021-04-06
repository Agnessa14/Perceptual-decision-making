function SVM_object_decoding_full_experiment_removed_scene(subject,task) 
%SVM_OBJECT_DECODING_FULL_EXPERIMENT_REMOVED_SCENE Perform object decoding (average of pairwise object decoding) using the SVM classifier. 
%Script only for subject 10 (and 15), who had only 3 correct trials for scene 47,
%and therefore needs to have a different number of trials.
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

%% Define the required variables
numConditionsAll = 60;
[~, trials_per_condition] = min_number_trials(timelock_triggers, numConditionsAll); %minimum number of trials per scene
removed_condition = find(trials_per_condition==min(trials_per_condition));
low_minnumtrials = min(trials_per_condition);
numTrials = min(trials_per_condition(trials_per_condition>low_minnumtrials));
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=100; 

%exclude trials from removed scene 
included_conditions = find(trials_per_condition>=numTrials);
numConditionsIncluded = numel(included_conditions);

%Preallocate 
decodingAccuracy=NaN(numPermutations,numConditionsIncluded,numConditionsIncluded,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditionsAll,timelock_triggers,numTrials,timelock_data);
    data = data(included_conditions,:,:,:); 
    
    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); 

    disp('Binning data into pseudotrials');
    numTrialsPerBin = round(numTrials/6);
    [pseudoTrials,numPTs] = create_pseudotrials(numTrialsPerBin,data);
   
    %only get the upper diagonal
    for condA=1:numConditionsIncluded-1 %1:59
        for condB = condA+1:numConditionsIncluded %2:60
            for timePoint = 1:numTimepoints 
                disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                    ', timepoint ->',num2str(timePoint)]);
                
                % L-1 pseudo trials go to training set, the Lth to testing set
                if numPTs == 2 %when only 1 PT for each condition
                    training_data=[squeeze(pseudoTrials(condA,1:end-1,:,timePoint))' ; squeeze(pseudoTrials(condB,1:end-1,:,timePoint))']; %(numbins-1)x63x1 each
                else
                    training_data=[squeeze(pseudoTrials(condA,1:end-1,:,timePoint)) ; squeeze(pseudoTrials(condB,1:end-1,:,timePoint))]; %(numbins-1)x63x1 each
                end
                
                testing_data=[squeeze(pseudoTrials(condA,end,:,timePoint))' ; squeeze(pseudoTrials(condB,end,:,timePoint))']; %1x63x1 each
         
                % class labels
                labels_train=[ones(1,numPTs-1) 2*ones(1,numPTs-1)]; %one label for each pseudotrial
                labels_test= [1 2]; 
                
                train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
                model=svmtrain_01(labels_train',training_data,train_param_str); 
                [~, accuracy, ~,] = svmpredict(labels_test',testing_data,model);
                decodingAccuracy(perm,condA,condB,timePoint)=accuracy(1);                
            end 
        end 
    end 
    toc
end

%% Add NaN to the removed scene
decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
DA_1 = [decodingAccuracy_avg(1:removed_condition-1,:,:);NaN(1,numConditionsAll-1,200);decodingAccuracy_avg(removed_condition:end,:,:)];
DA_2 = [DA_1(:,1:removed_condition-1,:),NaN(numConditionsAll,1,200),DA_1(:,removed_condition:end,:)];
decodingAccuracy_avg = DA_2;

%% Save the decoding accuracy
save(sprintf('/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/%s/svm_decoding_accuracy_%s.mat',subname,task_name),'decodingAccuracy_avg');

%% Save the number of trials used 
save(fullfile(results_dir,subname,sprintf('num_trials_included_%s.mat',task_name)),'numTrials');

