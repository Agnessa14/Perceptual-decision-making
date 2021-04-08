function time_generalization_within_task_object_decoding_full_experiment(subject,task) 
%TIME_GENERALIZATION_WITHIN_TASK_OBJECT_DECODING_FULL_EXPERIMENT Perform time-generalized (on each pair of timepoints) object decoding (average of pairwise object decoding) using the SVM classifier. 
%
%Input: subject ID, task (1=categorization, 2=fixation)
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
numConditions = 60;
[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); %minimum number of trials per scene
numTimepoints = size(timelock_data,3); 
numPermutations=100; 

%Preallocate 
decodingAccuracy=NaN(numPermutations,numConditions,numConditions,numTimepoints,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_triggers,numTrials,timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); 

    disp('Binning data into pseudotrials');
    numTrialsPerBin = round(numTrials/6);
    [pseudoTrials,numPTs] = create_pseudotrials(numTrialsPerBin,data);
   
    %only get the upper diagonal
    for condA=1:numConditions-1 %1:59 max
        for condB = condA+1:numConditions %2:60 max
            for timePoint1 = 1:numTimepoints 
                for timePoint2 = 1:numTimepoints
                    disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                        ', timepoints ->',num2str(timePoint1), ', and ->', num2str(timePoint2)]);

                    % L-1 pseudo trials go to training set, the Lth to testing set
                    if numPTs == 2 %when only 1 PT for each condition
                        training_data=[squeeze(pseudoTrials(condA,1:end-1,:,timePoint1))' ; squeeze(pseudoTrials(condB,1:end-1,:,timePoint2))']; %(numbins-1)x63x1 each
                    else
                        training_data=[squeeze(pseudoTrials(condA,1:end-1,:,timePoint1)) ; squeeze(pseudoTrials(condB,1:end-1,:,timePoint2))]; %(numbins-1)x63x1 each
                    end

                    testing_data=[squeeze(pseudoTrials(condA,end,:,timePoint1))' ; squeeze(pseudoTrials(condB,end,:,timePoint2))']; %1x63x1 each

                    % class labels
                    labels_train=[ones(1,numPTs-1) 2*ones(1,numPTs-1)]; %one label for each pseudotrial
                    labels_test= [1 2]; 
                    
                    %train & test the model
                    train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
                    model=svmtrain_01(labels_train',training_data,train_param_str); 
                    [~, accuracy, ~,] = svmpredict(labels_test',testing_data,model);
                    decodingAccuracy(perm,condA,condB,timePoint1,timePoint2)=accuracy(1); 
                end
            end 
        end 
    end 
    toc
end

%% Save the decoding accuracy
timeg_decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
save(fullfile(results_dir,subname,sprintf('time_generalized_svm_decoding_accuracy_%s.mat',task_name)),'timeg_decodingAccuracy_avg');

end