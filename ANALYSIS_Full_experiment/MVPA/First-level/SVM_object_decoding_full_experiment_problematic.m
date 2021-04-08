function SVM_object_decoding_full_experiment_problematic(subject,task) 
%SVM_OBJECT_DECODING_FULL_EXPERIMENT_PROBLEMATIC Perform object decoding 
%(average of pairwise object decoding) using the SVM classifier. For the
%fixation task of subjects 1-4, where triggers were not properly recorded,
%and samples 22 and 56 had to be excluded. 
%
%Input: subject ID
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
[numTrials, trials_per_condition] = min_number_trials(timelock_triggers, numConditionsAll); %minimum number of trials per scene
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=100; 

%Deal with the excluded samples
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
   
    %only get the lower diagonal
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

%% Save the decoding accuracy
decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations

% add NaNs in lieu of condition 22
added_matrix_22_row = NaN(1,58,numTimepoints);
added_matrix_22_column = NaN(59,1,numTimepoints);
decodingAccuracy_with_22 = cat(1,decodingAccuracy_avg(1:21,:,:),added_matrix_22_row,...
    decodingAccuracy_avg(22:end,:,:));
decodingAccuracy_with_22 = cat(2,decodingAccuracy_with_22(:,1:21,:),added_matrix_22_column,...
    decodingAccuracy_with_22(:,22:end,:));

% add NaNs in lieu of condition 56
added_matrix_56_row = NaN(1,59,numTimepoints);
added_matrix_56_column = NaN(60,1,numTimepoints);
decodingAccuracy_with_56 = cat(1,decodingAccuracy_with_22(1:55,:,:),added_matrix_56_row,...
    decodingAccuracy_with_22(56:end,:,:));
decodingAccuracy_with_56 = cat(2,decodingAccuracy_with_56(:,1:55,:),added_matrix_56_column,...
    decodingAccuracy_with_56(:,56:end,:));

%rename 
decodingAccuracy_avg = decodingAccuracy_with_56;
save(sprintf('/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/%s/svm_decoding_accuracy_%s.mat',subname,task_name),'decodingAccuracy_avg');

