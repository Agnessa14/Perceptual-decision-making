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
numConditions = 60;
[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); %minimum number of trials per scene
numTimepoints = numResampledTps; 
numPermutations=50; 

%Preallocate 
% decodingAccuracy=NaN(numPermutations,numConditions,numConditions,numTimepoints,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_triggers,numTrials,downsampled_timelock_data);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); 

    disp('Binning data into pseudotrials');
    numTrialsPerBin = round(numTrials/6);
    [pseudoTrials,numPTs] = create_pseudotrials(numTrialsPerBin,data);
   
    %only get the upper diagonal
    for condA=1:numConditions-1 
        for condB = condA+1:numConditions%1:numConditions 
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

%% Save the decoding accuracy
timeg_decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
save(fullfile(results_dir,subname,sprintf('time_generalized_svm_object_decoding_%s.mat',task_name)),'timeg_decodingAccuracy_avg');

end