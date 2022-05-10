function SVM_object_decoding_full_experiment_searchlight(subject,task)
%SVM_OBJECT_DECODING_FULL_EXPERIMENT Perform object decoding (average of pairwise object decoding) using the SVM classifier.
%
%Input: subject ID, task (1=categorization, 2=fixation)
%
%Output: NxNxP vector of accuracies in %, where N is the number of conditions and
%P is the number of timepoints.
%
%Author: Agnessa Karapetian (2021), modified by Muthukumar Pandaram (2022)

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
data_dir = fullfile('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/',subname);
load(fullfile(data_dir,sprintf('timelock_%s',task_name))); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)));

%only keep the trials with a positive RT & correct response
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); %actual data

%% Define the required variables
times = -195:100:705;
time_2_idx = (times/5)+40;

numConditions = 60;
[numTrials, ~] = min_number_trials(timelock_triggers, numConditions); %minimum number of trials per scene
numPermutations=100;
numChannels = 63; %regardless of missing channels
chanIdx = 1:numChannels; % select all channels

%Preallocate
% decodingAccuracy=NaN(numPermutations,numConditions,numConditions,numel(times),numChannels);

%Preallocate: for memory purposes, two matrices for 50 permutations each
decodingAccuracy_1 = NaN(round(numPermutations/2),numConditions,numConditions,numel(times),numChannels);
decodingAccuracy_2 = NaN(round(numPermutations/2),numConditions,numConditions,numel(times),numChannels);

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
    
    % Modifing neighbourhood and pseudotrial matrix due to missing channels
    [neighbourhoods, missing_channel_ids, pseudoTrials] = correct_neighbourhood_map(timelock.label,pseudoTrials);
    
    %only get the upper diagonal
    for condA=1:numConditions-1 %1:59
        for condB = condA+1:numConditions %2:60
            for timePoint = 1:numel(time_2_idx)
                t = time_2_idx(timePoint);
                disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                    ', timepoint ->',num2str(timePoint)]);
                for iChan = chanIdx
                    if ~ismember(iChan,missing_channel_ids) || isempty(missing_channel_ids)                        
                        neighbours = neighbourhoods(iChan , :);
                        neighbours = neighbours(~isnan(neighbours));
                        
                        % L-1 pseudo trials go to training set, the Lth to testing set
                        training_data=[squeeze(pseudoTrials(condA,1:end-1,neighbours,t)) ; squeeze(pseudoTrials(condB,1:end-1,neighbours,t))]; %(numbins-1)xnum_channelsx1 each
                        testing_data=[squeeze(pseudoTrials(condA,end,neighbours,t))' ; squeeze(pseudoTrials(condB,end,neighbours,t))']; %1x63x1 each
                        
                        % class labels
                        labels_train=[ones(1,numPTs-1) 2*ones(1,numPTs-1)]; %one label for each pseudotrial
                        labels_test= [1 2];
                        
                        train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; %look up the parameters online if needed
                        model=svmtrain_01(labels_train',training_data,train_param_str);
                        [~, accuracy, ~,] = svmpredict(labels_test',testing_data,model);
%                         decodingAccuracy(perm,condA,condB,timePoint,iChan)=accuracy(1);
                        if perm <= 50
                            decodingAccuracy_1(perm,condA,condB,timePoint,iChan)=accuracy(1);
                        else
                            decodingAccuracy_2(perm-50,condA,condB,timePoint,iChan)=accuracy(1);
                        end

                    end
                end
            end
        end
    end
    toc
end

%% Save the decoding accuracies
% decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations

decodingAccuracy_avg_both = NaN(2,numConditions,numConditions,numel(times),numChannels);
decodingAccuracy_avg_both(1,:,:,:,:) = squeeze(mean(decodingAccuracy_1,1));
decodingAccuracy_avg_both(2,:,:,:,:) = squeeze(mean(decodingAccuracy_2,1));
decodingAccuracy_avg = squeeze(mean(decodingAccuracy_avg_both,1));
save(fullfile(results_dir,subname,sprintf('svm_searchlight_object_decoding_accuracy_%s.mat',task_name)),'decodingAccuracy_avg');    

end


