function SVM_object_decoding_PDM(subject,subset) 
%SVM_OBJECT_DECODING_PDM Perform object decoding (average of pairwise object decoding) using the SVM classifier. 
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
timelock_data = timelock;
timelock_data.trialinfo = timelock_data.trialinfo(behav.RT>0 & behav.points==1); %triggers
timelock_data.trial = timelock_data.trial(behav.RT>0 & behav.points==1,:,:); %actual data
timelock_data.sampleinfo = timelock_data.sampleinfo(behav.RT>0 & behav.points==1,:);

%% Define the required variables

numConditions = 60;
[numTrials, ~] = min_number_trials(timelock_data.trialinfo, numConditions); %minimum number of trials per scene
numTrials = floor(subset*numTrials);
numTimepoints = size(timelock_data.trial,3); %number of timepoints
numPermutations=100; 

%Preallocate 
decodingAccuracy=NaN(numPermutations,numConditions,numConditions,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    [data,~] = create_data_matrix_subset(numConditions,timelock_data.trialinfo,numTrials,timelock_data.trial);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data); 

    disp('Binning data into pseudotrials');
    numTrialsPerBin = 2;
    [pseudoTrials,numPTs] = create_pseudotrials(numTrialsPerBin,data);
   
    %only get the lower diagonal
    for condA=1:numConditions-1 %1:65
        for condB = condA+1:numConditions %2:66
            for timePoint = 1:numTimepoints 
                disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                    ', timepoint ->',num2str(timePoint)]);
                
                % L-1 pseudo trials go to training set, the Lth to testing set
                if numPTs == 2
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
save(sprintf('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/%s/subset_%s_svm_decoding_accuracy',subname,num2str(subset)),'decodingAccuracy_avg');

