function SVM_object_decoding_PDM(subject) 
%SVM_OBJECT_DECODING_PDM Perform object decoding (average of pairwise object decoding) using the SVM classifier. 
%
%Input: subject ID
%
%Output: NxNxP vector of accuracies in %, where N is the number of conditions and
%P is the number of timepoints. 
%
%% Set-up prereqs
%add paths
addpath(genpath('/home/agnek95/CoSMoMVPA'));
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
[numTrials, ~] = minNumberTrials(timelock_data.trialinfo, numConditions); %minimum number of trials per scene
numElectrodes = size(timelock_data.trial,2); %number of electrodes
numTimepoints = size(timelock_data.trial,3); %number of timepoints
numPermutations=100; 

%Preallocate 
data = NaN(numConditions,numTrials,numElectrodes,numTimepoints);
decodingAccuracy=NaN(numPermutations,numConditions,numConditions,numTimepoints);
    
%% Decoding
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = createDataMatrix(numConditions,timelock_data.trialinfo,numTrials,data,timelock_data.trial);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data,numTimepoints,numConditions,numElectrodes,numTrials); 

    disp('Binning data into pseudotrials');
    numTrialsPerBin = 10;
    [pseudoTrials,numPTs] = createPseudotrials(numTrialsPerBin,numTrials,numConditions,numElectrodes,numTimepoints,data);
   
    %only get the lower diagonal
    for condA=1:numConditions-1 %1:65
        for condB = condA+1:numConditions %2:66
            for timePoint = 1:numTimepoints 
                disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                    ', timepoint ->',num2str(timePoint)]);
                % to train a SVM, you need to create a training & testing set
                % L-1 pseudo trials go to testing set, the Lth to training set
                training_data=[squeeze(pseudoTrials(condA,1:end-1,:,timePoint)) ; squeeze(pseudoTrials(condB,1:end-1,:,timePoint))]; %(numbins-1)x63x200 each
                testing_data=[squeeze(pseudoTrials(condA,end,:,timePoint))' ; squeeze(pseudoTrials(condB,end,:,timePoint))']; %1x63x200 each            
                
                % class labels
                labels_train=[ones(1,numPTs-1) 2*ones(1,numPTs-1)]; %they are L times 1, and L times 2 in a row (you just label any classes as 1 and 2, as the SVM at any time point only sees 2 classes)
                labels_test= [1 2]; % you will test on one pseudotrial from condA, and one from cond B
                
                %parameters for the SVM: s0: C-SVC; t0: linear kernel; 
                %b0: not train for probability estimates; 
                %c1: (cost) C parameter of C-SVC is 1; -q no output
                train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; 
                model=svmtrain_01(labels_train',training_data,train_param_str); 
                [~, accuracy, ~,] = svmpredict(labels_test',testing_data,model);
                decodingAccuracy(perm,condA,condB,timePoint)=accuracy(1);                

            end % time point
        end % cond A
    end %cond B
    toc
end % permutation

%% Save the decoding accuracy
decodingAccuracy_avg = squeeze(mean(decodingAccuracy,1)); %average over permutations
save(sprintf('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/%s/svm_decoding_accuracy',subname),'decodingAccuracy_avg');

%% Helper functions

%Determine number of trials persa object to find the minimum number of trials
function [N, numTrials] = minNumberTrials(conditionsList,numConditions)
numTrials = [];
for cond = 1:numConditions
    numTrials(cond) = numel(find(conditionsList==cond));
end
N = min(numTrials);

%Create the matrix of data
function [dataMatrix] = createDataMatrix(numConditions,triggers,minNumTrials,dataMatrix,data) 
%Create the matrix containing voltage values for each object, electrode and
%timepoint and the matrix containing behavioural data for each condition and trial
for c = 1:numConditions
    trialsMat = find(triggers==c);
    trialsMat = trialsMat(randperm(numel(trialsMat))); %randomize
    trialsMat = trialsMat(1:minNumTrials);
    dataMatrix(c,:,:,:) = data(trialsMat,:,:); 
end

%Create the matrix with pseudotrials
function [dataMatrixPT,numPseudotrials] = createPseudotrials(numTrialsBin,minNumTrials,numConditions,numElectrodes,numTimepoints,dataMatrix) 
numPseudotrials=round(minNumTrials/numTrialsBin); %must be an integer; number of pseudo trials
dataMatrixPT=NaN(numConditions,numPseudotrials,numElectrodes,numTimepoints); %pre-allocate memory
for step=1:numPseudotrials-1 %average by steps
    trial_selector=(1+(step-1)*numTrialsBin):(numTrialsBin+(step-1)*numTrialsBin); %select trials to be averaged
    dataMatrixPT(:,step,:,:)= nanmean(dataMatrix(:,trial_selector,:,:),2); %assign pseudo trial to dataMatrixPT
end
dataMatrixPT(:,numPseudotrials,:,:) = nanmean(dataMatrix(:,(1+(numPseudotrials-1)*numTrialsBin):end,:,:),2);


