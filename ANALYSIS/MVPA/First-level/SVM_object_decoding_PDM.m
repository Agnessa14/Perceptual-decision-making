function SVM_object_decoding_PDM(subn,taskType) %task 1 = categorization, task 2 = fixation cross

%% Change directory to the current subject and add paths
if subn < 10
    subn = ['0',num2str(subn)];
else
    subn = num2str(subn);
end
cd(fullfile('/scratch/agnek95/PDM/DATA/DATA/',subn)); 

%Add the MVNN toolbox to the path
addpath('/home/agnek95/OR/TOOLBOX/MVNN/MEG_SVM_decoding_MVNN'); 

%Add the libsvm toolbox to the path attached
addpath(genpath('/home/agnek95/OR/ANALYSIS/DECODING/libsvm'));

%Add the folder with the MVNN function to the path
addpath('/home/agnek95/OR/ANALYSIS/DECODING');

%Add the extractfield.m folder to the path
addpath('/home/agnek95/SMST/PDM/ANALYSIS/PREPROCESSING');

%Add all of the subject's subfolders to the path
addpath(genpath(pwd));

%% Compare the triggers from the behavioural data structure to the triggers in the preprocessed EEG file. We need this to determine when each block (task) starts and ends

%Load the preprocessed data structure 
load('timelock.mat');

%load behavioural data
mainDirectoryBeh = dir(fullfile('/scratch/agnek95/PDM/DATA/DATA/',subn,'/Behavioural/','*.mat')); %create directory for the beh file

%Behavioural data is not properly ordered when the directory is created -
%need to fix this

filenames = extractfield(mainDirectoryBeh,'name');
for m = 1:numel(mainDirectoryBeh)
    parsedName = split(filenames{m},'_');
    mainDirectoryBeh(m).blockNum = str2double(parsedName{4});     
end

%order by block number
T = struct2table(mainDirectoryBeh); % convert the struct array to a table
sortedT = sortrows(T, 'blockNum'); % sort the table by 'DOB
mainDirectoryBeh = table2struct(sortedT);     
filenamesSorted = extractfield(mainDirectoryBeh,'name');

%create the list of triggers, reaction times and points from the behavioural data
triggers = [];
RT = [];
points = [];
for f = 1:numel(mainDirectoryBeh)
    load(filenamesSorted{f});
    triggers = [triggers; blockData.samples];
    RT = [RT; blockData.rt];
    points = [points; blockData.points];
end

%create a matrix indicating which task each trial belongs to
task = NaN(numel(triggers),1);
for b = 1:numel(mainDirectoryBeh)
    parsedName = split(filenamesSorted{b},'_');
    load(filenamesSorted{b});
    trialsBlock = numel(blockData.samples);
    if contains(parsedName{6},'Categorization') 
        task((b-1)*trialsBlock + 1:b*trialsBlock) = 1;
    elseif contains(parsedName{6}, 'Fixation') 
        task((b-1)*trialsBlock + 1:b*trialsBlock) = 2;
    end
end

%compare preprocessed eeg triggers with the behavioural triggers, remove the extra behavioural triggers
%remove all paperclip trials from the behavioural triggers and from the task array
paperclip = find(triggers==999);
triggers(paperclip) = [];
task(paperclip) = [];
RT(paperclip) = [];
points(paperclip) = [];
prepTriggers = timelock.trialinfo;
for i = 1:numel(prepTriggers)
    while prepTriggers(i) ~= triggers(i)
        disp(i);
        triggers(i) = [];
        task(i) = [];
        RT(i) = [];
        points(i) = [];
    end
end

%Double-check that eeg and behavioural triggers are equal
if isequal(round(prepTriggers),round(triggers)) %instead of equality or ismembertol
    disp('all good!');
else
    disp('oh no');
    keyboard;
end

%% Partition the data into the two tasks (fixation cross and categorization), including error trials 
Y = timelock.trial; %num_trials x num_electrodes x num_datapoints
catY = Y(task == 1,:,:); %categorizationData
fixY = Y(task == 2,:,:);%fixationCrossData

%Trigger information
catTriggers = triggers(task == 1);
fixTriggers = triggers(task == 2);

%Reaction times
catRT = RT(task == 1); 
fixRT = RT(task == 2);

%% Define the required variables
%Non-responses - remove from the RT and the eeg data matrices
if taskType == 1
    noResp = isnan(catRT);
    catY(noResp,:,:) = [];
    catTriggers(noResp) = [];
elseif taskType == 2
    noResp = isnan(fixRT);
    fixY(noResp,:,:) = [];
    fixTriggers(noResp) = [];
end

numConditions = 66;
M = numConditions; %number of scenes
if taskType == 1
    [N, numb] = minNumberTrials(catTriggers, M); %minimum number of trials per scene
    O = size(catY,2); %number of electrodes
    T = size(catY,3); %number of timepoints
elseif taskType == 2
    [N, numb] = minNumberTrials(fixTriggers, M); %minimum number of trials per scene
    O = size(fixY,2); %number of electrodes
    T = size(fixY,3); %number of timepoints
end
minNumArray = numb;

%Preallocate the data and RT matrices
C = NaN(M,N,O,T);

%Number of permutations
numPermutations=100; 

%Preallocate the decoding accuracy matrix
DA=NaN(numPermutations,M,M,T);
    
%% Decoding
for perm = 1:numPermutations
    tic
    %Create data matrix: num conditions x num trials x num electrodes x num timepoints
    %and RT matrix: num conditions x num trials
    disp('Creating the data matrix');
    if taskType == 1
        [C] = createDataMatrix(M,catTriggers,N,C,catY);
    elseif taskType == 2
        [C] = createDataMatrix(M,fixTriggers,N,C,fixY);
    end

    %Perform MVNN
    disp('Performing MVNN');
    C = multivariate_noise_normalization(C,T,M,O,N); 
    
    %Bin eeg data into pseudotrials
    disp('Binning data into pseudotrials');
    [pseudoTC,numPTs] = createPseudotrials(5,N,M,O,T,C);
   
    %only get the lower diagonal
    for condA=1:M-1 %1:65
        for condB = condA+1:M %2:66
            for timePoint = 1:T %all time points are independent
                disp(['Running the classification: 1st sample ->', num2str(condA), ', 2nd sample ->',num2str(condB),...
                    ', timepoint ->',num2str(timePoint)]);
                % to train a SVM, you need to create a training & testing set
                % L-1 pseudo trials go to testing set, the Lth to training set
                training_data=[squeeze(pseudoTC(condA,1:end-1,:,timePoint)) ; squeeze(pseudoTC(condB,1:end-1,:,timePoint))]; %3x63x200 each
                testing_data=[squeeze(pseudoTC(condA,end,:,timePoint))' ; squeeze(pseudoTC(condB,end,:,timePoint))']; %1x63x200 each            
                
                % you have to define class labels for the trials you will train and test
                labels_train=[ones(1,numPTs-1) 2*ones(1,numPTs-1)]; %they are L times 1, and L times 2 in a row (you just label any classes as 1 and 2, as the SVM at any time point only sees 2 classes)
                labels_test= [1 2]; % you will test on one pseudotrial from condA, and one from cond B
                
                %DTH
                %parameters for the SVM: s0: C-SVC; t0: linear kernel; 
                %b0: not train for probability estimates; 
                %c1: (cost) C parameter of C-SVC is 1; -q no output
                train_param_str= '-s 0 -t 0 -b 0 -c 1 -q'; 
                model=svmtrain_01(labels_train',training_data,train_param_str); 
                [~, accuracy, ~,] = svmpredict(labels_test',testing_data,model);
                               
                % output is predicted_labels, accuracy (in%) and decision values; you need accuracy only for now, you need to save this. One possibility would be in a huge matrix DA that has the dimensions of all loops it is in
                DA(perm,condA,condB,timePoint)=accuracy(1);

            end % time point
        end % cond A
    end %cond B
    toc
end % permutation

%Plot the decoding accuracy
DA_end = squeeze(mean(DA,1)); %average over permutations
if taskType == 1
    taskTStr = 'categorization';
elseif taskType == 2
    taskTStr = 'fixation';
end
save(['/home/agnek95/SMST/PDM/RESULTS/',subn, '/mvnn_svm_decoding_accuracy_pdm_', taskTStr],'DA_end');


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
    %Find all the trials for the condition c (object)
    trialsMat = find(triggers==c);

    %Randomize the trials
    trialsMat = trialsMat(randperm(numel(trialsMat)));

    %Take the minimum number of trials
    trialsMat = trialsMat(1:minNumTrials);

    %Create the matrix containing voltage values for each object, electrode and timepoint
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


