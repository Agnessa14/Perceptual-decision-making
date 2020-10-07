function dth_SVM_pilot_2(subn)
%DTH_SVM_PILOT_2 Performs the distance-to-hyperplane analysis using SVM on
%a balanced dataset. 
%
%Input: subject ID (integer)
%
%Output: 
%   - NxP matrix of decision values, where N is the number of conditions
%   and P is the number of timepoints.
%   - Nx1 vector of RTs. 
%   - Nx1 vector of minimum trial #s. 

%% Set-up prereqs
%add paths
addpath(genpath('/home/agnek95/CoSMoMVPA'));
addpath(genpath('/scratch/agnek95/PDM/DATA/DATA_PILOT_2'));
addpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS/OTHER');
addpath(genpath('/home/agnek95/SMST/PDM/ANALYSIS/'));
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
decodingAccuracy=NaN(numPermutations,numConditions,numConditions,numTimepoints);

%% Running the MVPA
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,timelock_data.trialinfo,numTrials,timelock_data.trial);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);
    
    disp('Permuting the conditions')
    data = data(randperm(size(data,1)));
    
    
    
    
end

   
    
