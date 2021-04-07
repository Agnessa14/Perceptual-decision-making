function dth_pseudotrials_SVM_full_experiment_removed_scene(subject,task)
%DTH_PSEUDOTRIALS_SVM_FULL_EXPERIMENT_REMOVED_SCENE Performs the distance-to-hyperplane analysis using SVM on
%a balanced dataset. Instead of creating pseudoconditions out of scenes,
%the trials from across conditions are lumped into pseudotrials.  For the
%categorization task of subject 10 who only had 3 correct trials in the categorization task: removing all trials for scene 47. 
%
%
%Input: subject ID (integer), task (1=categorization, 2=fixation)
%
%Output: 
%   - NxP matrix of decision values, where N is the number of conditions
%   and P is the number of timepoints.
%   - Nx1 vector of RTs. 
%   - Nx1 vector of minimum trial #s. 

%% Add paths
%toolboxes and helper functions
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
addpath('/home/agnek95/OR/TOOLBOX/MVNN/MEG_SVM_decoding_MVNN'); %MVNN toolbox
addpath(genpath('/home/agnek95/OR/ANALYSIS/DECODING/libsvm')); %libsvm toolbox
addpath('/home/agnek95/OR/ANALYSIS/DECODING'); %MVNN function
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

subname = get_subject_name(subject);
task_name = get_task_name(task);

%check if there's a directory for that subject, otherwise create one
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
if ~isfolder(fullfile(results_dir,subname))
    mkdir(results_dir,subname);
end

%data and results
data_dir = fullfile('/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT/',subname);
results_dir = fullfile(results_dir,subname);
addpath(genpath(data_dir));
addpath(results_dir);

%% Prepare data
%load data
load(fullfile(data_dir,sprintf('timelock_%s',task_name))); %eeg
load(fullfile(data_dir,sprintf('preprocessed_behavioural_data_%s',task_name)));

%only keep the trials with a positive RT & correct response
timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); 
timelock_data = timelock.trial(behav.RT>0 & behav.points==1,:,:); 

%% Define the required variables
numConditionsAll = 60;
num_categories = 2; %categories to decode
numTimepoints = size(timelock_data,3); %number of timepoints
numPermutations=100; 
[~, trials_per_condition] = min_number_trials(timelock_triggers, numConditionsAll); %minimum number of trials per scene
removed_condition = find(trials_per_condition==min(trials_per_condition));
low_minnumtrials = min(trials_per_condition);
numTrials = min(trials_per_condition(trials_per_condition>low_minnumtrials));

%exclude trials from scene 47
included_conditions = find(trials_per_condition>=numTrials);
numConditionsIncluded = numel(included_conditions);

%Preallocate 
decisionValues = NaN(numPermutations,numConditionsIncluded,numTimepoints);
if removed_condition<=30
    num_conditions_artificial = 29;
    num_conditions_natural = 30; 
else
    num_conditions_artificial = 30;
    num_conditions_natural = 29; 
end

%% Running the MVPA
rng('shuffle');
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditionsAll,timelock_triggers,numTrials,timelock_data);
    data = data(included_conditions,:,:,:); 
    
    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);

    disp('Split into artificial and natural');
    data_artificial = data(1:num_conditions_artificial,:,:,:); 
    data_natural = data(num_conditions_artificial+1:end,:,:,:);
       
    disp('Training set: Reshape by taking the trials from all conditions for each category');
    size_data_artificial = size(data_artificial);
    size_data_natural = size(data_natural);
    data_artificial_reshaped = reshape(data_artificial,...
        [size_data_artificial(1)*size_data_artificial(2),size_data_artificial(3),size_data_artificial(4)]);
    data_natural_reshaped = reshape(data_natural,...
        [size_data_natural(1)*size_data_natural(2),size_data_natural(3),size_data_natural(4)]);
    
    disp('Permute the trials')
    data_artificial_permuted = data_artificial_reshaped(randperm(size(data_artificial_reshaped,1)),:,:);
    data_natural_permuted = data_natural_reshaped(randperm(size(data_natural_reshaped,1)),:,:);
    
    disp('Take the same number of trials for both categories')
    num_trials_natural = size(data_natural_permuted,1);
    data_artificial_permuted = data_artificial_permuted(1:num_trials_natural,:,:);
    
    disp('Put both categories into one matrix');
    data_both_categories = NaN([num_categories,size(data_artificial_permuted)]);
    data_both_categories(1,:,:,:) = data_artificial_permuted;
    data_both_categories(2,:,:,:) = data_natural_permuted;
    
    disp('Split into pseudotrials: Artificial scenes');
    numTrialsPerBin = 20; %try different combinations of bins/numTrialsPerBin
    [bins,numBins] = create_pseudotrials(numTrialsPerBin,data_both_categories);
        
    disp('Testing set: Average over trials');
    data_artificial_avg = squeeze(mean(data_artificial,2));
    data_natural_avg = squeeze(mean(data_natural,2));
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,:,:,t)); squeeze(bins(2,:,:,t))];  %train on all pseudotrials
        testing_data  = [squeeze(data_artificial_avg(:,:,t)); squeeze(data_natural_avg(:,:,t))]; %test on all conditions 
       
        labels_train  = [ones(numBins,1);2*ones(numBins,1)]; 
        labels_test   = [ones(num_conditions_artificial,1);2*ones(num_conditions_natural,1)]; 
        
        disp('Train the SVM');
        train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
        model=svmtrain_01(labels_train,training_data,train_param_str); 
        
        disp('Test the SVM');
        [~, ~, decision_values] = svmpredict(labels_test,testing_data,model);  
        
        disp('Putting the decision values into the big matrix');
        decisionValues(perm,:,t) = abs(decision_values);
    end
    toc
end


%% Add NaN to the removed scene
decisionValues_Avg = squeeze(mean(decisionValues,1));
DV_1 = [decisionValues_Avg(1:removed_condition-1,:);NaN(1,numTimepoints);decisionValues_Avg(removed_condition:end,:)];
decisionValues_Avg = DV_1;

%% Save the decision values
save(fullfile(results_dir,sprintf('dth_pseudotrials_svm_decisionValues_%s',task_name)),'decisionValues_Avg');

%% Get the average (over trials) reaction time for each condition
RT_per_condition = NaN(numConditionsAll,1);
RT_correct = behav.RT(behav.RT > 0 & behav.points == 1);

for c = 1:numConditionsAll
    RT_per_condition(c) = mean(RT_correct(timelock_triggers==c));
end

%Add NaN to the removed scene
RT_per_condition(removed_condition) = NaN;
filename_RT = sprintf('RTs_correct_trials_%s.mat',task_name);
save(fullfile(results_dir,filename_RT),'RT_per_condition');

end
   
    