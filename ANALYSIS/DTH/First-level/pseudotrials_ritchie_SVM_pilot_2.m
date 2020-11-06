function pseudotrials_ritchie_SVM_pilot_2(subject)
%PSEUDOTRIALS_RITCHIE_DTH_SVM_PILOT_2 Performs the distance-to-hyperplane analysis using SVM on
%a balanced dataset. Instead of creating pseudoconditions out of scenes,
%the trials from across conditions are lumped into pseudotrials.
%
%Input: subject ID (integer)
%
%Output: 
%   - NxP matrix of decision values, where N is the number of conditions
%   and P is the number of timepoints.
%   - Nx1 vector of RTs. 
%   - Nx1 vector of minimum trial #s. 

%% Add paths
%toolboxes and helper functions
addpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS/OTHER');
addpath(genpath('/home/agnek95/SMST/PDM/ANALYSIS/'));
addpath('/home/agnek95/OR/TOOLBOX/MVNN/MEG_SVM_decoding_MVNN'); %MVNN toolbox
addpath(genpath('/home/agnek95/OR/ANALYSIS/DECODING/libsvm')); %libsvm toolbox
addpath('/home/agnek95/OR/ANALYSIS/DECODING'); %MVNN function
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

subname = get_subject_name(subject);

%% Prepare data
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/ritchie_subject_%s',subname);
addpath(genpath(data_dir));
load(fullfile(data_dir,sprintf('s%s_PCA_S1_50hz.mat',subname))); %eeg
 
%take only the needed data 
triggers = data.TrialList(:,1);
task = data.TrialList(:,4);
trials_final = [];
for t = 1:size(data.TrialList,1)
    if task(t) == 1 %only take the active task and the correct trials
        trials_final = [trials_final;t];
    end
end
triggers_final = triggers(trials_final);
data_final = permute(data.class_dat(trials_final,:,:),[1 3 2]);

%% Define the required variables
numConditions = 24;
num_categories = 2; %categories to decode
num_conditions_per_category = numConditions/num_categories;

[numTrials, ~] = min_number_trials(triggers_final, numConditions); 
numTimepoints = size(data_final,3);
numPermutations=1; 

%Preallocate 
decisionValues = NaN(numPermutations,numConditions,numTimepoints);

%% Running the MVPA
for perm = 1:numPermutations
    tic   
    disp('Creating the data matrix');
    data = create_data_matrix(numConditions,triggers_final,numTrials,data_final);

    disp('Performing MVNN');
    data = multivariate_noise_normalization(data);

    disp('Split into artificial and natural');
    data_natural = data(1:num_conditions_per_category,:,:,:); %in Ritchie's dataset, the natural conditions come first
    data_artificial = data(num_conditions_per_category+1:end,:,:,:);
       
    disp('Training set: Reshape by taking the trials from all conditions for each category');
    size_data_category = size(data_natural);
    data_natural_reshaped = reshape(data_natural,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    data_artificial_reshaped = reshape(data_artificial,...
        [size_data_category(1)*size_data_category(2),size_data_category(3),size_data_category(4)]);
    
    disp('Permute the trials')
    data_natural_permuted = data_natural_reshaped(randperm(size(data_natural_reshaped,1)),:,:);
    data_artificial_permuted = data_artificial_reshaped(randperm(size(data_artificial_reshaped,1)),:,:);
       
    disp('Put both categories into one matrix');
    data_both_categories = NaN([num_categories,size(data_natural_permuted)]);
    data_both_categories(1,:,:,:) = data_natural_permuted;
    data_both_categories(2,:,:,:) = data_artificial_permuted;
    
    disp('Split into pseudotrials');
    numTrialsPerBin = 20; %try different combinations of bins/numTrialsPerBin
    [bins,numBins] = create_pseudotrials(numTrialsPerBin,data_both_categories);
    
    disp('Testing set: Average over trials');
    data_natural_avg = squeeze(mean(data_natural,2));
    data_artificial_avg = squeeze(mean(data_artificial,2));
    data_testing_both = NaN([num_categories,size(data_natural_avg)]);
    data_testing_both(1,:,:,:) = data_natural_avg;
    data_testing_both(2,:,:,:) = data_artificial_avg;
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,:,:,t)); squeeze(bins(2,:,:,t))];  %train on all pseudotrials
        testing_data  = [squeeze(data_testing_both(1,:,:,t)); squeeze(data_testing_both(2,:,:,t))]; %test on all conditions 
       
        labels_train  = [ones(numBins,1);2*ones(numBins,1)]; 
        labels_test   = [ones(num_conditions_per_category,1);2*ones(num_conditions_per_category,1)]; 
        
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

%% Save the decision values
decisionValues_Avg = squeeze(mean(decisionValues,1));
results_dir = fullfile('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/',subname);
save(fullfile(results_dir,'true_pseudotrials_svm_ritchie_decisionValues'),'decisionValues_Avg');

end
   
    
