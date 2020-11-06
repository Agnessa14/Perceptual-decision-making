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
numPermutations=100; 

%Preallocate 
decisionValues_Artificial=NaN(numPermutations,num_conditions_per_category,numTimepoints);
decisionValues_Natural = NaN(numPermutations,num_conditions_per_category,numTimepoints);

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
       
    disp('Average over trials');
    data_natural_avg = squeeze(mean(data_natural,2));
    data_artificial_avg = squeeze(mean(data_artificial,2));
  
    disp('Permute the conditions (scenes)');
    conditions_order = randperm(num_conditions_per_category)';
    data_natural_avg = data_natural_avg(conditions_order,:,:);
    data_artificial_avg = data_artificial_avg(conditions_order,:,:);
    
    disp('Put both categories into one matrix');
    data_both_categories = NaN([num_categories,size(data_artificial_avg)]);
    data_both_categories(1,:,:,:) = data_natural_avg;
    data_both_categories(2,:,:,:) = data_artificial_avg;
    
    disp('Split into bins of scenes');
    numScenesPerBin = 4;
    [bins,numBins] = create_pseudotrials(numScenesPerBin,data_both_categories);
    
    for t = 1:numTimepoints
        disp('Split into training and testing');
        training_data = [squeeze(bins(1,:,:,t)); squeeze(bins(2,:,:,t))];  %train on all pseudoconditions (bins) 
        testing_data  = [squeeze(data_both_categories(1,:,:,t)); squeeze(data_both_categories(2,:,:,t))]; %test on all conditions 
       
        labels_train  = [ones(numBins,1);2*ones(numBins,1)]; 
        labels_test   = [ones(num_conditions_per_category,1);2*ones(num_conditions_per_category,1)]; 
        
        disp('Train the SVM');
        train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
        model=svmtrain_01(labels_train,training_data,train_param_str); 
        
        disp('Test the SVM');
        [~, ~, decision_values] = svmpredict(labels_test,testing_data,model);  
        
        disp('Putting the decision values into the big matrix');
        for c = 1:num_conditions_per_category
            condition = conditions_order(c);
            decisionValues_Natural(perm,condition,t) = abs(decision_values(c));
            decisionValues_Artificial(perm,condition,t) = abs(decision_values(c+num_conditions_per_category));
        end
    end
    toc
end

%% Save the decision values
decisionValues_Natural_Avg = squeeze(mean(decisionValues_Natural,1)); %avg over permutations
decisionValues_Artificial_Avg = squeeze(mean(decisionValues_Artificial,1)); %avg over permutations
decisionValues_Avg = [decisionValues_Natural_Avg;decisionValues_Artificial_Avg];
results_dir = fullfile('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/',subname);
save(fullfile(results_dir,'witherrortrials_svm_ritchie_decisionValues'),'decisionValues_Avg');

end
   
    
