function pilot_2_dth_lda(subject)
%PILOT_2_DTH_LDA Obtain distances to the hyperplane from time-series data, 
%using a balanced dataset (same number of trials for each sample).
%
%Input: subject ID
%
%Output: NxP decision values (same rank as distances, which is all we
%need). N is the number of conditions and P is the number of timepoints. 
%
%
%% Set-up prereqs
%add paths
addpath(genpath('/home/agnek95/CoSMoMVPA'));
addpath(genpath('/scratch/agnek95/PDM/DATA/DATA_PILOT_2'));
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS/'));
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

% reset citation list - is this necessary
cosmo_check_external('-tic');

%subject name string
subname = get_subject_name(subject);
 
%% Prepare data
%load eeg and behavioural data
data_dir = sprintf('/scratch/agnek95/PDM/DATA/DATA_PILOT_2/%s/',subname);
load(fullfile(data_dir,'timelock')); %eeg
load(fullfile(data_dir,'preprocessed_behavioural_data'));

%only keep the trials with a positive RT & correct response
tl = timelock;
tl.trialinfo = tl.trialinfo(behav.RT>0 & behav.points==1);
tl.trial = tl.trial(behav.RT>0 & behav.points==1,:,:);
tl.sampleinfo = tl.sampleinfo(behav.RT>0 & behav.points==1,:);

%normalize EEG data: each channel of each trial
normalized_data = NaN(size(tl.trial));
for trial = 1:size(tl.trialinfo,1)
    for channel = 1:size(tl.trialinfo,2)
        normalized_data(trial,channel,:) = normalize(squeeze(tl.trial(trial,channel,:)));
    end
end       
    
% convert to cosmomvpa struct
ds_tl = cosmo_meeg_dataset(tl); %returns struct with field samples (numtrials x numelectrodes*numtimepoints)
samples_all = ds_tl.samples;
trialinfo_all = ds_tl.sa.trialinfo;

%% Set up some variables
num_permutations = 100;
num_conditions = 60;
num_timepoints = size(tl.trial,3);
decision_values = NaN(num_permutations,num_conditions,num_timepoints);
last_category1_sample = 30;

%% Create a balanced set
%collect the trials and number of trials for each condition
trials_condition = NaN(num_conditions,60); %max possible
num_trials_condition = NaN(num_conditions,1);

for c = 1:num_conditions
    trials_c = find(trialinfo_all==c);
    trials_condition(c,1:numel(trials_c)) = trials_c;
    num_trials_condition(c) = numel(trials_c);
end

%take the minimum number of trials
min_num_trials = min(num_trials_condition(num_trials_condition>0));

%% Loop   
for p = 1:num_permutations 
    
    %randomly permute the trials
    randomized_trials = NaN(size(trials_condition));
    for c = 1:num_conditions
        permuted = randperm(num_trials_condition(c));
        randomized_trials(c,1:numel(permuted)) = trials_condition(c,permuted);
    end
    
    %take the minimum number of trials of each of the 60 samples
    selected_trials = reshape(randomized_trials(:,1:min_num_trials)',[],1);
    selected_trials(isnan(selected_trials)) = [];
    
    %put back into the ds structure
    ds_tl.sa.trialinfo = trialinfo_all(selected_trials);
    ds_tl.samples = samples_all(selected_trials,:);
    
    %separate into artificial and natural trials
    art_trials = 1:numel(selected_trials)/2; %first half of the selected trials
    nat_trials = numel(selected_trials)/2 + 1:numel(selected_trials); %second half of the selected trials
    
    % set the targets (labels for the classifier): 1 = artificial, 2 = natural
    ds_tl.sa.targets = NaN(size(ds_tl.sa.trialinfo));
    for t = 1:size(ds_tl.sa.targets)
        if ds_tl.sa.trialinfo(t) <= last_category1_sample
            ds_tl.sa.targets(t) = 1; %artificial
        else
            ds_tl.sa.targets(t) = 2; %natural
        end
    end
    
    %split into training and testing: create chunks - 1 = training and 2 = testing
    training = [art_trials(1:round(numel(art_trials)/2));nat_trials(1:round(numel(nat_trials)/2))];
    testing = [art_trials(round(numel(art_trials)/2) + 1:end); nat_trials(round(numel(nat_trials)/2) + 1:end)];
    
    ds_tl.sa.chunks = NaN(numel(ds_tl.sa.targets),1);
    ds_tl.sa.chunks(training) = 1;
    ds_tl.sa.chunks(testing) = 2;
    
    % just to check everything is ok
    cosmo_check_dataset(ds_tl);

    %define searchlight parameters            
    time_radius=0; %on each timepoint
    time_nbrhood=cosmo_interval_neighborhood(ds_tl,'time',...
                                                'radius',time_radius);
    % print some info
    nbrhood_nfeatures=cellfun(@numel,time_nbrhood.neighbors);
    fprintf('Features have on average %.1f +/- %.1f neighbors\n', ...
                mean(nbrhood_nfeatures), std(nbrhood_nfeatures));

    % only keep features with at least 10 neighbors
    center_ids=find(nbrhood_nfeatures>10);

    % Use the cosmo_cross_validation_measure and set its parameters
    % (classifier and partitions) in a measure_args struct.
    measure = @modified_cosmo_decisionvalues_measure;
    measure_args = struct();

    % Define which classifier to use, using a function handle.
    measure_args.classifier = @modified_cosmo_classify_lda;
    partitions=cosmo_nchoosek_partitioner(ds_tl,'half');
    measure_args.partitions = partitions;
  
    %print stuff
    fprintf('There are %d partitions\n', numel(measure_args.partitions.train_indices));
    fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            measure_args.partitions.train_indices)));
    fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            measure_args.partitions.test_indices)));

    % run searchlight
    sl_tl_ds=cosmo_searchlight(ds_tl,time_nbrhood,measure,measure_args,'center_ids',center_ids);
    

    %avg over trials for each condition
    absolute_value_dvs = abs(sl_tl_ds.samples); 
    t = 1:num_timepoints;
    for c = 1:num_conditions 
        trials_c = ds_tl.sa.trialinfo == c;
        decision_values(p,c,:) = arrayfun(@(x) squeeze(mean(absolute_value_dvs(trials_c,x),1)),t);
    end  
end

%% Save
%decision values
averaged_decision_values = squeeze(mean(decision_values,1)); 
save_path = fullfile('/home/agnek95/SMST/PDM_PILOT_2/RESULTS/',subname);
save(fullfile(save_path,'dth_lda'),'averaged_decision_values');    

%% RTs per condition
RT_per_condition = NaN(num_conditions,1);
RT_correct = behav.RT(behav.RT > 0 & behav.points == 1);

for c = 1:num_conditions
    RT_per_condition(c) = mean(RT_correct(trialinfo_all==c));
end

save(fullfile(save_path,'RTs_correct_answers'),'RT_per_condition');    

%% Number of trials per condition
save(fullfile(save_path,'num_trials_per_condition'),'num_trials_condition'); 

