function bootstrap_peak_latency_dth(task_distance,task_RT)
%BOOTSTRAP_PEAK_LATENCY Apply bootstrapping to calculate the 95% confidence
%interval of peak distance-to-hyperplane latency.
%
%Input: task of EEG distances (1=categorization,2=distraction), task of RT
%(same), scene conditions ('artificial','natural' or 'both')
%
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
% results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
                                    %%%%% SETUP %%%%%
%% Load all subjects' data
task_name_distance = get_task_name(task_distance);

for c = 1:3 %artificial, natural, both
    switch c
        case 1
            conditions = 'artificial';
        case 2 
            conditions = 'natural';
        case 3
            conditions = 'both';
    end
    filename = sprintf('all_subjects_dth_%s',task_name_distance,conditions);
    if task_distance~=task_RT
        filename = sprintf('all_subjects_dth_%s_distance_%s',task_name_distance,conditions);
    end
    load(fullfile(results_avg,filename));
    
end

                                  %%%%% BOOTSTRAPPING %%%%%

%% 0) Prepare required variables
num_bootstrap_samples = 1000;
num_datasets = size(decoding_accuracies_all,1);
peak_latency_all_samples = NaN(num_bootstrap_samples,1);

%% 1) Create the bootstrap samples
rng(0,'twister');
max_dataset = num_datasets; %for the random dataset index generator
min_dataset = 1;
for bs = 1:num_bootstrap_samples
    datasets = NaN(size(decoding_accuracies_all));
    for d = 1:num_datasets
        idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
        datasets(d,:) = decoding_accuracies_all(idx,:);
    end
    avg_datasets = squeeze(mean(datasets,1));
    peak_latency_all_samples(bs) = find(avg_datasets==max(avg_datasets),1);
end

%% 2) Get mean bootstrapped peak latency
peak_latency = round(mean(peak_latency_all_samples));

%% 3) Get 95% confidence interval
confidence_interval = NaN(2,1);
confidence_interval(1) = prctile(peak_latency_all_samples,2.5);
confidence_interval(2) = prctile(peak_latency_all_samples,97.5);




