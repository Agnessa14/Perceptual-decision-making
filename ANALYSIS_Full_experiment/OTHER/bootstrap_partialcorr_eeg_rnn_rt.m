function bootstrap_partialcorr_eeg_rnn_rt(conditions)
%BOOTSTRAP_PARTIALCORR_EEG_RNN_RT Calculate the peak latencies and confidence intervals
%for the partial correlations of EEG-RT and RNN-RT.
%
%Input: conditions ('artificial','natural' or 'both')
%
%Output: peak latency and confidence interval
%
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
addpath(genpath('/home/agnek95/helper_functions')); %correlate from M. Hebart
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';

%% Load the EEG data (distances to the hyperplane for all subjects), rCNN RTs and human RTs 
load(fullfile(results_avg_dir,'distances_all_subjects.mat'),'distances_all_subjects');
load('/scratch/agnek95/PDM/DATA/RNN_RTs/RNN_RTs_entropy_threshold_0.02.mat','data');
rCNN_RTs = data;
load(fullfile(results_avg_dir,'RT_all_subjects_5_35_categorization.mat'),'RTs');
human_RTs = nanmedian(RTs,1);

%Select the appropriate scenes
if strcmp(conditions,'artificial')
    conds = 1:30;
elseif strcmp(conditions,'natural')
    conds = 31:60;
elseif strcmp(conditions,'both')
    conds = 1:60;
end
selected_EEG_distances = distances_all_subjects(:,conds,:);
selected_human_RTs = human_RTs(conds)';
selected_rCNN_RTs = rCNN_RTs(conds);

%Define some variables
num_subjects = size(distances_all_subjects,1);
num_timepoints = size(distances_all_subjects,3);

%% Calculate the partial correlations
size_data = [num_subjects,num_timepoints];
r_eeg_unique = NaN(size_data);
r_rcnn_unique = NaN(size_data);
r_eeg_full = NaN(size_data);
r_rcnn_full = NaN(size_data);
shared_rt_eeg = NaN(size_data);
shared_rt_rcnn = NaN(size_data);

for sub = 1:num_subjects

    RTs_median_sub = selected_human_RTs;
    selected_rCNN_RTs_sub = selected_rCNN_RTs;
    distances_eeg_sub = squeeze(selected_EEG_distances(sub,:,:)); 

    %remove excluded scene if needed
    if any(isnan(distances_eeg_sub))
        excluded_scene = find(isnan(distances_eeg_sub(:,1)));
        distances_eeg_sub(excluded_scene,:) = [];
        RTs_median_sub(excluded_scene) = [];
        selected_rCNN_RTs_sub(excluded_scene) = [];
    end
        
    for t = 1:num_timepoints
        distances_eeg_sub_t = squeeze(distances_eeg_sub(:,t)); 
        
        %unique = semipartial
        r_eeg_u = correlate([RTs_median_sub, distances_eeg_sub_t*-1, selected_rCNN_RTs_sub'],'type','spearman','method','partialcorr');
        r_eeg_unique(sub,t) = r_eeg_u(2,1);
        r_rcnn_u = correlate([RTs_median_sub, selected_rCNN_RTs_sub', distances_eeg_sub_t*-1],'type','spearman','method','partialcorr');
        r_rcnn_unique(sub,t) = r_rcnn_u(2,1);
        
        %full 
        r_eeg_full(sub,t) = corr(RTs_median_sub,distances_eeg_sub_t*-1,'type','spearman');
        r_rcnn_full(sub,t) = corr(RTs_median_sub,selected_rCNN_RTs_sub','type','spearman');
        
        %shared = full (normal) -semipartial
        shared_rt_eeg(sub,t) = r_eeg_full(sub,t).^2-r_eeg_unique(sub,t).^2;
        shared_rt_rcnn(sub,t) = r_rcnn_full(sub,t).^2-r_rcnn_unique(sub,t).^2;   %shared_rt_eeg=shared_rt_rcnn     
    end 
end


                                                %%%%% BOOTSTRAPPING %%%%%
%% 0) Prepare required variables
num_bootstrap_samples = 1000;
num_datasets = num_subjects;
peak_latency_all_samples = NaN(num_bootstrap_samples,1);     

%% 1) Create the bootstrap samples
rng(0,'twister');
max_dataset = num_datasets; %for the random dataset index generator
min_dataset = 1;

%do for 1) unique eeg 2) shared eeg and rcnn rt
bootstrap_stats.confidence_interval = NaN(2,2);
bootstrap_stats.peak_latency_true = NaN(2,1);
bootstrap_stats.peak_latency_bs = NaN(2,1);

for p = 1:2 
    if p == 1
        data = r_eeg_unique(:,41:end); %starting from 1ms - discard the baseline
    elseif p == 2
        data = shared_rt_eeg(:,41:end);
    end
    
    for bs = 1:num_bootstrap_samples
        datasets = NaN(size(data));
        for d = 1:num_datasets
            idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
            datasets(d,:) = data(idx,:);
        end
        avg_datasets = squeeze(mean(datasets,1));
        peak_latency_all_samples(bs) = find(avg_datasets==max(avg_datasets),1)*5;
    end

    %% 2.0) Get mean bootstrapped peak latency - just to know
    bootstrap_stats.peak_latency_bs(p) = round(mean(peak_latency_all_samples));

    %% 2.1) Get true peak latency
    bootstrap_stats.peak_latency_true(p) = find(mean(data,1)==max(mean(data,1)),1)*5; %peak latency of avg-across-subjects curve
 
    %% 3) Get 95% confidence interval
    bootstrap_stats.confidence_interval(p,1) = prctile(peak_latency_all_samples,2.5);
    bootstrap_stats.confidence_interval(p,2) = prctile(peak_latency_all_samples,97.5);
 
end    
%Save results
filename_save = sprintf('bootstrap_stats_partial_corr_eeg_rcnn_rt_%s_scenes',conditions);
save(fullfile(results_avg_dir,filename_save),'bootstrap_stats');

end


