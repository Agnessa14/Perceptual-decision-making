function [peak_latency, confidence_interval] = bootstrap_peak_latency(decoding_accuracies_all)
%BOOTSTRAP_PEAK_LATENCY Apply bootstrapping to calculate the 95% confidence
%interval of peak decoding latency.
%
%Input: NxP matrix of decoding accuracies, where N is the number or
%subjects and P is the number of timepoints.
%
%Output: peak latency (in ms), confidence interval (onset and offset, in
%ms)
%

%% Paths
% addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
% results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';

                                    %%%%% SETUP %%%%%
%% Pre-allocate + setup for loading
% task_name = get_task_name(task);
% sorted_subjects = sort(subjects);
% numTimepoints = 200;
% decoding_accuracies_all = NaN(sorted_subjects(end),numTimepoints);
% if strcmp(analysis,'object_decoding')
%     filename = 'svm_decoding_accuracy';
% elseif strcmp(analysis,'category_decoding')
%     filename = 'svm_artificial_vs_natural_decoding_accuracy';
% end
% 
% %% Load the subject-specific decoding curves 
% for subject = subjects
%     subname = get_subject_name(subject);
%     subject_results_dir = fullfile(results_dir,subname);
%     load(fullfile(subject_results_dir,sprintf('%s_%s.mat',filename,task_name)));
%     if strcmp(analysis,'object_decoding')
%         decoding_accuracies_all(subject,:) = squeeze(nanmean(nanmean(decodingAccuracy_avg,1),2)); %decodingAccuracy_avg; average over conditions
%     elseif strcmp(analysis,'category_decoding')
%         decoding_accuracies_all(subject,:) = decodingAccuracy_avg;
%     end
% end 
% 
% %Remove nan subjects
% decoding_accuracies_all = decoding_accuracies_all(~isnan(decoding_accuracies_all(:,1)),:); %if timepoint 1 is empty, assume there is no data for that subject

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




