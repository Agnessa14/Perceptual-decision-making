function [peak_latency, confidence_interval] = bootstrap_peak_latency(subjects,task,analysis)
%BOOTSTRAP_PEAK_LATENCY Apply bootstrapping to calculate the 95% confidence
%interval of peak decoding latency.
%
%Input: subject IDs, task (1=categorization,2=distraction), analysis
%('object_decoding' or 'category_decoding')
%
%Output: peak latency (in ms), confidence interval (onset and offset, in
%ms)
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
% results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Pre-allocate + setup for loading
task_name = get_task_name(task);
sorted_subjects = sort(subjects);
numTimepoints = 200;
decoding_accuracies_all_subjects = NaN(sorted_subjects(end),numTimepoints);
if strcmp(analysis,'object_decoding')
    filename = 'svm_decoding_accuracy';
    var_name = 'decodingAccuracy_avg';
end

%% Load the subject-specific decoding curves 
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,sprintf('%s_%s.mat',filename,task_name)));
    decoding_accuracies_all_subjects(subject,:) = squeeze(nanmean(nanmean(eval(var_name),1),2)); %decodingAccuracy_avg; average over conditions
end 







