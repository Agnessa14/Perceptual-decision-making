function bootstrap_decoding_peak_difference(subjects,analysis)
%BOOTSTRAP_DECODING_PEAK_DIFFERENCE Calculate the confidence intervals of the difference between peak
%latencies of categorization and distraction decoding.
%
%Input: subject IDs, analysis ('object_decoding' or 'category_decoding')
%
%Output: confidence interval for the difference
%
%Author: Agnessa Karapetian, 2022
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate
numConditions = 60;
numTimepoints = 200;
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects_cat = NaN(max(subjects),numConditions,numConditions,numTimepoints);
    decoding_accuracies_all_subjects_dis = NaN(max(subjects),numConditions,numConditions,numTimepoints);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects_cat = NaN(max(subjects),numTimepoints);
    decoding_accuracies_all_subjects_dis = NaN(max(subjects),numTimepoints);
end

%% Loop: collect results from all subjects
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        cat_filename = 'svm_decoding_accuracy_categorization.mat';
        dis_filename = 'svm_decoding_accuracy_fixation.mat';
        all_dimensions = repmat({':'},1,3); %conditions x conditions x timepoints
    elseif strcmp(analysis,'category_decoding')
        cat_filename = 'cross_validated_dth_pseudotrials_svm_decodingAccuracy_categorization.mat';
        dis_filename = 'cross_validated_dth_pseudotrials_svm_decodingAccuracy_fixation.mat';
        all_dimensions = {':'}; % timepoints
    end
    load(fullfile(subject_results_dir,cat_filename),'decodingAccuracy_avg');
    decoding_accuracies_all_subjects_cat(subject,all_dimensions{:}) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,dis_filename),'decodingAccuracy_avg');
    decoding_accuracies_all_subjects_dis(subject,all_dimensions{:}) = decodingAccuracy_avg;
end

%% For stats matrix: num subjects x num timepoints (average over conditions)
if strcmp(analysis,'object_decoding')
    for_stats_cat = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects_cat,2),3));
    for_stats_dis = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects_dis,2),3));
elseif strcmp(analysis,'category_decoding')
    for_stats_cat = decoding_accuracies_all_subjects_cat;
    for_stats_dis = decoding_accuracies_all_subjects_dis;
end

%% Difference between tasks
for_stats_data_cat = for_stats_cat(subjects,:)-50;
for_stats_data_dis = for_stats_dis(subjects,:)-50;

%%% Boostrap %%%
                              
%% 0) Prepare required variables
bootstrap_peak_difference.num_bootstrap_samples = 1000;
num_datasets = numel(subjects);
peak_latency_all_samples = NaN(bootstrap_peak_difference.num_bootstrap_samples,2);
max_dataset = num_datasets; %for the random dataset index generator
min_dataset = 1;

%% 1) Create the bootstrap samples & calculate the peak difference
rng(0,'twister');
%for each bootstrap sample, create a new dataset (N=30) and
%calculate the peak latency
for bs = 1:bootstrap_peak_difference.num_bootstrap_samples
    
    datasets_cat = NaN(size(for_stats_data_cat));
    datasets_dis = NaN(size(for_stats_data_dis));
    
    for d = 1:num_datasets
        idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
        datasets_cat(d,:) = for_stats_data_cat(idx,:);
        datasets_dis(d,:) = for_stats_data_dis(idx,:);
    end
    
    avg_datasets_cat = squeeze(mean(datasets_cat,1));
    avg_datasets_dis = squeeze(mean(datasets_dis,1));
    
    peak_latency_all_samples(bs,1) = (find(avg_datasets_cat==max(avg_datasets_cat),1)-40)*5;
    peak_latency_all_samples(bs,2) = (find(avg_datasets_dis==max(avg_datasets_dis),1)-40)*5;
    
end

%% 2) Calculate difference between peaks for bootstrap samples 
peak_latency_difference = abs(squeeze(peak_latency_all_samples(:,1))-squeeze(peak_latency_all_samples(:,2)));

%% 3) Get 95% confidence interval for peak difference
bootstrap_peak_difference.CI_diff_cat_vs_dis = NaN(2,1);
bootstrap_peak_difference.CI_diff_cat_vs_dis(1) = prctile(peak_latency_difference,2.5);
bootstrap_peak_difference.CI_diff_cat_vs_dis(2) = prctile(peak_latency_difference,97.5);

%% Save the matrix
filename = sprintf('bootstrap_peak_difference_%s_subjects_%d_%d_10k_perms.mat',analysis,subjects(1),subjects(end));
save(fullfile(results_avg_dir,filename),'bootstrap_peak_difference'); 

end