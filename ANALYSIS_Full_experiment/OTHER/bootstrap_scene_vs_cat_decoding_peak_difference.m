function bootstrap_scene_vs_cat_decoding_peak_difference(subjects,task)
%BOOTSTRAP_DECODING_PEAK_DIFFERENCE Calculate the confidence intervals of the difference between peak
%latencies of categorization and distraction decoding.
%
%Input: subject IDs, task (1 for categorization, 2 for distraction)
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
task_name = get_task_name(task);
numConditions = 60;
numTimepoints = 200;
scene_decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints);
cat_decoding_accuracies_all_subjects = NaN(max(subjects),numTimepoints);

%% Loop: collect results from all subjects
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    
    scene_dec_filename = sprintf('svm_decoding_accuracy_%s.mat',task_name);
    cat_dec_filename = sprintf('cross_validated_dth_pseudotrials_svm_decodingAccuracy_%s.mat',task_name);
    all_dimensions_scene = repmat({':'},1,3);  %conditions x conditions x timepoints
    all_dimensions_cat = {':'}; % timepoints

    load(fullfile(subject_results_dir,scene_dec_filename),'decodingAccuracy_avg');
    scene_decoding_accuracies_all_subjects(subject,all_dimensions_scene{:}) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,cat_dec_filename),'decodingAccuracy_avg');
    cat_decoding_accuracies_all_subjects(subject,all_dimensions_cat{:}) = decodingAccuracy_avg;
end

%% For stats matrix: num subjects x num timepoints (average over conditions)
for_stats_scene = squeeze(nanmean(nanmean(scene_decoding_accuracies_all_subjects,2),3));
for_stats_cat = cat_decoding_accuracies_all_subjects;
for_stats_data_scene = for_stats_scene(subjects,:)-50;
for_stats_data_cat = for_stats_cat(subjects,:)-50;

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
    
    datasets_scene = NaN(size(for_stats_data_scene));
    datasets_cat = NaN(size(for_stats_data_cat));
    
    for d = 1:num_datasets
        idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
        datasets_scene(d,:) = for_stats_data_scene(idx,:);
        datasets_cat(d,:) = for_stats_data_cat(idx,:);
    end
    
    avg_datasets_scene = squeeze(mean(datasets_scene,1));
    avg_datasets_cat = squeeze(mean(datasets_cat,1));
    
    peak_latency_all_samples(bs,1) = (find(avg_datasets_scene==max(avg_datasets_scene),1)-40)*5;
    peak_latency_all_samples(bs,2) = (find(avg_datasets_cat==max(avg_datasets_cat),1)-40)*5;
    
end

%% 2) Calculate difference between peaks for bootstrap samples 
peak_latency_difference = abs(squeeze(peak_latency_all_samples(:,1))-squeeze(peak_latency_all_samples(:,2)));

%% 3) Get 95% confidence interval for peak difference
bootstrap_peak_difference.CI_diff_scene_vs_cat = NaN(2,1);
bootstrap_peak_difference.CI_diff_scene_vs_cat(1) = prctile(peak_latency_difference,2.5);
bootstrap_peak_difference.CI_diff_scene_vs_cat(2) = prctile(peak_latency_difference,97.5);

%% Save the matrix
filename = sprintf('bootstrap_scene_vs_cat_peak_difference_%s_subjects_%d_%d.mat',task_name,subjects(1),subjects(end));
save(fullfile(results_avg_dir,filename),'bootstrap_peak_difference'); 

end