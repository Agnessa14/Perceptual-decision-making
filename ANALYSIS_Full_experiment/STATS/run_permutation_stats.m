function significant_timepoints = run_permutation_stats(subjects,task,analysis,for_stats,nperm,cluster_threshold,significance_threshold) 
%RUN_PERMUTATION_STATS Call on the permutation stats script to perform
%cluster based permutation tests. 
%
%Input: subject IDs, task (1=categorization,2=distraction,3=cross-task), analysis
%('object_decoding', 'category_decoding', 'time_object_decoding', 'time_category_decoding' or 'rsa_time_object')
%
%Output: SxP matrix containing decoding accuracies, where S is the number
%of subjects and P is the number of timepoints.


%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
if task < 3
    task_name = get_task_name(task);
elseif task == 3
    task_name = 'cross_task';
end

%% Parameters for the stats script  
tail = 'right'; 

%% Run the stats script 
[significant_timepoints, pvalues] = permutation_cluster_1sample_alld(for_stats, nperm, cluster_threshold, significance_threshold,tail); 

%% Save 
save(fullfile(results_avg_dir,sprintf('significant_timepoints_subjects_%d_%d_%s_task_%s',...
    subjects(1),subjects(end),task_name,analysis)),'significant_timepoints');
save(fullfile(results_avg_dir,sprintf('pvalues_subjects_%d_%d_%s_task_%s',...
    subjects(1),subjects(end),task_name,analysis)),'pvalues');
