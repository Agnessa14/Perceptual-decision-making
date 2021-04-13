function significant_timepoints = run_permutation_stats(subjects,task,analysis) 
%RUN_PERMUTATION_STATS Call on the permutation stats script to perform
%cluster based permutation tests. 
%
%Input: subject IDs, task (1=categorization,2=distraction), analysis
%('object_decoding' or 'category_decoding')
%
%Output: SxP matrix containing decoding accuracies, where S is the number
%of subjects and P is the number of timepoints.


%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
task_name = get_task_name(task);

%% Create the for_stats matrix: decoding accuracies for each subject and each timepoint
for_stats = all_subjects_for_stats(subjects,task,analysis);
disp('Creating for_stats_matrix');

%% Parameters for the stats script  
nperm = 10000;
cluster_threshold = 0.05;
significance_threshold = 0.05;
tail = 'right'; 

%% Run the stats script 
[significant_timepoints, pvalues] = permutation_cluster_1sample_alld(for_stats, nperm, cluster_threshold, significance_threshold,tail); 

%% Save 
save(fullfile(results_avg_dir,sprintf('significant_timepoints_%d_subjects_%s_task_%s',...
    numel(subjects),task_name,analysis)),'significant_timepoints');
save(fullfile(results_avg_dir,sprintf('pvalues_%d_subjects_%s_task_%s',...
    numel(subjects),task_name,analysis)),'pvalues');
