function significantVarWeight = weighted_run_permutation_stats(subjects,task_distance,task_RT,analysis,category,for_stats,nperm,cluster_threshold,significance_threshold,tail) 
%WEIGHTED_RUN_PERMUTATION_STATS Call on the weighted permutation stats script to perform
%cluster based permutation tests. 
%
%Input: subject IDs, task_distance & task_RT (1=categorization,2=distraction,3=cross-task), analysis
%('random_dth'), category
%('average','both','artificial','natural'),for_stats (NxP matrix of DTH-RT
%correlations, N being the number of subjects), nperm (number of
%permutations), cluster_threshold & signficance_threshold, tail ('left' or
%'both') for one-tailed or two-tailed test
%
%Output: 1xP vector of significances (1 or 0), where P is the number of
%timepoints.
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
task_distance_name = get_task_name(task_distance);

%% Run the stats script 
[significantVarWeight,pValWeight,significantVarSize,pValSize] = permutation_cluster_1sample_weight_alld(for_stats,nperm,cluster_threshold,significance_threshold,tail); 

%% Save 
permutation_stats.SignificantMaxClusterSize = significantVarSize;
permutation_stats.SignificantMaxClusterWeight = significantVarWeight;
permutation_stats.pValueClusterSize = pValSize;
permutation_stats.pValueWeight = pValWeight;
permutation_stats.info.num_permutations = nperm;
permutation_stats.info.cluster_th = cluster_threshold;
permutation_stats.info.significance_th = significance_threshold;
permutation_stats.info.tail = tail;

if task_distance==task_RT
    filename = 'dth_permutation_stats';
else
    filename = 'dth_permutation_stats_crosstask';
end

save(fullfile(results_avg_dir,sprintf('%s_%d_%d_%s_task_%s_%s',...
    filename,subjects(1),subjects(end),task_distance_name,analysis,category)),'permutation_stats');
end
