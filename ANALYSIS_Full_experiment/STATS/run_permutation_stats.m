function [significant_timepoints,pvalues] = run_permutation_stats(for_stats,nperm,cluster_threshold,significance_threshold,tail) 

%RUN_PERMUTATION_STATS Call on the permutation stats script to perform
%cluster based permutation tests. 
%
%Input: subject IDs, task (1=categorization,2=distraction,3=cross-task), analysis
%('object_decoding', 'category_decoding', 'time_object_decoding', 'time_category_decoding' or 'rsa_time_object')
%
%Output: SxP matrix containing decoding accuracies, where S is the number
%of subjects and P is the number of timepoints


%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));


%% Run the stats script 
[significant_timepoints, pvalues] = permutation_cluster_1sample_alld(for_stats, nperm, cluster_threshold, significance_threshold,tail); 

