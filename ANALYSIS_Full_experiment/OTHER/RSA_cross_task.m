function RSA_cross_task(subjects,with_stats) 
%RSA_CROSS_TASK Perform representational similarity analysis on SVM object decoding from both tasks. 
%
%Input: subject IDs, with stats (1) or without (0)
%
%Output: 
% 1)NxNxP representational dissimilarity matrices (1-Pearson's
% coefficient), one for each task
% 2)PxP matrix of 1-Spearman's correlations, where P is the number of timepoints. 
% 3)RSA time-time plot based on 2)
% 4)If with stats, NxN plot of statistical significance
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate
numConditions = 60;
numTimepoints = 200;
rdm_all_subjects_cat = NaN(max(subjects),numConditions,numConditions,numTimepoints);
rdm_all_subjects_dis = NaN(max(subjects),numConditions,numConditions,numTimepoints);

%% Loop: collect results from all subjects 
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    cat_filename = 'rdm_pearson_categorization.mat'; %svm_decoding_accuracy
    dis_filename = 'rdm_pearson_fixation.mat'; %svm_decoding_accuracy

    load(fullfile(subject_results_dir,cat_filename));
    rdm_all_subjects_cat(subject,:,:,:) = rdm_avg;
    load(fullfile(subject_results_dir,dis_filename));
    rdm_all_subjects_dis(subject,:,:,:) = rdm_avg;
end   

%% Average over subjects and if method 2, perform RSA
rdm_cat = squeeze(nanmean(rdm_all_subjects_cat,1)); %avg over subjects
rdm_dis = squeeze(nanmean(rdm_all_subjects_dis,1)); 
[rdm_rsa,rdm_cat_flattened,rdm_dis_flattened] = representational_SA(rdm_cat,rdm_dis,numTimepoints);

%% Plot
h = pcolor(rdm_rsa); 
set(gcf, 'Position', get(0, 'Screensize'));
set(h, 'EdgeColor', 'none');
axis square;
cbar = colorbar;
ylabel(cbar,'1-Spearman''s coefficient');
xlabel('Timepoints: Distraction task');
ylabel('Timepoints: Categorization task');
title(sprintf('Time-generalized RSA of scene processing in categorization and distraction tasks (N=%d)',numel(subjects)));

% Save RDMs, RSA results and figures
%Matrices
save(fullfile(results_avg_dir,sprintf('average_rdm_categorization_method_%d_subjects_%d_%d',method,subjects(1),subjects(end))),'rdm_cat');
save(fullfile(results_avg_dir,sprintf('average_rdm_fixation_method_%d_subjects_%d_%d',method,subjects(1),subjects(end))),'rdm_dis');
save(fullfile(results_avg_dir,sprintf('rsa_cross_task_method_%d_subjects_%d_%d',method,subjects(1),subjects(end))),'rdm_rsa');

%Figures
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_method_%d_subjects_%d_%d',method,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_method_%d_subjects_%d_%d.png',method,subjects(1),subjects(end)))); %save as png
close(gcf);

%% Plot stats if needed
if with_stats
    num_perms = 1000;
    
    %load the stats if they have been run already
    filename_sign = fullfile(results_avg_dir,sprintf('permutation_stats_rsa_cross_task_subjects_%d_%d.mat',...
        subjects(1),subjects(end)));
    if exist(filename_sign,'file')
        load(filename_sign);
        significant_timepoints = permutation_stats.SignificantMaxClusterSize; 
    else
        [~,significant_timepoints] = rsa_weighted_cluster_perm_stats(subjects,rdm_cat_flattened,rdm_dis_flattened,rdm_rsa,'left',1,num_perms);
    end
    
    %plot the stats - different plot from the RSA
    h = pcolor(significant_timepoints); 
    set(gcf, 'Position', get(0, 'Screensize'));
    set(h, 'EdgeColor', 'none');
    axis square;
    cbar = colorbar;
    ylabel(cbar,'Significance');
    xlabel('Timepoints: Distraction task');
    ylabel('Timepoints: Categorization task');
    title(sprintf('Statistical analysis of time-generalized RSA of scene processing in categorization and distraction tasks (N=%d)',numel(subjects)));  
    
    %save the plot
    saveas(gcf,fullfile(results_avg_dir,sprintf('stats_rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_avg_dir,sprintf('stats_rsa_cross_task_subjects_%d_%d.png',subjects(1),subjects(end)))); 
    close(gcf);
end 

end
    
