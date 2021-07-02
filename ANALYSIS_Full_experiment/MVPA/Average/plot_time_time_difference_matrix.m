function plot_time_time_difference_matrix(subjects,analysis,with_stats)
%PLOT_TIME_TIME_DIFFERENCE_MATRIX Plot the difference between categorization & distraction within-task 
%time-generalization, averaged over all participants.
%
%Input: subject IDs (e.g., 1:13), analysis('object_decoding' or
%'category_decoding'), with or without stats (1/0)
%
%Output: 2D heatmap of difference between decoding accuracies. 
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load the data
numTimepoints = 50;
filename_for_stats_cat = fullfile(results_avg_dir,sprintf('for_stats_timegen_svm_%s_subjects_%d_%d_categorization.mat',analysis,...
    subjects(1),subjects(end)));

filename_for_stats_fix = fullfile(results_avg_dir,sprintf('for_stats_timegen_svm_%s_subjects_%d_%d_fixation.mat',analysis,...
    subjects(1),subjects(end)));

if exist(filename_for_stats_cat,'file') && exist(filename_for_stats_fix,'file')
    load(filename_for_stats_cat,'for_stats');
    for_stats_cat = for_stats;
    load(filename_for_stats_fix,'for_stats');
    for_stats_fix = for_stats;
else
    error('Run ''plot_timetime_decoding_full_experiment'' first.')
end

%% Plot the difference
for_stats_diff = for_stats_cat - for_stats_fix;
diff_matrix = squeeze(nanmean(for_stats_diff,1));
h = pcolor(diff_matrix); 
set(h, 'EdgeColor', 'none');
axis square;
hold on;
plot_title = sprintf('Difference in time-generalized %s between categorization and distraction tasks (N=%d)',analysis,numel(subjects));
title_bool = 0;
if title_bool==1
    title(plot_title);
end
cbar = colorbar;
ylabel('Timepoints trained on');
xlabel('Timepoints tested on');
ylabel(cbar,'%');
caxis([-5 20]);
xticks(0:5:50);
xticklabels(-200:100:800);
yticklabels(-200:100:800);
yticks(0:5:50);
xline(10,'--','Color','w');
yline(10,'--','Color','w');

set(gca,'FontName','Arial','FontSize',11)
%Plot a black line over the diagonal for cross-task
plot(1:numTimepoints,1:numTimepoints,'LineWidth',2.5,'Color','k');

if with_stats

    stats_decoding.num_perms = 1000;
    stats_decoding.qvalue = 0.01;
    stats_decoding.tail = 'right';
    filename = fullfile(results_avg_dir,...
        sprintf('stats_fdr_timetime_diff_%s_subjects_%d_%d.mat',analysis,subjects(1),subjects(end)));
    if exist(filename,'file')
        load(filename,'stats_decoding');
    else
        [stats_decoding.significant_timepoints,stats_decoding.pvalues,...
            stats_decoding.crit_p,stats_decoding.adjusted_pvalues]...
            = fdr_permutation_cluster_1sample_alld(for_stats-50,...
            stats_decoding.num_perms,stats_decoding.tail,stats_decoding.qvalue);
        save(filename,'stats_decoding');
    end
    
    %plot contour
    contour(stats_decoding.significant_timepoints,1,'LineColor','w','LineWidth',2);

end

%% Save the matrix plot
save(fullfile(results_avg_dir,sprintf('diff_matrix_timegen_svm_%s_subjects_%d_%d',analysis,subjects(1),subjects(end))),'diff_matrix');
saveas(gcf,fullfile(results_avg_dir,sprintf('diff_matrix_timegen_svm_%s_subjects_%d_%d',analysis,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('diff_matrix_timegen_svm_%s_subjects_%d_%d.png',analysis,subjects(1),subjects(end)))); %save as png
close(gcf)



