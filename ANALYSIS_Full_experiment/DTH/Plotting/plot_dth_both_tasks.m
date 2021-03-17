function plot_dth_both_tasks(subjects) %distance art, distance nat, RT
%PLOT_DTH_BOTH_TASKS Plot the fixed effects distance-to-hyperplane results for both
%tasks.
%
%Input: subjects' ID (e.g., 1:13)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%Load the correlations
load(fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_rt_both_categories_%d_subjects_categorization.mat',numel(subjects))));
correlation_both_categorization = correlation_dth_rt_both;
load(fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_rt_both_categories_%d_subjects_fixation.mat',numel(subjects))));
correlation_both_fixation = correlation_dth_rt_both;

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
plot(correlation_both_categorization,'LineWidth',2);
hold on;
plot(correlation_both_fixation,'LineWidth',2);
hold on;

%Plotting parameters
title =  sprintf('Correlation between the distance to hyperplane and reaction time in 60 scenes (N=%d)',numel(subjects));
legend_plot = {'Scene categorization','Distraction'};
xticks(0:10:200);
plotting_parameters(title,legend_plot,40,12,'best','Spearman''s coefficient'); %[0.7 0.85 0.1 0.01]

%% Save
saveas(gcf,fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_both_tasks',numel(subjects)))); 
saveas(gcf,fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_both_tasks.svg',numel(subjects))));

end