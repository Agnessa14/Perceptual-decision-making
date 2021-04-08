function plot_dth_both_tasks(subjects,with_cross_task) %distance art, distance nat, RT
%PLOT_DTH_BOTH_TASKS Plot the fixed effects distance-to-hyperplane results for both
%tasks.
%
%Input: subjects' ID (e.g., 1:13), add the cross task curves (1 yes, 0 no)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%Load the correlations
load(fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_categorization.mat',numel(subjects))));
corr_both_categorization = dth_results.corr_both_categories;
clear dth_results;
load(fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_fixation.mat',numel(subjects))));
corr_both_distraction = dth_results.corr_both_categories;
clear dth_results;

if with_cross_task   
    load(fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_cross_task_categorization_distances.mat',numel(subjects))));
    corr_both_cross_task_categorization_dist = dth_results.corr_both_categories;
    clear dth_results;
    load(fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_cross_task_fixation_distances.mat',numel(subjects))));
    corr_both_cross_task_distraction_dist = dth_results.corr_both_categories;
    clear dth_results;
end

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
plot(corr_both_categorization,'LineWidth',2);
hold on;
plot(corr_both_distraction,'LineWidth',2);
hold on;

if with_cross_task
    plot(corr_both_cross_task_categorization_dist,'LineWidth',2)
    hold on;
    plot(corr_both_cross_task_distraction_dist,'LineWidth',2)
    hold on;
end

%Plotting parameters
plot_title =  sprintf('Correlation between the distance to hyperplane and reaction time in 60 scenes (N=%d)',numel(subjects));
if with_cross_task
    legend_plot = {'Scene categorization','Distraction',...
        'Cross-task: scene categorization distances with distraction RTs',...
        'Cross-task: distraction distances with categorization RTs'};
else
    legend_plot = {'Scene categorization','Distraction'};
end
xticks(0:10:200);
plotting_parameters(plot_title,legend_plot,40,12,'best','Spearman''s coefficient'); %[0.7 0.85 0.1 0.01]

%% Save
if with_cross_task
    saveas(gcf,fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_both_tasks_with_crosstask',numel(subjects)))); 
    saveas(gcf,fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_both_tasks_with_crosstask.svg',numel(subjects))));
else
    saveas(gcf,fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_both_tasks',numel(subjects)))); 
    saveas(gcf,fullfile(results_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_both_tasks.svg',numel(subjects))));
end
end