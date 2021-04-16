function plot_RT_full_experiment(subjects)
%PLOT_RT_FULL_EXPERIMENT Create a scatterplot of reaction times.
%
%Returns a scatterplot of reaction times in the categorization task vs
%reaction times in the distraction task, for each participant. 
%Different colours represent different participants and different shapes
%represent different categories (artificial vs natural).
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load RTs for each subject
artificial_conditions = 1:30;
natural_conditions = 31:60;
artRTCat_all = NaN(max(subjects),numel(artificial_conditions));
natRTCat_all = NaN(max(subjects),numel(natural_conditions));
artRTDis_all   = NaN(max(subjects),numel(artificial_conditions));
natRTDis_all   = NaN(max(subjects),numel(natural_conditions));

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,'RTs_correct_trials_categorization.mat'));
    artRTCat_all(subject,:) = RT_per_condition(artificial_conditions);
    natRTCat_all(subject,:) = RT_per_condition(natural_conditions);
    load(fullfile(results_dir,subname,'RTs_correct_trials_fixation.mat'));
    artRTDis_all(subject,:) = RT_per_condition(artificial_conditions);
    natRTDis_all(subject,:) = RT_per_condition(natural_conditions);
end

%Average over subjects
artRTCat = nanmean(artRTCat_all,1)*1000; %in ms
natRTCat = nanmean(natRTCat_all,1)*1000;

artRTDis = nanmean(artRTDis_all,1)*1000; 
natRTDis = nanmean(natRTDis_all,1)*1000;

%% Plot
num_subjects = numel(subjects);
figure;
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
c = linspace(1,numel(artificial_conditions),numel(natural_conditions));
scatter(artRTCat,artRTDis,150,c,'s','filled');
hold on;
scatter(natRTCat,natRTDis,150,c,'filled');
legend({'Artificial scenes','Natural scenes'},'FontSize',16,'Location','Best');
b = colorbar;
ylabel(b,'Scene ID','FontSize',12);
xlim([440 530]);
ylim([440 530]);
xlabel('RT categorization task (ms)');
ylabel('RT distraction task (ms)');
title(sprintf('Reaction times for each scene across subjects (N=%d)',num_subjects),'FontSize',20);

%% Save
saveas(gcf,fullfile(results_avg_dir,sprintf('RT_plot_cat_vs_fix_both_tasks_%d_subjects',num_subjects)));
saveas(gcf,fullfile(results_avg_dir,sprintf('RT_plot_cat_vs_fix_both_tasks_%d_subjects.svg',num_subjects)));
close(gcf);
end