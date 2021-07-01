function boxplot_RT_full_experiment(subjects,with_stats)
%BOXPLOT_RT_FULL_EXPERIMENT Create a boxplot of reaction times.
%
%Returns a scatterplot of reaction times in the categorization task vs
%reaction times in the distraction task, for each participant. 
%Different colours represent different participants and different shapes
%represent different categories (artificial vs natural).
%
%Author: Agnessa Karapetian, 2021
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
corr_random_effects = NaN(max(subjects),1);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,'RTs_correct_trials_categorization.mat'));
    artRTCat_all(subject,:) = RT_per_condition(artificial_conditions);
    natRTCat_all(subject,:) = RT_per_condition(natural_conditions);
    cat_rt_all_sub = RT_per_condition;
    load(fullfile(results_dir,subname,'RTs_correct_trials_fixation.mat'));
    artRTDis_all(subject,:) = RT_per_condition(artificial_conditions);
    natRTDis_all(subject,:) = RT_per_condition(natural_conditions);
    dis_rt_all_sub = RT_per_condition;
    corr_random_effects(subject) = corr(cat_rt_all_sub,dis_rt_all_sub,'type','Pearson');
end

%Average over subjects
artRTCat = nanmean(artRTCat_all,1)*1000; %in ms
natRTCat = nanmean(natRTCat_all,1)*1000;

artRTDis = nanmean(artRTDis_all,1)*1000; 
natRTDis = nanmean(natRTDis_all,1)*1000;

%% Boxplot
cat_RT_all = [artRTCat natRTCat]';
dis_RT_all = [artRTDis natRTDis]';
data = [artRTCat natRTCat artRTDis natRTDis];
groups = [ones(30,1); 2*ones(30,1); 3*ones(30,1); 4*ones(30,1)];
labels = {'Artificial', 'Natural', 'Artificial', 'Natural'};
boxplot(data,groups,'Labels',labels);
% xlabel('RT (ms)'); %figure out how to put xlabel 'Categorization' and
% 'Distraction'
ylabel('RT (ms)');
title_bool = 0;
if title_bool == 1
    title(sprintf('Reaction times for each scene across subjects (N=%d)',num_subjects));
end
set(gca,'FontName','Arial');

%% Stats
if with_stats
    stats_behav.num_perms = 1000;
    stats_behav.alpha = 0.05;
    stats_behav.tail = 'both';
    filename = fullfile(results_avg_dir,...
        sprintf('stats_fdr_rt_correlation_tasks_subjects_%d_%d.mat',subjects(1),subjects(end)));
    if exist('filename','file')
        load(filename,'stats_behav');
    else
        corr_rt_all = NaN(stats_behav.num_perms,1);
        corr_rt_all(1) = corr_rt; %ground truth
        for perm = 2:stats_behav.num_perms
            permuted_cat_RT = cat_RT_all(randperm(numel(cat_RT_all)));
            corr_rt_all(perm) = corr(permuted_cat_RT,dis_RT_all,'type','Pearson');
        end
        
        all_p_values = (stats_behav.num_perms+1 - tiedrank(abs(corr_rt_all))) / stats_behav.num_perms;                 
        stats_behav.pvalue = all_p_values(1);
        save(filename,'stats_behav');
    end

end


%% Save
saveas(gcf,fullfile(results_avg_dir,sprintf('boxplot_cat_vs_fix_subjects_%d_%d',subjects(1),subjects(end))));
saveas(gcf,fullfile(results_avg_dir,sprintf('boxplot_cat_vs_fix_subjects_%d_%d.svg',subjects(1),subjects(end))));
close(gcf);

end