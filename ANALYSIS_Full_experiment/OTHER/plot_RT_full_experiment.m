function plot_RT_full_experiment(subjects,with_stats)
%PLOT_RT_FULL_EXPERIMENT Create a scatterplot of reaction times.
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

                        %%% Scatterplot: correlation dis/cat %%%
%% Correlation
% corr_random_avg = squeeze(nanmean(corr_random_effects,1));

cat_RT_all = [artRTCat natRTCat]';
dis_RT_all = [artRTDis natRTDis]';

corr_rt = corr(cat_RT_all,dis_RT_all,'type','Pearson');

%% Plot
% num_subjects = numel(subjects);
% figure;
% set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
% c = linspace(1,numel(artificial_conditions),numel(natural_conditions));
% scatter(artRTCat,artRTDis,150,c,'s','filled');
% hold on;
% scatter(natRTCat,natRTDis,150,c,'filled');
% legend({'Artificial scenes','Natural scenes'},'FontSize',16,'Location','Best');
% b = colorbar;
% ylabel(b,'Scene ID','FontSize',12);
% xlim([430 530]);
% ylim([430 530]);
% xlabel('RT categorization task (ms)');
% ylabel('RT distraction task (ms)');
% title(sprintf('Reaction times for each scene across subjects (N=%d)',num_subjects),'FontSize',16);

num_subjects = numel(subjects);
figure;
color_plot = [0.05 0.6 0.13];
scatter(cat_RT_all,dis_RT_all,[],color_plot,'filled');
hold on
xlim([470 515]);
l = lsline;
l.Color = 'k';
l.LineWidth = 2;
lsline;
xlabel('RT categorization task (ms)');
ylabel('RT distraction task (ms)');
title_bool = 0;
if title_bool == 1
    title(sprintf('Reaction times for each scene across subjects (N=%d)',num_subjects));
end
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);

%% Stats
if with_stats
    rng('shuffle');
    stats_behav.num_perms = 1000;
    stats_behav.alpha = 0.05;
    stats_behav.tail = 'both';
    filename = fullfile(results_avg_dir,...
        sprintf('stats_fdr_rt_correlation_tasks_subjects_%d_%d.mat',subjects(1),subjects(end)));
    if exist(filename,'file')
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

    %check significance
    if stats_behav.pvalue < stats_behav.alpha
        text(472,467,['R = ',num2str(corr_rt),', p = ',num2str(stats_behav.pvalue)],'FontSize',font_size); 
    else
        text(472,467,['R = ',num2str(corr_rt)],'FontSize',font_size); 
    end

end

if ~with_stats
    text(472,467,['R = ',num2str(corr_rt)]); 
end

%% Save
saveas(gcf,fullfile(results_avg_dir,sprintf('scatter_RT_plot_cat_vs_fix_subjects_%d_%d',subjects(1),subjects(end))));
saveas(gcf,fullfile(results_avg_dir,sprintf('scatter_RT_plot_cat_vs_fix_subjects_%d_%d.svg',subjects(1),subjects(end))));
close(gcf);

end