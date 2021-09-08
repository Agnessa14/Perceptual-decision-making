function plot_RT_rnn_vs_human(subjects,with_stats)
%PLOT_RT_RNN_VS_HUMAN Create a scatterplot of reaction times.
%
%Returns a scatterplot of human reaction times in the categorization task vs
%RNN reaction times. 
%
%Author: Agnessa Karapetian, 2021
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load RTs 
numConditions = 60;
entropy_thresh = '0.09';

%RNN
load(sprintf('/scratch/agnek95/PDM/DATA/RNN_RTs/reaction_time_entropy_th_%s_model_22.08.mat',...
    entropy_thresh),'data');
RT_rnn_all = data';

%Human
RT_human_all = NaN(max(subjects),numConditions);

for subject = subjects    
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,'RTs_correct_trials_categorization.mat'),'RT_per_condition');
    RT_human_all(subject,:) = RT_per_condition; 
end
RT_human_all_med = nanmedian(RT_human_all,1)';


                        %%% Scatterplot: correlation rnn/human %%%
%% Correlation
corr_rt = corr(RT_human_all_med,RT_rnn_all,'type','Spearman');

%% Plot
num_subjects = numel(subjects);
figure;
color_plot = [0.05 0.6 0.13];
scatter(normalize(RT_human_all_med),normalize(RT_rnn_all),[],color_plot,'filled');
hold on
l = lsline;
l.Color = 'k';
l.LineWidth = 2;
lsline;
xlabel('Human RT (normalized)');
ylabel('RNN RT (normalized)');
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
        sprintf('stats_fdr_rnn_human_rt_correlation_subjects_%d_%d.mat',subjects(1),subjects(end)));
    if exist(filename,'file')
        load(filename,'stats_behav');
    else
        corr_rt_all = NaN(stats_behav.num_perms,1);
        corr_rt_all(1) = corr_rt; %ground truth
        for perm = 2:stats_behav.num_perms
            permuted_human_RT = RT_human_all_med(randperm(numel(RT_human_all_med)));
            corr_rt_all(perm) = corr(permuted_human_RT,RT_rnn_all,'type','Pearson');
        end
        
        all_p_values = (stats_behav.num_perms+1 - tiedrank(abs(corr_rt_all))) / stats_behav.num_perms;                 
        stats_behav.pvalue = all_p_values(1);
        save(filename,'stats_behav');
    end

    %check significance
    if stats_behav.pvalue < stats_behav.alpha
        text(5,3,['R = ',num2str(corr_rt),', p = ',num2str(stats_behav.pvalue)],'FontSize',font_size); 
    else
        text(5,3,['R = ',num2str(corr_rt)],'FontSize',font_size); 
    end

end

if ~with_stats
    text(5,3,['R = ',num2str(corr_rt)]); 
end

%% Save
saveas(gcf,fullfile(results_avg_dir,sprintf('rt_rnn_human_scatterplot_subjects_%d_%d',subjects(1),subjects(end))));
saveas(gcf,fullfile(results_avg_dir,sprintf('rt_rnn_human_scatterplot_subjects_%d_%d.svg',subjects(1),subjects(end))));
close(gcf);

end