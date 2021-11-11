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
entropy_thresh = '0.07';
model_name = 'model_02.11_2';

%RNN
load(sprintf('/scratch/agnek95/PDM/DATA/RNN_RTs/reaction_time_entropy_th_%s_%s.mat',...
    entropy_thresh,model_name),'data');
RT_rnn_all = data';

%Human
% RT_human_all = NaN(max(subjects),numConditions);
load(fullfile(results_avg_dir,'RTs_shuffled_across_subjects.mat'),'data');
RTs_eeg = data;

%Take second half of all subjects (shuffled) & get median
num_subjects = numel(subjects);
RT_human_all_med = nanmedian(RTs_eeg((num_subjects/2)+1:end,:),1)';


                        %%% Scatterplot: correlation rnn/human %%%
%% Correlation
corr_rt = corr(RT_human_all_med,RT_rnn_all,'type','Spearman');

%% Plot
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
xlim([-3 3]);

%% Stats
if with_stats
    rng('shuffle');
    stats_behav.num_perms = 1000;
    stats_behav.alpha = 0.05;
    stats_behav.tail = 'both';
    filename = fullfile(results_avg_dir,...
        sprintf('cv_%s_stats_fdr_rnn_human_rt_correlation_subjects_%d_%d.mat',model_name,subjects(1),subjects(end)));
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
        text(-3,1,['R = ',num2str(corr_rt),', p = ',num2str(stats_behav.pvalue)],'FontSize',font_size); 
    else
        text(-3,1,['R = ',num2str(corr_rt)],'FontSize',font_size); 
    end

end

if ~with_stats
    text(5,3,['R = ',num2str(corr_rt)]); 
end

%% Save
if strcmp(model_name,'model_02.11')
    savename = '1';
elseif strcmp(model_name,'model_02.11_2')
    savename = '2';
else
    error('Unknown model name');
end
saveas(gcf,fullfile(results_avg_dir,sprintf('cv_%s_rt_rnn_human_scatterplot_subjects_%d_%d',savename,subjects(1),subjects(end))));
saveas(gcf,fullfile(results_avg_dir,sprintf('cv_%s_rt_rnn_human_scatterplot_subjects_%d_%d.svg',savename,subjects(1),subjects(end))));
close(gcf);

end