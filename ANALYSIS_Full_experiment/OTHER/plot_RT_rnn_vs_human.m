function plot_RT_rnn_vs_human
%PLOT_RT_RNN_VS_HUMAN Create a scatterplot of reaction times.
%
%Returns a scatterplot of human reaction times in the categorization task vs
%RNN reaction times. 
%
%Input: subject IDs (e.g., 1:13), conditions ('all', 'artificial' or
%'natural'), with or without stats (1/0)
%
%Author: Agnessa Karapetian, 2021
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load the correlations and the RTs for the noise ceiling calculation
load(fullfile(results_avg_dir,'02.11_2_rnn/Model_RDM_redone','correlation_RT_human_RNN_cross-validated.mat'),'data');
load(fullfile(results_avg_dir,'RT_all_subjects_5_35_categorization.mat'),'RTs');
correlations_RT = data;
num_subjects = size(correlations_RT,1);
RT_noise_ceiling_lower = NaN(3,1);
xticklabels = cell(3,1);
h = NaN(3,3);
location_x = NaN(3,1);
label_x = NaN(3,1);
figure;

for c = 1:3
    if c == 1
        xticklabels{c} = 'All scenes';
        conds = 1:60;
        location_x(c) = 0.5;
        label_x(c) = 0.5;
    elseif c==2
        xticklabels{c} = 'Man-made scenes';
        conds = 1:30;
        location_x(c) = 1.5;
        label_x(3) = 1.5;
    elseif c==3
        xticklabels{c} = 'Natural scenes';
        conds = 31:60;
        location_x(c) = 1;
        label_x(2) = 1;
    end

    corr_conditions = correlations_RT(:,c);
    average_corr = mean(corr_conditions,1);
    
    %Calculate noise ceiling
    noise_ceil_temp = NaN(1,num_subjects);
    for subject = 1:num_subjects
        conds_bool = ~isnan(RTs(subject,conds));
        noise_ceil_temp(subject) = corr(nanmean(RTs(1:end~=subject,conds(conds_bool)),1)',RTs(subject,conds(conds_bool))','type','Pearson');
    end
    RT_noise_ceiling_lower(c) = mean(noise_ceil_temp);
    
    %Plot all subjects, average and noise ceiling
    color_subjects = [0.75 0.75 0.75];   
    color_avg = [0 0.5 0.5];
    color_noise_ceiling = 'k';
    x = repmat(location_x(c),1,num_subjects);
    h(c,1) = scatter(x,corr_conditions,[],color_subjects,'filled','SizeData',100,'DisplayName','Single subject');
    hold on;
    h(c,2) = scatter(location_x(c),RT_noise_ceiling_lower(c),[],color_noise_ceiling,'filled','SizeData',100,'DisplayName','Noise ceiling');
    h(c,3) = scatter(location_x(c),average_corr,[],color_avg,'filled','SizeData',200,'DisplayName','Average');   
end

%% Plot parameters
% xlabel('Corr human-RNN RTs (Pearson''s r)');
% ylabel('Conditions');
xlim([0.2 3]);
legend(h(1,:),'Location','southeast');
% set(gca,'xtick',location_x,'xticklabel',xticklabels);
set(gca,'xtick',label_x,'xticklabel','');
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);

%% Save
%RT plots
saveas(gcf,fullfile(results_avg_dir,'corr_RT_human_RNN.fig'));
saveas(gcf,fullfile(results_avg_dir,'corr_RT_human_RNN.svg'));
close(gcf);

end