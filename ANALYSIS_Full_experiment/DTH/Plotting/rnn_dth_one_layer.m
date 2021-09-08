function rnn_dth_one_layer(with_stats)
%RNN_DTH_ONE_LAYER Performs distance-to-hyperplane  using
%the distances-to-hyperplane using RNN dats. Plots the last layer.
%
%Input: with stats (1 or 0)
%
%Saves the plot and correlation values.
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_avg_dir));

%% DTH
numConditions = 60;
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;  
entropy_thresh = '0.09';

%Load RTs
load(sprintf('/scratch/agnek95/PDM/DATA/RNN_RTs/reaction_time_entropy_th_%s_model_22.08.mat',...
    entropy_thresh),'data');
RT = data';
RT_art = RT(artificial_conditions);
RT_nat = RT(natural_conditions);
    
%load distances
load('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS/rnn_distances_all_scenes_model_22.08.mat','data');
distances = data;
     
%% Correlate each subject's distances with the median RT
% num_layers=size(distances,1);
num_timepoints_rnn=size(distances,2);
all_conditions = [artificial_conditions natural_conditions];

%only take the last layer
for tp = 1:num_timepoints_rnn
    correlation_art = arrayfun(@(x) ...
        corr(squeeze(distances(end,x,artificial_conditions)),...
        RT_art,'type','Spearman'), 1:num_timepoints_rnn);
    correlation_nat = arrayfun(@(x) ...
        corr(squeeze(distances(end,x,natural_conditions)),...
        RT_nat,'type','Spearman'), 1:num_timepoints_rnn);
    correlation_both = arrayfun(@(x) ...
        corr(squeeze(distances(end,x,all_conditions)),...
        RT,'type','Spearman'), 1:num_timepoints_rnn);
end

%% Plot
figure(abs(round(randn*10)));
cmap_1 = cool;
cmap_2 = summer;
color_art = cmap_1(200,:); %purple
color_nat = cmap_2(100,:); %green
color_both = 'k';
plot(correlation_art,'LineWidth',2,'Color',color_art);
hold on;
plot(correlation_nat,'LineWidth',2,'Color',color_nat);
plot(correlation_both,'LineWidth',2,'Color',color_both);

%% Plot stats if needed
if with_stats   
    %define some variables for the stats and the plot
    permutation_stats.num_perms = 1000;
    permutation_stats.cluster_th = 0.05;
    permutation_stats.significance_th = 0.05;
    permutation_stats.tail = 'left';
    for c = 1:3
        if c == 1
            category = 'artificial';            
            plot_location = -0.7;
            color = color_art;
        elseif c == 2
            category = 'natural';            
            plot_location = -0.75;
            color = color_nat;
        elseif c == 3
            category = 'both'; 
            plot_location = -0.8;
        end

        %Check if stats already exist, otherwise run the stats script
        filename_sign = 'rnn_layer7_dth_permutation_stats';   
        filename = fullfile(results_avg_dir,sprintf('%s_%s.mat',filename_sign,category));
        if exist(filename,'file')
            load(filename,'permutation_stats');
        else
            [permutation_stats.SignificantMaxClusterWeight,permutation_stats.pValWeight,...
                permutation_stats.SignificantMaxClusterSize,permutation_stats.pValSize] = ...
                permutation_cluster_1sample_weight_alld(for_stats,permutation_stats.num_perms,...
                permutation_stats.cluster_th,permutation_stats.significance_th,permutation_stats.tail); 
        save(filename,'permutation_stats');
        end
        
        %Plot
        if c < 3
            st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
            st(st==0) = NaN;
            plot(st,'*','Color',color); 
            hold on;
        end
    end
end 

%% Plotting parameters
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
legend_plot = {'Artificial scenes','Natural scenes','All scenes'}; 
% ylim([-0.5 0.3]);
legend_bool = 0;
title_bool = 0;
plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Spearman''s r'); 

%% Save correlations and figures
dth_results.corr_both_categories = avg_corr_both;
dth_results.corr_artificial = avg_corr_art;
dth_results.corr_natural = avg_corr_nat;

save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
save(fullfile(save_path,'rnn_dth_layer7.mat'),'dth_results');
saveas(gcf,fullfile(save_path,'rnn_dth_layer7.fig'));
saveas(gcf,fullfile(save_path,'rnn_dth_layer7.svg'));

close(gcf);
end









