function rnn_dth_all_layers(with_stats)
%RNN_DTH_ALL_LAYERS Performs distance-to-hyperplane  using
%the distances-to-hyperplane using RNN dats. Plots all layers.
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
addpath(genpath('/home/agnek95/PyColormap4Matlab'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_avg_dir));

%% DTH
entropy_thresh = '0.09';

%Load RTs
load(sprintf('/scratch/agnek95/PDM/DATA/RNN_RTs/reaction_time_entropy_th_%s_model_22.08.mat',...
    entropy_thresh),'data');
RT = data';
    
%load distances
load('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS/rnn_distances_all_scenes_model_22.08.mat','data');
distances = data;
     
%% Correlate each subject's distances with the median RT
num_layers=size(distances,1);
num_timepoints_rnn=size(distances,2);
% all_conditions = [artificial_conditions natural_conditions];
correlation_both = NaN(num_layers,num_timepoints_rnn);
for l = 1:num_layers
    for tp = 1:num_timepoints_rnn
        correlation_both(l,:) = arrayfun(@(x) ...
            corr(squeeze(distances(l,x,:)),...
            RT,'type','Spearman'), 1:num_timepoints_rnn);
    end
end

%% Plot
figure;
% cmap_both = autumn;
colormap_plot = 'OrRd';
cmap = getPyPlot_cMap(colormap_plot,num_timepoints_rnn);

% legend_plot = cell(num_layers,1);
% colorplot = cmap_both(1:30:num_layers*30,:);

for l=1:num_layers
    plot(correlation_both(l,:),'LineWidth',2,'Color',cmap(l,:));
    hold on;
%     legend_plot{l} = sprintf('Layer %s',num2str(l));
end
% legend(legend_plot,'Location','best');
      
%% Plot stats if needed
if with_stats   
    %define some variables for the stats and the plot
    stats.num_perms = 1000;
    stats.tail = 'left';
    stats.alpha = 0.1;
    plot_location = -0.45:-0.05:-0.75;

    %Check if stats already exist, otherwise run them 
    filename_sign = 'rnn_all_layers_dth_stats';   
    filename = fullfile(results_avg_dir,sprintf('%s_both.mat',filename_sign));
    if exist('filename','file')
        load(filename,'stats');
    else
        corr_all = NaN(stats.num_perms,num_layers,num_timepoints_rnn);
        corr_all(1,:,:) = correlation_both; %ground truth
        for perm = 2:stats.num_perms
            for l = 1:num_layers
                permuted_RT = RT(randperm(numel(RT)));
                corr_all(perm,l,:) = arrayfun(@(x) corr(squeeze(distances(end,x,:)),permuted_RT,'type','Spearman'),1:num_timepoints_rnn);
            end
        end
        all_p_values = NaN(stats.num_perms,num_layers,num_timepoints_rnn);
        if strcmp(stats.tail,'left')
            corr_all_stat = corr_all*-1;
        else %two-tailed
            corr_all_stat = abs(corr_all);
        end
        for l = 1:num_layers
            for t=1:num_timepoints_rnn
                all_p_values(:,l,t) = (stats.num_perms+1 - tiedrank(corr_all_stat(:,l,t))) / stats.num_perms;
            end
        end
        stats.pvalue = squeeze(all_p_values(1,:,:)); %ground truth
        for l=1:num_layers
            stats.SignificantVariables(l,:) = stats.pvalue(l,:)<stats.alpha;
        end
        save(filename,'stats');
   end
   for l = 1:num_layers
        %Plot
        st = (stats.SignificantVariables(l,:)*plot_location(l)); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',cmap(l,:)); 
        hold on;
   end
%             [fdr_stats.SignificantVariables,fdr_stats.crit_p,~,fdr_stats.adjusted_pvalues] = fdr_bh(fdr_stats.pvalue,fdr_stats.alpha,'pdep');

end 

%% Plotting parameters
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);xlim([1,8]);
xticks(1:8);
ylim([-0.8,0.4]);
yticks(-0.8:0.2:0.4);
xlabel('Timepoint');
ylabel('Spearman''s r');

% legend_plot = {'Artificial scenes','Natural scenes','All scenes'}; 
% legend_bool = 0;
% title_bool = 0;
% plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Spearman''s r'); 

%% Save correlations and figures
dth_results.corr_both_categories = correlation_both;

save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
save(fullfile(save_path,'rnn_dth_all_layers.mat'),'dth_results');
saveas(gcf,fullfile(save_path,'rnn_dth_all_layers.fig'));
saveas(gcf,fullfile(save_path,'rnn_dth_all_layers.svg'));

close(gcf);
end









