function rnn_dth_distances_RNN(subjects,with_stats)
%RNN_DTH_DISTANCES_RNN Performs the cross-modal distance-to-hyperplane analysis using
%the distances-to-hyperplane from RNN and the median RT across
%subjects.
%
%Input: subjects' ID (e.g., 1:13), add stats (1 with, 0 without)
%
%Correlates the decision values of each subject with reaction times averaged over
%participants of each condition (60 scenes), averaged, at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%% Get the distances 
numTimepoints = 200;
numConditions = 60;
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;  
entropy_thresh = '0.02';
model_name = 'model_02.11_2';

%load distances
load(sprintf('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS/rnn_distances_all_scenes_%.mat',model_name),'data');
distances = data;

%load RTs & get median
RTs = NaN(max(subjects),numConditions);
RTs_art = NaN(max(subjects),numConditions/2);
RTs_nat = NaN(size(RTs_art));
% for subject = subjects    
%     subname = get_subject_name(subject);
%     load(fullfile(results_dir,subname,'RTs_correct_trials_categorization.mat'),'RT_per_condition');
%     RTs(subject,:) = RT_per_condition; 
%     RTs_art(subject,:) = RT_per_condition(artificial_conditions);
%     RTs_nat(subject,:) = RT_per_condition(natural_conditions);
% end
load(fullfile(results_avg_dir,'RT_all_subjects_5_35_categorization.mat'),'RTs');
RT = nanmedian(RTs,1);
RT_art = nanmedian(RTs_art,1);
RT_nat = nanmedian(RTs_nat,1);

%% Correlate each subject's distances with the median RT
t = 1:numTimepoints;
size_corr = [max(subjects),numTimepoints];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);

%Correlate over all layers
num_layers=size(distances,1);
num_timepoints_rnn=size(distances,2);
size_corr = [num_layers,num_timepoints_rnn];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);
all_conditions = [artificial_conditions natural_conditions];

for l = 1:num_layers
    for tp = 1:num_timepoints_rnn
        correlation_art(l,:) = arrayfun(@(x) ...
            corr(squeeze(distances(l,x,artificial_conditions)),...
            RT_art','type','Spearman'), 1:num_timepoints_rnn);
        correlation_nat(l,:) = arrayfun(@(x) ...
            corr(squeeze(distances(l,x,natural_conditions)),...
            RT_nat','type','Spearman'), 1:num_timepoints_rnn);
        correlation_both(l,:) = arrayfun(@(x) ...
            corr(squeeze(distances(l,x,all_conditions)),...
            RT','type','Spearman'), 1:num_timepoints_rnn);
    end
end


%% Plot
cmap_1 = cool;
cmap_2 = summer;
cmap_3 = autumn;

for c = 1:3
    figure;
    legend_plot = cell(num_layers,1);

    if c==1
        corrplot = correlation_art;
        colorplot = cmap_1(1:30:num_layers*30,:);
    elseif c==2
        corrplot = correlation_nat;
        colorplot = cmap_2(1:30:num_layers*30,:);
    elseif c==3
        corrplot = correlation_both;
        colorplot = cmap_3(1:30:num_layers*30,:);
    end
    for l=1:num_layers
        plot(corrplot(l,:),'LineWidth',2,'Color',colorplot(l,:,:));
        hold on;
        legend_plot{l} = sprintf('Layer %s',num2str(l));
    end
    legend(legend_plot,'Location','best');
end  


%% Plot stats if needed
if with_stats   
    %define some variables for the stats and the plot
    analysis = 'random_dth';
    permutation_stats.num_perms = 1000;
    permutation_stats.cluster_th = 0.05;
    permutation_stats.significance_th = 0.05;
    permutation_stats.tail = 'left';
    for c = 1:3
        if c == 1
            category = 'artificial';            
            plot_location = -0.22;
            color = color_art;
            for_stats = correlation_art(subjects,:);
        elseif c == 2
            category = 'natural';            
            plot_location = -0.24;
            color = color_nat;
            for_stats = correlation_nat(subjects,:);
        elseif c == 3
            category = 'both'; 
            plot_location = -0.26;
            color = 'k';
            for_stats = correlation_both(subjects,:);
        end

        %Check if stats already exist, otherwise run the stats script
        distances_str = 'rnn';
        filename_sign = sprintf('cv_%s_rnn_dth_rnn_distances_permutation_stats_crosstask',model_name);        
        filename = fullfile(results_avg_dir,sprintf('%s_%d_%d_distances_%s_%s_%s.mat',filename_sign,...
            subjects(1),subjects(end),distances_str,analysis,category));
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
        st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;
    end
end 

%% Plotting parameters
% if task_distance==task_RT
%     if task_distance==1
%         task_title = 'scene categorization';
%     elseif task_distance==2
%         task_title = 'distraction';
%     end
%     plot_title =  sprintf('Correlation between the distance to hyperplane and reaction time in a %s task (N=%d)',task_title,numel(subjects));
% elseif task_distance==1 && task_RT==2
%     plot_title =  sprintf('Correlation between the distance to hyperplane from a categorization task and reaction time from a distraction task (N=%d)',numel(subjects));
% elseif task_distance==2 && task_RT==1
%     plot_title =  sprintf('Correlation between the distance to hyperplane from a distraction task and reaction time from a categorization task (N=%d)',numel(subjects));
% end
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
% legend_plot = {'Artificial scenes','Natural scenes','Both'}; 
% legend(legend_plot);
ylim([-0.3 0.3]);
yticks(-0.3:0.1:0.3);
xticks(0:20:200);
xticklabels(-200:100:800);  
xline(40,'--');
% legend_bool = 0;
% title_bool = 0;
% plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Spearman''s coefficient'); 

%% Save correlations and figures
dth_results.corr_both_categories = avg_corr_both;
dth_results.corr_artificial = avg_corr_art;
dth_results.corr_natural = avg_corr_nat;

save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
file_name = 'distance_rnn_rt_eeg';

for layer=1:num_layers
    save(fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_layer_%d.mat',subjects(1),subjects(end),file_name,model_name,layer)),'dth_results');
    saveas(gcf,fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_layer_%d.svg',subjects(1),subjects(end),file_name,model_name,layer))); 
    saveas(gcf,fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_layer_%d.fig',subjects(1),subjects(end),file_name,model_name,layer))); 
end

close(gcf);

end



