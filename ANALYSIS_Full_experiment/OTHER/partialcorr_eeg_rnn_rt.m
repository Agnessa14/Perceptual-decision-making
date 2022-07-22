function partialcorr_eeg_rnn_rt(conditions,with_stats,with_rcnn)
%PARTIALCORR_EEG_RNN_RT Calculate the partial correlations of EEG-RT and
%RNN-RT.
%
%Input: conditions ('artificial','natural' or 'both'), with/without stats
%(1/0), with or without average rCNN unique variance (1/0)
%
%Output: results from the GLM analysis.
%
%Calculates and plots the partial correlations (unique and shared) associated
%with the relationships between EEG and human RTs, and RNN RTs and human RTs.
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
addpath(genpath('/home/agnek95/helper_functions')); %correlate from M. Hebart
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';

%% Load the EEG data (distances to the hyperplane for all subjects), rCNN RTs and human RTs 
load(fullfile(results_avg_dir,'distances_all_subjects.mat'),'distances_all_subjects');
load('/scratch/agnek95/PDM/DATA/RNN_RTs/RNN_RTs_entropy_threshold_0.02.mat','data');
rCNN_RTs = data;
load(fullfile(results_avg_dir,'RT_all_subjects_5_35_categorization.mat'),'RTs');
human_RTs = nanmedian(RTs,1);

%Select the appropriate scenes
if strcmp(conditions,'artificial')
    conds = 1:30;
elseif strcmp(conditions,'natural')
    conds = 31:60;
elseif strcmp(conditions,'both')
    conds = 1:60;
end
selected_EEG_distances = distances_all_subjects(:,conds,:);
selected_human_RTs = human_RTs(conds)';
selected_rCNN_RTs = rCNN_RTs(conds);

%Define some variables
num_subjects = size(distances_all_subjects,1);
num_timepoints = size(distances_all_subjects,3);
% num_scenes = numel(conds);

%% Calculate the partial correlations
size_data = [num_subjects,num_timepoints];
r_eeg_unique = NaN(size_data);
r_rcnn_unique = NaN(size_data);
r_eeg_full = NaN(size_data);
r_rcnn_full = NaN(size_data);
shared_rt_eeg = NaN(size_data);
shared_rt_rcnn = NaN(size_data);

for sub = 1:num_subjects
%     num_scenes_sub = num_scenes;
    RTs_median_sub = selected_human_RTs;
    selected_rCNN_RTs_sub = selected_rCNN_RTs;
    distances_eeg_sub = squeeze(selected_EEG_distances(sub,:,:)); 

    %remove excluded scene if needed
    if any(isnan(distances_eeg_sub))
        excluded_scene = find(isnan(distances_eeg_sub(:,1)));
        distances_eeg_sub(excluded_scene,:) = [];
        RTs_median_sub(excluded_scene) = [];
        selected_rCNN_RTs_sub(excluded_scene) = [];
%         num_scenes_sub = num_scenes-1;
    end
        
    for t = 1:num_timepoints
        distances_eeg_sub_t = squeeze(distances_eeg_sub(:,t)); 
        
        %unique = semipartial
        r_eeg_u = correlate([RTs_median_sub, distances_eeg_sub_t*-1, selected_rCNN_RTs_sub'],'type','spearman','method','partialcorr');
        r_eeg_unique(sub,t) = r_eeg_u(2,1);
        r_rcnn_u = correlate([RTs_median_sub, selected_rCNN_RTs_sub', distances_eeg_sub_t*-1],'type','spearman','method','partialcorr');
        r_rcnn_unique(sub,t) = r_rcnn_u(2,1);
        
        %full 
        r_eeg_full(sub,t) = corr(RTs_median_sub,distances_eeg_sub_t*-1,'type','spearman');
        r_rcnn_full(sub,t) = corr(RTs_median_sub,selected_rCNN_RTs_sub','type','spearman');
        
        %shared = full (normal) -semipartial
        shared_rt_eeg(sub,t) = r_eeg_full(sub,t).^2-r_eeg_unique(sub,t).^2;
        shared_rt_rcnn(sub,t) = r_rcnn_full(sub,t).^2-r_rcnn_unique(sub,t).^2;   %shared_rt_eeg=shared_rt_rcnn     
    end 
end

%Plot colors
color_shared = [0.5 0.7 1];
color_rcnn =  [1 0.5 0.8];
color_eeg =   [0.8 0.9 0.3]; %[0.466 0.674 0.188];

%Average over subjects
partial_corr.eeg = squeeze(mean(r_eeg_unique,1));
partial_corr.rcnn = squeeze(mean(r_rcnn_unique,1));
partial_corr.shared = squeeze(mean(sign(shared_rt_eeg).*sqrt(abs(shared_rt_eeg)),1)); %put on quadratic scale

%% Plot
legend_names = cell(3,1);
for c = [1,3] % don't plot rcnn rt unjique variance - looks weird %EEG, rCNN and shared 
    switch c
        case 1
            for_stats_data = r_eeg_unique;
            plot_location = -0.15;
            color_data = color_eeg;
            data = partial_corr.eeg;
            legend_entry = 'EEG unique variance';
            variable_name = 'eeg_unique';
        case 2
            for_stats_data = r_rcnn_unique;
            plot_location = -0.16;
            color_data = color_rcnn;
            data = partial_corr.rcnn;
            legend_entry = 'rCNN RTs unique variance';
            variable_name = 'rcnn_unique';
        case 3
            for_stats_data = shared_rt_eeg;
            plot_location = -0.17;
            color_data = color_shared;
            data = partial_corr.shared;
            legend_entry = 'Shared variance';
            variable_name = 'shared';
    end
    
    %Stats
    if with_stats
        %Load stats if they exist
        filename = fullfile(results_avg_dir,sprintf('commonality_analysis_stats_fdr_subjects_%s_%s.mat',conditions,variable_name));
        if exist(filename,'file')
            load(filename,'stats_comm_a');
        else
            stats_comm_a.num_perms = 10000;
            stats_comm_a.tail = 'right';
            stats_comm_a.qvalue = 0.05;
            [stats_comm_a.significant_timepoints,stats_comm_a.pvalues,...
                stats_comm_a.crit_p, stats_comm_a.adjusted_pvalues]...
                = fdr_permutation_cluster_1sample_alld(for_stats_data,...
                stats_comm_a.num_perms,stats_comm_a.tail,stats_comm_a.qvalue);
            save(filename,'stats_comm_a');
        end

        %1) Plot significant timepoints
        st = (stats_comm_a.significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color_data); 
        hold on;
    end
    
    %2) Calculate and plot error bars
    stdDM = std(for_stats_data); 
    err = stdDM/sqrt(size(for_stats_data,1)); %standard deviation/sqrt of num subjects  

    %plot as a shaded area
    top_curve = data + err;
    bottom_curve = data - err;
    x2 = [1:num_timepoints, fliplr(1:num_timepoints)];
    shaded_area = [top_curve, fliplr(bottom_curve)];
    fill(x2, shaded_area, color_data,'FaceAlpha',0.5,'LineStyle','none');
    hold on;

    %3) Plot the curves
    legend_names{c} = plot(data,'LineWidth',2,'Color',color_data,'DisplayName',legend_entry);
end   

%% Plot the average rCNN RT-human RT unique variance as a dot
if with_rcnn
    scatter(200,mean(partial_corr.rcnn),100,color_rcnn,'filled');
end
%% Plot parameters
set(gca,'FontName','Arial','FontSize',12);
xticks(0:20:200);
xticklabels(-200:100:800);   
ylim([-0.2,0.5]);
onset_time = 40;
xline(onset_time,'--');
ylabel('Commonality coefficient');% ylabel('R^{2}');
xlabel('Time (ms)');
% legend([legend_names{1},legend_names{2},legend_names{3}],'Location','best');

%Save results & plot 
filename_save = sprintf('partial_corr_eeg_rcnn_rt_%s_scenes_modif',conditions);
if with_rcnn
    filename_fig = sprintf('with_rcnn_%s',filename_save);
end
save(fullfile(results_avg_dir,filename_save),'partial_corr');
saveas(gcf,fullfile(results_avg_dir,sprintf('%s.fig',filename_fig))); 
saveas(gcf,fullfile(results_avg_dir,sprintf('%s.svg',filename_fig))); 
close(gcf); 

end


