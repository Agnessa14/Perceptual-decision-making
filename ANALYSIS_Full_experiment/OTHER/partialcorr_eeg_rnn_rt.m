function partialcorr_eeg_rnn_rt(conditions,with_stats)
%PARTIALCORR_EEG_RNN_RT Calculate the partial correlations of EEG-RT and
%RNN-RT.
%
%Input: conditions ('artificial','natural' or 'both'), with/without stats
%(1/0)
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
num_scenes = numel(conds);

%% Calculate the partial correlations
r2_eeg_human_rt = NaN(num_subjects,num_timepoints);
r2_rcnn_rt_human_rt = NaN(num_subjects,num_timepoints);
r2_rcnn_eeg_full = NaN(num_subjects,num_timepoints);
common_all = NaN(num_subjects,num_timepoints);

for sub = 1:num_subjects
    num_scenes_sub = num_scenes;
    RTs_median_sub = selected_human_RTs;
    selected_rCNN_RTs_sub = selected_rCNN_RTs;
    distances_eeg_sub = squeeze(selected_EEG_distances(sub,:,:)); 

    %remove excluded scene if needed
    if any(isnan(distances_eeg_sub))
        excluded_scene = find(isnan(distances_eeg_sub(:,1)));
        distances_eeg_sub(excluded_scene,:) = [];
        RTs_median_sub(excluded_scene) = [];
        selected_rCNN_RTs_sub(excluded_scene) = [];
        num_scenes_sub = num_scenes-1;
    end
        
    for t = 1:num_timepoints
        distances_eeg_sub_t = squeeze(distances_eeg_sub(:,t)); 

        r_eeg = correlate([RTs_median_sub distances_eeg_sub_t*-1],'type','spearman','method','partialcorr');
        r_rcnn = correlate([RTs_median_sub selected_rCNN_RTs_sub'],'type','spearman','method','partialcorr');
        r_full = correlate([RTs_median_sub distances_eeg_sub_t*-1 selected_rCNN_RTs_sub'],'type','spearman','method','partialcorr');
        
        r2_eeg_human_rt(sub,t) = r_eeg(2,1).^2;
        r2_rcnn_rt_human_rt(sub,t) = r_rcnn(2,1).^2;
        r2_rcnn_eeg_full(sub,t) = r_full(2,1).^2;
        
        common_all(sub,t) = r_rcnn(2,1).^2 + r_eeg(2,1).^2 - r_full(2,1).^2;

    end 
end


%% Plot
%Average over subjects
partial_corr.eeg = squeeze(mean(r2_eeg_human_rt,1));
partial_corr.rcnn = squeeze(mean(r2_rcnn_rt_human_rt,1));
partial_corr.full = squeeze(mean(r2_rcnn_eeg_full,1));
partial_corr.shared = squeeze(mean(common_all,1));

%Color parameters
color_shared = [.7 0 .7];
color_rcnn = [0 .3 .2];
color_eeg = [.5 .1 0];
color_full = [0.8 0.8 0.8];

%plot the full r^2 first as a shaded area
figure;
p_full = area(1:num_timepoints,sqrt(partial_corr.full)); 
hold on;
        
p_eeg = plot(sqrt(partial_corr.eeg),'LineWidth',2,'DisplayName','Partial correlation EEG','Color',color_eeg);
hold on;
p_rcnn = plot(sqrt(partial_corr.rcnn),'LineWidth',2,'DisplayName','Partial correlation rCNN RT','Color',color_rcnn);
p_shared = plot(sqrt(partial_corr.shared),'LineWidth',2,'DisplayName','Shared rCNN and EEG','Color',color_shared);
%Plot full as a shaded area

%Plot parameters
p_full.FaceColor = color_full;
p_full.EdgeColor = color_full;
p_full.DisplayName = 'Full';
set(gca,'FontName','Arial','FontSize',12);
xticks(0:20:200);
xticklabels(-200:100:800);        
onset_time = 40;
xline(onset_time,'--');
ylabel('R^{2}');
xlabel('Time (ms)');
legend([p_eeg,p_rcnn,p_shared,p_full],'Location','best');

%Save results & plot 
save(fullfile(results_avg_dir,sprintf('partial_corr_eeg_rcnn_rt_%s_scenes',conditions)),'partial_corr');
saveas(gcf,fullfile(results_avg_dir,sprintf('partial_corr_eeg_rcnn_rt_%s_scenes.fig',conditions))); 
saveas(gcf,fullfile(results_avg_dir,sprintf('partial_corr_eeg_rcnn_rt_%s_scenes.svg',conditions))); 
close(gcf); 

end


