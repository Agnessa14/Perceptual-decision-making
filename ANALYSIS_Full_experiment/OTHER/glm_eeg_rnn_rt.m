function glm_eeg_rnn_rt(conditions)
%GLM_EEG_RNN_RT Perform the GLM analysis to determine how much variance in
%human RT is explained by EEG data and rCNN RTs. 
%
%Input: conditions ('artificial','natural' or 'both')
%
%Output: results from the GLM analysis.
%
%Performs the GLM analysis on the selected scenes, using EEG data and rCNN RTs as regressors.
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
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

%% Run the GLM
dth_shared = NaN(num_subjects,num_timepoints);
dth_eeg_human_rt = NaN(num_subjects,num_timepoints);
dth_rcnn_rt_human_rt = NaN(num_subjects,num_timepoints);
for sub = 1:num_subjects
    num_scenes_sub = num_scenes;
    RTs_median = selected_human_RTs;
    selected_rCNN_RTs_sub = selected_rCNN_RTs;
    distances_eeg_sub = squeeze(selected_EEG_distances(sub,:,:)); 

    %remove excluded scene if needed
%     if sub == 6
%         keyboard;
%     end
    if any(isnan(distances_eeg_sub))
        excluded_scene = find(isnan(distances_eeg_sub(:,1)));
        distances_eeg_sub(excluded_scene,:) = [];
        RTs_median(excluded_scene) = [];
        selected_rCNN_RTs_sub(excluded_scene) = [];
        num_scenes_sub = num_scenes-1;
    end
        
    for t = 1:num_timepoints
        distances_eeg_sub_t = squeeze(distances_eeg_sub(:,t)); 
        
        [~,~,~,~,r2_eeg_human_rt] = regress(RTs_median, [ones(num_scenes_sub,1) distances_eeg_sub_t]); 
        [~,~,~,~,r2_rcnn_rt_human_rt] = regress(RTs_median, [ones(num_scenes_sub,1) selected_rCNN_RTs_sub']); 
        [~,~,~,~,r2_full] = regress(RTs_median, [ones(num_scenes_sub,1),distances_eeg_sub_t,selected_rCNN_RTs_sub']);

        dth_shared(sub,t) = r2_eeg_human_rt(1)+r2_rcnn_rt_human_rt(1)-r2_full(1); 
        dth_eeg_human_rt(sub,t) = r2_eeg_human_rt(1);
        dth_rcnn_rt_human_rt(sub,t) = r2_rcnn_rt_human_rt(1);
    end 
end


%% Plot
%Average over subjects
avg_dth_shared = squeeze(mean(dth_shared,1));
avg_dth_eeg_human_rt = squeeze(mean(dth_eeg_human_rt,1));
avg_dth_rcnn_rt_human_rt = squeeze(mean(dth_rcnn_rt_human_rt,1));

%Plot
p_sh = plot(avg_dth_shared,'LineWidth',2,'DisplayName','Shared variance');
hold on;
p_eeg = plot(avg_dth_eeg_human_rt,'LineWidth',2,'DisplayName','Unique variance EEG');
p_rCNN = plot(avg_dth_rcnn_rt_human_rt,'LineWidth',2,'DisplayName','Unique variance rCNN RT');
set(gca,'FontName','Arial','FontSize',12);
xticks(0:20:200);
xticklabels(-200:100:800);        
onset_time = 40;
xline(onset_time,'--');
ylabel('r squared');
xlabel('Time (ms)');
legend([p_sh,p_eeg,p_rCNN],'Location','best');

%Save glm results and plot
save()
saveas(gcf,fullfile(results_avg_dir,sprintf('glm_rt_eeg_rcnn_rt_%s_scenes.fig',conditions))); 
saveas(gcf,fullfile(results_avg_dir,sprintf('glm_rt_eeg_rcnn_rt_%s_scenes.svg',conditions))); 
close(gcf); 

end


