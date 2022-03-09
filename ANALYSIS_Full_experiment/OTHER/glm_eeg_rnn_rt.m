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
glm.avg_dth_shared = squeeze(mean(dth_shared,1));
glm.avg_dth_eeg_human_rt = squeeze(mean(dth_eeg_human_rt,1));
glm.avg_dth_rcnn_rt_human_rt = squeeze(mean(dth_rcnn_rt_human_rt,1));

%Timecourse plot
color_shared = [0.3 0 .5];
color_rcnn = [0 0 .7];
color_eeg = [.7 .1 0];
figure;
p_sh = plot(glm.avg_dth_shared,'LineWidth',2,'DisplayName','Shared variance','Color',color_shared);
hold on;
p_rCNN = plot(glm.avg_dth_rcnn_rt_human_rt,'LineWidth',2,'DisplayName','Unique variance rCNN RT','Color',color_rcnn);
p_eeg = plot(glm.avg_dth_eeg_human_rt,'LineWidth',2,'DisplayName','Unique variance EEG','Color',color_eeg);
set(gca,'FontName','Arial','FontSize',12);
xticks(0:20:200);
xticklabels(-200:100:800);        
onset_time = 40;
xline(onset_time,'--');
ylabel('r squared');
xlabel('Time (ms)');
% legend([p_sh,p_eeg,p_rCNN],'Location','best');
%Save plot
saveas(gcf,fullfile(results_avg_dir,sprintf('glm_rt_eeg_rcnn_rt_%s_scenes.fig',conditions))); 
saveas(gcf,fullfile(results_avg_dir,sprintf('glm_rt_eeg_rcnn_rt_%s_scenes.svg',conditions))); 
close(gcf); 

%Bar plot: at peak
peak_t = find(glm.avg_dth_shared==max(glm.avg_dth_shared));
glm.peak_shared = glm.avg_dth_shared(peak_t);
glm.peak_eeg_rt = glm.avg_dth_eeg_human_rt(peak_t);
glm.peak_rcnn_rt = glm.avg_dth_rcnn_rt_human_rt(peak_t); %always the same
peak_r2 = [glm.peak_shared glm.peak_eeg_rt glm.peak_rcnn_rt];
figure;
b = bar(peak_r2);
%Assign colours
b.FaceColor = 'flat';
b.CData(1,:) = color_shared;
b.CData(2,:) = color_rcnn;
b.CData(3,:) = color_eeg;
%labels
xlabels = {'Shared variance'; 'rCNN RT unique variance'; 'EEG unique variance'};
set(gca,'xticklabel',xlabels)
ylabel('rsquared')

%Save glm results and bar plot
saveas(gcf,fullfile(results_avg_dir,sprintf('glm_peak_rt_eeg_rcnn_rt_%s_scenes.fig',conditions)));
saveas(gcf,fullfile(results_avg_dir,sprintf('glm_peak_rt_eeg_rcnn_rt_%s_scenes.svg',conditions)));
close(gcf);

end


