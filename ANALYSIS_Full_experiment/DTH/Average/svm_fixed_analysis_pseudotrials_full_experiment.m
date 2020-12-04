function svm_fixed_analysis_pseudotrials_full_experiment(subjects,task) %distance art, distance nat, RT
%SVM_FIXED_ANALYSIS_PSEUDOTRIALS_FULL_EXPERIMENT Performs the distance-to-hyperplane analysis using
%the svm classifier 60 scenes on normalized distances.
%
%Input: subjects' ID (e.g., 1:13), task (1=categoization, 2=distraction);
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
addpath(genpath(results_dir));
task_name = get_task_name(task);

%% Get the distances from all subjects
numTimepoints = 200;
numConditions = 60;
distances = NaN(numel(subjects),numConditions,numTimepoints);
RTs = NaN(numel(subjects),numConditions);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_name)));
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name)));
    if ismember(subject,1:4) && task==2
        decisionValues_Avg = cat(1,decisionValues_Avg,NaN(2,numTimepoints));
        RT_per_condition = cat(1,RT_per_condition,NaN(2,1));
    end
    distances(subject,:,:) = decisionValues_Avg;   
    RTs(subject,:) = normalize(RT_per_condition);
end

%% Get the median RTs and mean distances of all subjects for each condition 
included_conditions = find(~isnan(RTs(subjects(1),:))); %can be any subject
numConditionsIncluded = numel(included_conditions);
medianRT = nanmedian(RTs(subjects,included_conditions),1);
mean_distances = squeeze(nanmean(distances(subjects,included_conditions,:),1)); %avg over subjects

%% Correlate DTH and RT
t = 1:numTimepoints;
correlation_dth_rt_both = arrayfun(@(x) corr(mean_distances(:,x),medianRT','type','Spearman'),t);
correlation_dth_rt_art  = arrayfun(@(x) corr(mean_distances(1:numConditionsIncluded/2,x),medianRT(1:numConditionsIncluded/2)','type','Spearman'),t);
correlation_dth_rt_nat  = arrayfun(@(x) corr(mean_distances((numConditionsIncluded/2)+1:end,x),medianRT((numConditionsIncluded/2)+1:end)','type','Spearman'),t);
correlation_dth_rt_avg = mean([correlation_dth_rt_art;correlation_dth_rt_nat],1);

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
plot(correlation_dth_rt_both,'LineWidth',2);
hold on;
plot(correlation_dth_rt_art,'LineWidth',2);
hold on;
plot(correlation_dth_rt_nat,'LineWidth',2);
hold on;
plot(correlation_dth_rt_avg,'LineWidth',2);
hold on;

%Plotting parameters
if task==1
    task_title = 'scene categorization';
elseif task==2
    task_title = 'distraction';
end
title =  sprintf('Correlation between the distance to hyperplane and reaction time in a %s task (N=%d)',task_title,numel(subjects));
legend_plot = {'All scenes','Artificial scenes','Natural scenes',...
    'Average of artificial and natural scenes','Stimulus onset'};
plotting_parameters(title,legend_plot,40)

%% Save
%correlations
save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
save(fullfile(save_path,sprintf('pseudotrials_SVM_DTH_rt_both_categories_%d_subjects_%s.mat',numel(subjects),task_name)),'correlation_dth_rt_both');
save(fullfile(save_path,sprintf('pseudotrials_SVM_DTH_rt_artificial_%d_subjects_%s.mat',numel(subjects),task_name)),'correlation_dth_rt_art');
save(fullfile(save_path,sprintf('pseudotrials_SVM_DTH_rt_natural_%d_subjects_%s.mat',numel(subjects),task_name)),'correlation_dth_rt_nat');
save(fullfile(save_path,sprintf('pseudotrials_SVM_DTH_rt_avg_%d_subjects_%s.mat',numel(subjects),task_name)),'correlation_dth_rt_avg');

%figures
saveas(gcf,fullfile(save_path,sprintf('pseudotrials_SVM_DTH_%d_subjects_%s',numel(subjects),task_name))); 
saveas(gcf,fullfile(save_path,sprintf('pseudotrials_SVM_DTH_%d_subjects_%s.svg',numel(subjects),task_name)));

end