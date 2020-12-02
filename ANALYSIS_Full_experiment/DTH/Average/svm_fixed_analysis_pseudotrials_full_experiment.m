function svm_fixed_analysis_pseudotrials_full_experiment(subjects,task) %distance art, distance nat, RT
%SVM_FIXED_ANALYSIS_PSEUDOTRIALS_FULL_EXPERTIMENT Performs the distance-to-hyperplane analysis using
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
    distances(subject,:,:) = decisionValues_Avg;   
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name)));
    RTs(subject,:) = normalize(RT_per_condition);
end

%% Get the median RTs of all subjects for each condition
%Normalize and get median
medianRT = nanmedian(RTs,1);

%% Correlate DTH and RT
mean_distances = squeeze(nanmean(distances,1)); %avg over subjects
correlation_dth_rt_both = NaN(1,numTimepoints);
correlation_dth_rt_art = NaN(1,numTimepoints);
correlation_dth_rt_nat = NaN(1,numTimepoints);

for t = 1:numTimepoints
    correlation_dth_rt_both(t) = corr(mean_distances(:,t),medianRT','type','Spearman');
    correlation_dth_rt_art(t) = corr(mean_distances(1:numConditions/2,t),...
        medianRT(1:last_category1_sample)','type','Spearman');
    correlation_dth_rt_nat(t) = corr(mean_distances((numConditions/2)+1:end,t),...
        medianRT(first_category2_sample:end)','type','Spearman');
end
correlation_dth_rt_avg = mean([correlation_dth_rt_art;correlation_dth_rt_nat],1);

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
plot(correlation_dth_rt_art,'LineWidth',2);
hold on;
plot(correlation_dth_rt_nat,'LineWidth',2);
hold on;
plot(correlation_dth_rt_avg,'LineWidth',2);
hold on;

%Plotting parameters
title =  sprintf('Correlation between the distance to hyperplane and reaction time for 60 scenes in a %s task (N=%d)', numel(subjects),task_name);
legend_plot = {'Artificial scenes','Natural scenes','Average of artificial and natural scenes', ...
    'All scenes','Stimulus onset'};
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