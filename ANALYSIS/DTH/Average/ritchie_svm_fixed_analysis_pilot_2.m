function ritchie_svm_fixed_analysis_pilot_2(subjects,conditions) %conditions 1:12 for animate or 13:24 for inanimate
%RITCHIE_SVM_FIXED_ANALYSIS_PILOT_2 Performs the distance-to-hyperplane analysis using
%the svm classifier on 24 objects (on Ritchie 2015 data).
%
%Input: subjects' ID (e.g., 1:13)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (24 objects), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS';
addpath(genpath(results_dir));

%% Get the distances from all subjects
numTimepoints = 36;
distances = NaN(numel(subjects),numel(conditions),numTimepoints);
RTs = NaN(numel(subjects),numel(conditions));

for subject = 1:numel(subjects)
    s = subjects(subject);
    subname = get_subject_name(s);
    load(fullfile(results_dir,subname,'witherrortrials_svm_ritchie_decisionValues.mat'));
    distances(subject,:,:) = decisionValues_Avg(conditions,:);   
    load(fullfile(results_dir,subname,'RT_active_mean.mat'));
    RTs(subject,:) = normalize(meanRT(conditions));
end

%% Get the median RTs of all subjects for each condition
medianRT = median(RTs,1);

%% Correlate DTH and RT
avg_distance = squeeze(mean(distances,1)); %avg over subjects

correlation_type = 'Spearman';
timepoints = 1:numTimepoints;
correlation_dth_RT = arrayfun(@ (x) corr(medianRT',avg_distance(:,x),'type',correlation_type),timepoints);

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
plot(correlation_dth_RT,'LineWidth',2);
hold on;
analyze_ritchie_data(subjects,conditions);
title(sprintf('Correlation between reaction time and distance to hyperplane in 12 animate objects (N=%d)',numel(subjects)));
legend({'My script','Ground truth'},'FontSize',12)
xlabel('Timepoint')
ylabel('Spearman''s coefficient')

%% Save figures
save_path = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';
saveas(gcf,fullfile(save_path,sprintf('good_ritchie_myscript_vs_ground_animate_%d_subjects',numel(subjects))));
saveas(gcf,fullfile(save_path,sprintf('good_ritchie_myscript_vs_ground_animate_%d_subjects.svg',numel(subjects))));

end