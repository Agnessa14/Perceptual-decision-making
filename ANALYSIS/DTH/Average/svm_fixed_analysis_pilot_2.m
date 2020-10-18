function svm_fixed_analysis_pilot_2(subjects) %distance art, distance nat, RT
%SVM_FIXED_ANALYSIS_PILOT_2 Performs the distance-to-hyperplane analysis using
%the svm classifier 60 scenes on normalized distances.
%
%Input: subjects' ID (e.g., 1:13)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS';
addpath(genpath(results_dir));

%% Get the distances from all subjects
numTimepoints = 200;
numConditions = 60;
distances = NaN(numel(subjects),numConditions,numTimepoints);
RTs = NaN(numel(subjects),numConditions);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,'pseudotrials_decisionValues.mat'));
    distances(subject,:,:) = decisionValues_Avg;   
    load(fullfile(results_dir,subname,'RTs_correct_trials.mat'));
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
last_category1_sample = 30;
first_category2_sample = 31;

for t = 1:numTimepoints
    correlation_dth_rt_both(t) = corr(mean_distances(:,t),medianRT','type','Spearman');
    correlation_dth_rt_art(t) = corr(mean_distances(1:last_category1_sample,t),...
        medianRT(1:last_category1_sample)','type','Spearman');
    correlation_dth_rt_nat(t) = corr(mean_distances(first_category2_sample:end,t),...
        medianRT(first_category2_sample:end)','type','Spearman');
end
correlation_dth_rt_avg = mean([correlation_dth_rt_art;correlation_dth_rt_nat],1);

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen

%Artificial scenes
plot(correlation_dth_rt_art,'LineWidth',2);
hold on;

%Natural scenes
plot(correlation_dth_rt_nat,'LineWidth',2);
hold on;

%Average of artificial and natural
plot(correlation_dth_rt_avg,'LineWidth',2);
hold on;

% axis([0,200,-0.7,0.6]);
%Legend
legend_art_nat = {'Artificial scenes','Natural scenes','Average of artificial and natural scenes', ...
    'All scenes','Stimulus onset'};

%Both scenes in one plot
plotting_dth_one(correlation_dth_rt_both,...
    sprintf('Correlation between the distance to hyperplane and reaction time for 60 scenes in a categorization task (N=%d)', numel(subjects)),legend_art_nat);

%% Save
%correlations
save_path = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';
save(fullfile(save_path,sprintf('pt_SVM_DTH_rt_correlation_both_categories_%d_subjects',numel(subjects))),'correlation_dth_rt_both');
save(fullfile(save_path,sprintf('pt_SVM_DTH_rt_correlation_artificial_%d_subjects',numel(subjects))),'correlation_dth_rt_art');
save(fullfile(save_path,sprintf('pt_SVM_DTH_rt_correlation_natural_%d_subjects',numel(subjects))),'correlation_dth_rt_nat');
save(fullfile(save_path,sprintf('pt_SVM_DTH_rt_correlation_both_categories_avg_%d_subjects',numel(subjects))),'correlation_dth_rt_avg');

%figures
saveas(gcf,fullfile(save_path,sprintf('pt_SVM_DTH_artificial_natural_%d_subjects',numel(subjects))));
saveas(gcf,fullfile(save_path,sprintf('pt_SVM_DTH_artificial_natural_%d_subjects.svg',numel(subjects))));

end