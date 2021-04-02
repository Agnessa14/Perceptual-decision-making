function dth_all_distances_median_RT(subjects,task,with_stats) %distance art, distance nat, RT
%DTH_ALL_DISTANCES_MEDIAN_RT Performs the distance-to-hyperplane analysis using
%the distances-to-hyperplane for each subject and the median RT across subjects.
%
%Input: subjects' ID (e.g., 1:13), task (1=categoization, 2=distraction), add stats (1 with, 0 without)
%
%Correlates the decision values of each subject with reaction times averaged over
%participants of each condition (60 scenes), averaged, at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
addpath(genpath(results_dir));
task_name = get_task_name(task);

%% Get the distances from all subjects
numTimepoints = 200;
numConditions = 60;
sorted_subjects = sort(subjects); %order by ID
distances = NaN(sorted_subjects(end),numConditions,numTimepoints);
RTs = NaN(numel(subjects),numConditions);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_name)));
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name)));
    distances(subject,:,:) = decisionValues_Avg;   
    RTs(subject,:) = normalize(RT_per_condition);
end

%% Get the median RTs of all subjects for each condition 
medianRT = nanmedian(RTs,1);

%% Correlate each subject's distances with the median RT
t = 1:numTimepoints;
size_corr = [numel(subjects),numTimepoints];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);
correlation_avg = NaN(size_corr);

for subject = subjects
    correlation_art(subject,:) = arrayfun(@(x) corr(squeeze(distances(subject,1:numConditions/2,x)),medianRT,'type','Spearman'),t,'UniformOutput',false);
    correlation_nat(subject,:) = arrayfun(@(x) corr(squeeze(distances(subject,(numConditions/2)+1:end,x)),medianRT,'type','Spearman'),t,'UniformOutput',false);
    correlation_both(subject,:) = arrayfun(@(x) corr(squeeze(distances(subject,:,x)),medianRT,'type','Spearman'),t,'UniformOutput',false);
    correlation_avg(subject,:) = mean([squeeze(correlation_art(subject,:));squeeze(correlation_nat(subject,:))],1);
end

%% Average over participants
avg_corr_art = squeeze(nanmean(correlation_art,1));
avg_corr_nat = squeeze(nanmean(correlation_nat,1));
avg_corr_both = squeeze(nanmean(correlation_both,1));
avg_corr_avg = squeeze(nanmean(correlation_avg,1));


% correlation_dth_rt_art  = arrayfun(@(x) corr(mean_distances(1:numConditions/2,x),medianRT(1:numConditions/2)','type','Spearman'),t);
% correlation_dth_rt_nat  = arrayfun(@(x) corr(mean_distances((numConditions/2)+1:end,x),medianRT((numConditions/2)+1:end)','type','Spearman'),t);
% correlation_dth_rt_avg = mean([correlation_dth_rt_art;correlation_dth_rt_nat],1);






