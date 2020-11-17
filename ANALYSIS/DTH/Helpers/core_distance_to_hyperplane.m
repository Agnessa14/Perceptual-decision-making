function correlation_dth_rt = core_distance_to_hyperplane(subjects,subset,category,numTimepoints)
%DISTANCE_TO_HYPERPLANE Correlate the distance to hyperplane and RTs over
%time.
%Input: 
%-subject IDs (e.g., 1:13)
%-subset of data (e.g., 10%)
%-category of object (e.g., artificial)
%-number of timepoints to consider(e.g., 200)
%
%
%% Define the condition indices
if strcmp(category, 'artificial')
    conditions = 1:30;
elseif strcmp(category, 'natural')
    conditions = 31:60;
end

%% Preallocate
distances = NaN(numel(subjects),numel(conditions),numTimepoints);
RTs = NaN(numel(subjects),numel(conditions));

%% Get the RTs and distances of all subjects for each condition
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('subset_%s_decisionValues.mat',num2str(subset))));
    distances(subject,:,:) = decisionValues_Avg;   
    load(fullfile(results_dir,subname,sprintf('subset_%s_RTs_correct_trials.mat',num2str(subset))));
    RTs(subject,:) = normalize(RT_per_condition);
end

medianRT = nanmedian(RTs,1);
mean_distances = squeeze(nanmean(distances,1)); %avg over subjects

%% Correlate distances and RTs
correlation_dth_rt = arrayfun(@ (x) corr(mean_distances(conditions,x),...
        medianRT(conditions)','type','Spearman'), 1:numTimepoints);
end