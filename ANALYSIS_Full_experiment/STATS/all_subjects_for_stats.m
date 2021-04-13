function for_stats = all_subjects_for_stats(subjects,task,analysis)
%ALL_SUBJECTS_FOR_STATS Gather the subject x time matrices needed to perform cluster-based permutation tests.
%
%Input: subject IDs, task (1=categorization,2=distraction), analysis
%('object_decoding' or 'category_decoding')
%
%Output: SxP matrix containing decoding accuracies, where S is the number
%of subjects and P is the number of timepoints.

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
task_name = get_task_name(task);

%% Preallocate
numConditions = 60;
numTimepoints = 200;

%% Set up the figure for plotting
sorted_subjects = sort(subjects); %order by ID
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
    filename = sprintf('svm_decoding_accuracy_%s.mat',task_name);
    all_dimensions = repmat({':'},1,3);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects = NaN(sorted_subjects(end),numTimepoints);
    filename = sprintf('svm_artificial_vs_natural_decoding_accuracy_%s.mat',task_name);
    all_dimensions = ':';
end 

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,filename));
    decoding_accuracies_all_subjects(subject,all_dimensions) = decodingAccuracy_avg;  
end   

%% Average over Remove any NaN 
if strcmp(analysis,'object_decoding')
    for_stats = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects,2),3));
    for_stats = for_stats(~isnan(for_stats(:,1)),:); %for non-included subjects: if timepoint 1 is NaN, the subject is excluded
elseif strcmp(analysis,'category_decoding')
    for_stats = decoding_accuracies_all_subjects;
    for_stats = for_stats(~isnan(for_stats(:,1)),:); %for non-included subjects: if timepoint 1 is NaN, the subject is excluded
end
    %% Subtract 50 to get the difference with chance (needed for stats)
for_stats = for_stats-50;

%% Save
save(fullfile(results_avg_dir,sprintf('for_stats_%d_subjects_%s_task_%s',numel(subjects),task_name,analysis)),'for_stats');

end
