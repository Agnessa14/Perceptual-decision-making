function for_stats = all_subjects_for_stats(subjects,task,analysis)
%ALL_SUBJECTS_FOR_STATS Gather the subject x time matrices needed to perform cluster-based permutation tests.
%
%Input: subject IDs, task (1=categorization,2=distraction), analysis
%('object_decoding', 'category_decoding', 'time_object_decoding', 'time_category_decoding' or 'rsa_time_object')
%
%Output: SxP matrix containing decoding accuracies, where S is the number
%of subjects and P is the number of timepoints.
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
if task < 3
    task_name = get_task_name(task);
elseif task == 3
    if strcmp(analysis,'object_decoding') || strcmp(analysis,'time_object_decoding')
        task_name = 'crosstask';
    elseif strcmp(analysis,'category_decoding') || strcmp(analysis,'time_category_decoding')
        task_name = 'cross_task';
    end
end

%% Preallocate
numConditions = 60;
if strcmp(analysis,'object_decoding')
    numTimepoints = 200;
    decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints);
    filename = sprintf('svm_decoding_accuracy_%s.mat',task_name);
elseif strcmp(analysis,'category_decoding')
    numTimepoints = 200;
    decoding_accuracies_all_subjects = NaN(max(subjects),numTimepoints);
    filename = sprintf('svm_artificial_vs_natural_decoding_accuracy_%s.mat',task_name);
elseif strcmp(analysis,'time_object_decoding')
    numTimepoints = 50;
    decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints,numTimepoints);
    filename = sprintf('time_generalized_svm_object_decoding_%s.mat',task_name);
elseif strcmp(analysis,'time_category_decoding')
    numTimepoints = 50;
    decoding_accuracies_all_subjects = NaN(max(subjects),numTimepoints,numTimepoints);
    filename = sprintf('3PT_time_gen_svm_artificial_vs_natural_decoding_accuracy_%s.mat',task_name);
elseif strcmp(analysis,'rsa_time_object')
    numTimepoints = 200;
    decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints);
    filename = sprintf('rsa_pearson_%s.mat',task_name);
end 

dimensionality = size(size(decoding_accuracies_all_subjects),2);
all_dimensions = repmat({':'},1,dimensionality);

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,filename));
    if contains(analysis,'time')
        decoding_accuracies_all_subjects(subject,all_dimensions{:}) = timeg_decodingAccuracy_avg;  
    else
        decoding_accuracies_all_subjects(subject,all_dimensions{:}) = decodingAccuracy_avg;  
    end
end   

%% Average over conditions and remove any NaN subjects
if contains(analysis,'category_decoding')
    for_stats = decoding_accuracies_all_subjects;
else
    for_stats = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects,2),3));
end

size_for_stats = size(for_stats);
for_stats = for_stats(~isnan(for_stats(:,1)),:);
size_final = [size(for_stats,1) size_for_stats(2:end)];
for_stats = reshape(for_stats,size_final);

%% Subtract 50 to get the difference with chance (needed for stats)
if ~contains(analysis,'rsa') %only for analyses with SVM
    for_stats = for_stats-50;
end

%% Save
save(fullfile(results_avg_dir,sprintf('for_stats_subjects_%d_%d_%s_task_%s',subjects(1),subjects(end),task_name,analysis)),'for_stats');
end
