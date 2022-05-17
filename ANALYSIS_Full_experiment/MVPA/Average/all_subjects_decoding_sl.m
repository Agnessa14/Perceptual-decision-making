function all_subjects_decoding_sl(subjects,analysis,task)
%ALL_SUBJECTS_DECODING_SL (categorization and distraction). Collects decoding the
%data from subjects to put it in one matrix.
%
%Input: subject IDs, analysis ('svm_object_decoding', 'pearson_object_decoding' or 'category_decoding'),
%task (1 or 2)
%Output: matrix containing decoding accuracies for multiple subjects,
%timepoints and channels
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate
numConditions = 60;
numTimepoints = 10;
numTimepointsDTH = 200;
numChannels = 63;
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints,numChannels);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects = NaN(max(subjects),numTimepoints,numChannels);
    distances_all_subjects = NaN(max(subjects),numConditions,numTimepointsDTH,numChannels);
end
task_name = get_task_name(task);

%% Loop: collect results from all subjects
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'svm_object_decoding')
        filename = sprintf('svm_searchlight_object_decoding_accuracy_%s.mat',task_name);
        all_dimensions = repmat({':'},1,4); %conditions x conditions x timepoints x channels
    elseif strcmp(analysis,'pearson_object_decoding')
        filename = sprintf('rdm_pearson_searchlight_%s.mat',task_name);
        all_dimensions = repmat({':'},1,4); %conditions x conditions x timepoints x channels        
    elseif strcmp(analysis,'category_decoding')
        filename = sprintf('cross_validated_dth_pseudotrials_svm_decodingAccuracy_searchlight_%s.mat',task_name);
        all_dimensions = repmat({':'},1,3); % timepoints x channels 
    end
    filename_full = fullfile(subject_results_dir,filename);
    if ~exist(filename_full,'file')
        disp(subject);
    else
        if strcmp(analysis,'svm_object_decoding') || strcmp(analysis,'category_decoding')
            load(filename_full,'decodingAccuracy_avg');
            decoding_accuracies_all_subjects(subject,all_dimensions{:}) = decodingAccuracy_avg;
            filename_dist = sprintf('cross_validated_dth_pseudotrials_svm_decisionValues_searchlight_all_timepoints_%s.mat',task_name);
            load(fullfile(subject_results_dir,filename_dist),'decisionValues_Avg');
            distances_all_subjects(subject,all_dimensions{:}) = decisionValues_Avg;
        elseif strcmp(analysis,'pearson_object_decoding')
            load(filename_full,'rdm_avg');
            decoding_accuracies_all_subjects(subject,all_dimensions{:}) = rdm_avg;
        end
    end
end

decoding_accuracies = decoding_accuracies_all_subjects(subjects,all_dimensions{:});

%% Save the matrix
if strcmp(analysis,'category_decoding')
    distances = distances_all_subjects(subjects,all_dimensions{:});
    save(fullfile(results_avg_dir,sprintf('distances_all_subjects_all_timepoints_searchlight_%s.mat',task_name)),'distances'); 
    save(fullfile(results_avg_dir,sprintf('%s_all_subjects_all_timepoints_searchlight_%s.mat',analysis,task_name)),'decoding_accuracies'); 
else
    save(fullfile(results_avg_dir,sprintf('%s_all_subjects_searchlight_%s.mat',analysis,task_name)),'decoding_accuracies'); 
end