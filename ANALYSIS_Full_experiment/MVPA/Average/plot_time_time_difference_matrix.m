function plot_time_time_difference_matrix(subjects,analysis)
%PLOT_TIME_TIME_DIFFERENCE_MATRIX Plot the difference between categorization & distraction within-task 
%time-generalization, averaged over all participants.
%
%Input: subject IDs (e.g., 1:13), analysis('object_decoding' or 'category_decoding')
%
%Output: 2D heatmap of difference between decoding accuracies. 
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate & collect results from all subjects
numTimepoints = 50;
numConditions = 60;
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects_categorization = NaN(max(subjects),numConditions,numConditions,numTimepoints,numTimepoints);
    decoding_accuracies_all_subjects_distraction = NaN(max(subjects),numConditions,numConditions,numTimepoints,numTimepoints);
    all_dimensions = repmat({':'},1,4);  
    analysis_name = analysis(1:6);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects_categorization = NaN(max(subjects),numTimepoints,numTimepoints);
    decoding_accuracies_all_subjects_distraction = NaN(max(subjects),numTimepoints,numTimepoints);
    all_dimensions = repmat({':'},1,2);
    analysis_name = analysis(1:8);
end

%Loop over subjects
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        filename = 'time_generalized_svm_object_decoding';
    elseif strcmp(analysis,'category_decoding')
        filename = 'time_generalized_svm_object_decoding';
    end
    load(fullfile(subject_results_dir,sprintf('%s_categorization.mat',filename)));
    decoding_accuracies_all_subjects_categorization(subject,all_dimensions{:}) = timeg_decodingAccuracy_avg;
    load(fullfile(subject_results_dir,sprintf('%s_fixation.mat',filename)));
    decoding_accuracies_all_subjects_distraction(subject,all_dimensions{:}) = timeg_decodingAccuracy_avg; 
end   

%% Average over subjects and conditions and remove any NaN (for non-included subjects and conditions)
if strcmp(analysis,'object_decoding')
    avg_over_conditions_all_subjects_cat = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_categorization,1),2),3)); 
    avg_over_conditions_all_subjects_dis = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_distraction,1),2),3)); 
elseif strcmp(analysis,'category_decoding')
    avg_over_conditions_all_subjects_cat = squeeze(nanmean(decoding_accuracies_all_subjects_categorization,1)); 
    avg_over_conditions_all_subjects_dis = squeeze(nanmean(decoding_accuracies_all_subjects_distraction,1)); 
end

%% Plot the difference
diff_matrix = avg_over_conditions_all_subjects_cat - avg_over_conditions_all_subjects_dis;
h = pcolor(diff_matrix); 
set(gcf, 'Position', get(0, 'Screensize'));
set(h, 'EdgeColor', 'none');
axis square;
plot_title = sprintf('Difference in time-generalized %s decoding between categorization and distraction tasks (N=%d)',analysis_name,numel(subjects));
title(plot_title);
cbar = colorbar;
ylabel('Timepoints trained on');
xlabel('Timepoints tested on');
ylabel(cbar,'Decoding accuracy (%)');

%% Save the matrices and plot
save(fullfile(results_avg_dir,sprintf('timegen_svm_%s_decoding_subjects_%d_%d_categorization',...
    analysis_name,subjects(1),subjects(end))),'avg_over_conditions_all_subjects_cat');
save(fullfile(results_avg_dir,sprintf('timegen_svm_%s_decoding_subjects_%d_%d_fixation',...
    analysis_name,subjects(1),subjects(end))),'avg_over_conditions_all_subjects_dis');
save(fullfile(results_avg_dir,sprintf('diff_matrix_timegen_svm_%s_decoding_subjects_%d_%d',...
    analysis_name,subjects(1),subjects(end))),'diff_matrix');
saveas(gcf,fullfile(results_avg_dir,sprintf('diff_matrix_timegen_svm_%s_decoding_subjects_%d_%d',analysis_name,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('diff_matrix_timegen_svm_%s_decoding_subjects_%d_%d.png',analysis_name,subjects(1),subjects(end)))); %save as png
close(gcf)



