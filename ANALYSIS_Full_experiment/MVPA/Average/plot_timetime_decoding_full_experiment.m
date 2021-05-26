function plot_timetime_decoding_full_experiment(subjects,task,analysis,with_stats)
%PLOT_TIMETIME_DECODING_FULL_EXPERIMENT Plot the results from time-generalized object decoding, averaged over
%all participants and subject-specific.
%
%Input: subject IDs (e.g., 1:13), task (1=categorization, 2=distraction,3=cross-task),
%analysis('object_decoding' or 'category_decoding'), with stats (1) or
%without (0)
%
%Output: 2D heatmap of decoding accuracies for each pair of timepoints
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
    decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints,numTimepoints);
    all_dimensions = repmat({':'},1,4);
    analysis_name = analysis(1:6);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects = NaN(max(subjects),numTimepoints,numTimepoints);
    all_dimensions = repmat({':'},1,2);
    analysis_name = analysis(1:8);
end

%task name
if task < 3
    task_name = get_task_name(task);
    if task == 1
        task_name_title = task_name;
    else
        task_name_title = 'distraction';
    end
elseif task == 3
    if strcmp(analysis,'object_decoding')
        task_name = 'crosstask';
    elseif strcmp(analysis,'category_decoding')
        task_name = 'cross_task';
    end
end

%Loop over subjects
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        load(fullfile(subject_results_dir,sprintf('time_generalized_svm_object_decoding_%s.mat',task_name)));
    elseif strcmp(analysis,'category_decoding')
        load(fullfile(subject_results_dir,sprintf('3PT_time_gen_svm_artificial_vs_natural_decoding_accuracy_%s.mat',task_name)));
    end    
    decoding_accuracies_all_subjects(subject,all_dimensions{:}) = timeg_decodingAccuracy_avg;
end   

%% Average over subjects and conditions and remove any NaN (for non-included subjects and conditions)
if strcmp(analysis,'object_decoding')
    avg_over_conditions_all_subjects = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects,1),2),3)); 
elseif strcmp(analysis,'category_decoding')
    avg_over_conditions_all_subjects = squeeze(nanmean(decoding_accuracies_all_subjects,1)); 
end

%% Plot the average of all subjects
h = pcolor(avg_over_conditions_all_subjects); 
set(gcf, 'Position', get(0, 'Screensize'));
set(h, 'EdgeColor', 'none');
axis square;
if task<3
    plot_title = sprintf('Time-generalized %s decoding in a %s task (N=%d)',analysis_name,task_name_title,numel(subjects));
elseif task == 3 
    plot_title = sprintf('Time-generalized %s decoding across categorization and distraction tasks (N=%d)',analysis_name,numel(subjects));
end
title(plot_title);
cbar = colorbar;
if strcmp(analysis,'object_decoding')
    ylabel('Timepoints trained on');
    xlabel('Timepoints tested on');
elseif strcmp(analysis,'category_decoding')
    ylabel('Timepoints trained on: Categorization task');
    xlabel('Timepoints tested on: Distraction task');
end

ylabel(cbar,'Decoding accuracy (%)');

%% Save the plot
% saveas(gcf,fullfile(results_avg_dir,sprintf('timegen_svm_%s_decoding_subjects_%d_%d_%s',analysis_name,subjects(1),subjects(end),task_name))); %save as matlab figure
% saveas(gcf,fullfile(results_avg_dir,sprintf('timegen_svm_%s_decoding_subjects_%d_%d_%s.png',analysis_name,subjects(1),subjects(end),task_name))); %save as matlab figure
close(gcf)

%% Plot stats if needed
if with_stats
    num_perms = 1000;
    cluster_th = 0.1;
    significance_th = 0.1;
    %load the subjects x timepoints matrix
    filename_forstats = fullfile(results_avg_dir,sprintf('for_stats_subjects_%d_%d_%s_time_%s.mat',...
    subjects(1),subjects(end),task_name,analysis));
    if exist(filename_forstats,'file')
        load(filename_forstats);
    else
        for_stats = all_subjects_for_stats(subjects,task,sprintf('time_%s',analysis));
    end
    
    %load the stats if they have been run already
    filename_sign = fullfile(results_avg_dir,sprintf('significant_timepoints_subjects_%d_%d_%s_%s_time_%s.mat',...
    subjects(1),subjects(end),task_name,analysis));
    if exist(filename_sign,'file')
        load(filename_sign);
    else
        significant_timepoints = run_permutation_stats(subjects,task,sprintf('time_%s',analysis),for_stats,num_perms,cluster_th,significance_th);
    end
    
    %plot the stats - different plot from the RSA
    significant_timepoints(isnan(significant_timepoints))=0;
    h = pcolor(significant_timepoints); 
    set(gcf, 'Position', get(0, 'Screensize'));
    set(h, 'EdgeColor', 'none');
    axis square;
    cbar = colorbar;
    ylabel(cbar,'Significance');
    if strcmp(analysis,'object_decoding')
        ylabel('Timepoints trained on');
        xlabel('Timepoints tested on');
    elseif strcmp(analysis,'category_decoding')
        ylabel('Timepoints trained on: Categorization task');
        xlabel('Timepoints tested on: Distraction task');
    end

    title(sprintf('Statistical analysis of time-generalized RSA of scene processing in categorization and distraction tasks (N=%d)',numel(subjects)));  
    
    %save the plot
    saveas(gcf,fullfile(results_avg_dir,sprintf('stats_time_%s_%s_subjects_%d_%d',task_name,analysis,subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_avg_dir,sprintf('stats_time_%s_%s_subjects_%d_%d.png',task_name,analysis,subjects(1),subjects(end)))); 
    close(gcf);
end 

end

