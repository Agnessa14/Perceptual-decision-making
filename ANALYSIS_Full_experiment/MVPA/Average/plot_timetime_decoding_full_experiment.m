function plot_timetime_decoding_full_experiment(subjects,task,analysis,with_stats)
%PLOT_TIMETIME_DECODING_FULL_EXPERIMENT Plot the results from time-generalized object decoding, averaged over
%all participants.
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

%% Define some variables
numTimepoints = 50;
numConditions = 60;

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

%Check if the matrix of averaged accuracies exists
filename_data = fullfile(results_avg_dir,sprintf('timegen_svm_%s_subjects_%d_%d_%s.mat',analysis,...
    subjects(1),subjects(end),task_name));
filename_for_stats = fullfile(results_avg_dir,sprintf('for_stats_timegen_svm_%s_subjects_%d_%d_%s.mat',analysis,...
    subjects(1),subjects(end),task_name));
if exist(filename_data,'file') && exist(filename_for_stats,'file')
    load(filename_data,'avg_over_conditions_all_subjects');
    load(filename_for_stats,'for_stats');
else
    %Preallocate
    if strcmp(analysis,'object_decoding')
        decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints,numTimepoints);
        all_dimensions = repmat({':'},1,4);
    elseif strcmp(analysis,'category_decoding')
        decoding_accuracies_all_subjects = NaN(max(subjects),numTimepoints,numTimepoints);
        all_dimensions = repmat({':'},1,2);
    end
    
    %Loop over subjects
    for subject = subjects
        subname = get_subject_name(subject);
        subject_results_dir = fullfile(results_dir,subname);
        if strcmp(analysis,'object_decoding')
            filename = 'time_generalized_svm_object_decoding';
        elseif strcmp(analysis,'category_decoding')
            filename = 'time_gen_svm_artificial_vs_natural_decoding_accuracy';
        end
        
        %for cross task, add a prefix
        if task == 3
            filename = sprintf('2_models_%s',filename);
        end
        
        %load & fill the preallocated matrix
        load(fullfile(subject_results_dir,sprintf('%s_%s.mat',filename,task_name)),'timeg_decodingAccuracy_avg');
        decoding_accuracies_all_subjects(subject,all_dimensions{:}) = timeg_decodingAccuracy_avg;
    end
    
    %% Average over conditions (use for stats) and remove any NaN (for non-included conditions)
    if strcmp(analysis,'object_decoding')
        for_stats = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects(subjects,:,:,:,:),2),3));
    elseif strcmp(analysis,'category_decoding')
        for_stats = decoding_accuracies_all_subjects(subjects,:,:,:,:);
    end
    save(filename_for_stats,'for_stats');

    %% Average over subjects, remove NaNs
    avg_over_conditions_all_subjects = squeeze(mean(for_stats,1));
    
    %% Save the matrix
    save(filename_data,'avg_over_conditions_all_subjects');

end   

%% Plot the average of all subjects
h = pcolor(avg_over_conditions_all_subjects-50);
set(h, 'EdgeColor', 'none');
axis square;
hold on;
if task<3
    plot_title = sprintf('Time-generalized %s in a %s task (N=%d)',analysis,task_name_title,numel(subjects));
elseif task == 3
    plot_title = sprintf('Time-generalized %s across categorization and distraction tasks (N=%d)',analysis,numel(subjects));
end
title_bool = 0;
if title_bool==1
    title(plot_title);
end
cbar = colorbar;
if task < 3
    ylabel('Timepoints trained on');
    xlabel('Timepoints tested on');
elseif task == 3
    ylabel({'Timepoints trained on:';'Categorization task'});
    xlabel({'Timepoints tested on:'; 'Distraction task'});
end
ylabel(cbar,'Decoding accuracy-50 (%)');
caxis([-5 20]);
xticks(0:5:50);
xticklabels(-200:100:800);
yticks(0:5:50);
xline(10,'--','Color','w');
yline(10,'--','Color','w');
yticklabels(-200:100:800);
set(gca,'FontName','Arial','FontSize',18)

%Plot a black line over the diagonal for cross-task
if task == 3
    plot(1:numTimepoints,1:numTimepoints,'--','LineWidth',2.5,'Color','k');
end

%% Plot stats if needed
if with_stats

    stats_decoding.num_perms = 1000;
    stats_decoding.qvalue = 0.01;
    stats_decoding.tail = 'right';
    filename = fullfile(results_avg_dir,...
        sprintf('stats_fdr_timetime_%s_%s_subjects_%d_%d.mat',analysis,task_name,subjects(1),subjects(end)));
    if exist(filename,'file')
        load(filename,'stats_decoding');
    else
        [stats_decoding.significant_timepoints,stats_decoding.pvalues,...
            stats_decoding.crit_p,stats_decoding.adjusted_pvalues]...
            = fdr_permutation_cluster_1sample_alld(for_stats-50,...
            stats_decoding.num_perms,stats_decoding.tail,stats_decoding.qvalue);
        save(filename,'stats_decoding');
    end

    %plot the stats
    contour(stats_decoding.significant_timepoints,1,'LineColor','w','LineWidth',2);

end

%% Save the plot
saveas(gcf,fullfile(results_avg_dir,sprintf('timegen_svm_%s_subjects_%d_%d_%s',analysis,subjects(1),subjects(end),task_name))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('timegen_svm_%s_subjects_%d_%d_%s.png',analysis,subjects(1),subjects(end),task_name))); %save as matlab figure
close(gcf)

end

