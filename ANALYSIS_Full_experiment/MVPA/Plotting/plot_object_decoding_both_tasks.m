function plot_object_decoding_both_tasks(subjects,with_stats,analysis)
%PLOT_OBJECT_BOTH_TASKS Plot the results from object decoding, averaged over
%all participants for both tasks (categorization and distraction).
%
%Input: subject IDs, with_stats (1 plot with stats, 0 plot without),
%analysis ('object_decoding' or 'category_decoding')
%
%Output: curve of decoding accuracies per timepoint, for two tasks
%
%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate
numConditions = 60;
numTimepoints = 200;
sorted_subjects = sort(subjects);
decoding_accuracies_all_subjects_cat = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
decoding_accuracies_all_subjects_fix = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'svm_decoding_accuracy_categorization.mat'));
    decoding_accuracies_all_subjects_cat(subject,:,:,:) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,'svm_decoding_accuracy_fixation.mat'));
    decoding_accuracies_all_subjects_fix(subject,:,:,:) = decodingAccuracy_avg;
end   

%% Remove any NaN (for non-included subjects)
avg_over_conditions_all_subjects_cat = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_cat,1),2),3));
avg_over_conditions_all_subjects_fix = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_fix,1),2),3));

%% Plot the average of all subjects
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
color_cat = 'b';
color_fix = 'm';
plot(avg_over_conditions_all_subjects_cat,'Linewidth',3, 'Color', color_cat);
hold on;
plot(avg_over_conditions_all_subjects_fix,'Linewidth',3, 'Color', color_fix);
hold on;
title = sprintf('Object decoding per timepoint for %d subjects',numel(subjects));
onset_time = 40; 
xticks(0:10:200);
ylim([45 80]);

%% Plot stats if needed
if with_stats
    for task = 1:2
        if task == 1
            task_name = 'categorization';
            plot_name = avg_over_conditions_all_subjects_cat;
            plot_location = 47.5;
            color = color_cat;
        elseif task == 2
            task_name = 'fixation';
            plot_name = avg_over_conditions_all_subjects_fix;
            plot_location = 46;
            color = color_fix;
        end
        %error bars
        load(fullfile(results_avg_dir,sprintf('for_stats_%d_subjects_%s_task_%s',...
        numel(subjects),task_name,analysis)));
        stdDM = std(for_stats); %std(26x200)
        err = stdDM/sqrt(size(for_stats,1)); %standard deviation/sqrt of num subjects
        errorbar(plot_name, err, 'Color',color); %plot
        hold on;

        %significant timepoints
        load(fullfile(results_avg_dir,sprintf('significant_timepoints_%d_subjects_%s_task_%s',...
        numel(subjects),task_name,analysis)));
        st = (significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;
        
        %peak latency and 95% confidence interval 
        [peak_latency, CI] = bootstrap_peak_latency(subjects,task,analysis);
        quiver(peak_latency,80,0,-5,0,'Color',color,'ShowArrowHead','on','MaxHeadSize',1,'LineWidth',2) %kind of ugly arrow..check the mathworks page
        quiver(CI(1),80,0,-4,0,'Color',color,'ShowArrowHead','off','LineStyle',':','LineWidth',2); 
        quiver(CI(2),80,0,-4,0,'Color',color,'ShowArrowHead','off','LineStyle',':','LineWidth',2); 
    end
end 

legend_cell = {'Scene categorization','Distraction'}; %can figure out a way to add the arrowws and CIs to the legend
plotting_parameters(title,legend_cell{1:3},onset_time,12,'best','Decoding accuracy (%)'); %[0.75 0.7 0.1 0.1]

%% Save the plot
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_object_decoding_%d_subjects_both_tasks',numel(subjects)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_object_decoding_%d_subjects_both_tasks.svg',numel(subjects)))); %save as svg
close(gcf);    

end