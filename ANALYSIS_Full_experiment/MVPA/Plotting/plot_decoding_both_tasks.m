function plot_decoding_both_tasks(subjects,with_stats,analysis)
%PLOT_DECODING_BOTH_TASKS Plot the results from object decoding, averaged over
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
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects_cat = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
    decoding_accuracies_all_subjects_fix = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects_cat = NaN(sorted_subjects(end),numTimepoints);
    decoding_accuracies_all_subjects_fix = NaN(sorted_subjects(end),numTimepoints);
end

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        cat_filename = 'svm_decoding_accuracy_categorization.mat';
        fix_filename = 'svm_decoding_accuracy_fixation.mat';
        all_dimensions = repmat({':'},1,3); %conditions x conditions x timepoints
    elseif strcmp(analysis,'category_decoding')
        cat_filename = 'svm_artificial_vs_natural_decoding_accuracy_categorization.mat';
        fix_filename = 'svm_artificial_vs_natural_decoding_accuracy_fixation.mat';
        all_dimensions = ':'; % timepoints
    end
    load(fullfile(subject_results_dir,cat_filename));
    decoding_accuracies_all_subjects_cat(subject,all_dimensions{:}) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,fix_filename));
    decoding_accuracies_all_subjects_fix(subject,all_dimensions{:}) = decodingAccuracy_avg;
end   

%% Average over subjects + conditions and remove any NaN (for non-included subjects)
if strcmp(analysis,'object_decoding')
    avg_over_conditions_all_subjects_cat = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_cat,1),2),3));
    avg_over_conditions_all_subjects_fix = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_fix,1),2),3));
elseif strcmp(analysis,'category_decoding')
    avg_over_conditions_all_subjects_cat = squeeze(nanmean(decoding_accuracies_all_subjects_cat,1));
    avg_over_conditions_all_subjects_fix = squeeze(nanmean(decoding_accuracies_all_subjects_fix,1));
end

%% Plot the average of all subjects
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
color_cat = 'b';
color_fix = 'm';
plot(avg_over_conditions_all_subjects_cat,'Linewidth',3, 'Color', color_cat);
hold on
plot(avg_over_conditions_all_subjects_fix,'Linewidth',3, 'Color', color_fix);
if strcmp(analysis,'object_decoding')
    analysis_title = 'Object';
elseif strcmp(analysis,'category_decoding')
    analysis_title = 'Category';
end
plot_title = sprintf('%s decoding per timepoint for %d subjects',analysis_title,numel(subjects));
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
        filename_forstats = fullfile(results_avg_dir,sprintf('for_stats_subjects_%d_%d_%s_task_%s.mat',...
        subjects(1),subjects(end),task_name,analysis));
        if exist(filename_forstats,'file')
            load(filename_forstats);
        else
            for_stats = all_subjects_for_stats(subjects,task,analysis);
        end
        stdDM = std(for_stats); %std(26x200)
        err = stdDM/sqrt(size(for_stats,1)); %standard deviation/sqrt of num subjects
        errorbar(plot_name, err, 'Color',color); %plot
        hold on;

        %significant timepoints
        filename_sign = fullfile(results_avg_dir,sprintf('significant_timepoints_subjects_%d_%d_%s_task_%s.mat',...
        subjects(1),subjects(end),task_name,analysis));
        if exist(filename_sign,'file')
            load(filename_sign);
        else
            significant_timepoints = run_permutation_stats(subjects,task,analysis,for_stats);
        end
        st = (significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        
        %peak latency and 95% confidence interval 
        [peak_latency, CI] = bootstrap_peak_latency(subjects,task,analysis);
        if strcmp(analysis,'object_decoding')
            height = 80;
        elseif strcmp(analysis,'category_decoding')
            height = 77;
        end
        quiver(peak_latency,height,0,-4,0,'Color',color,'ShowArrowHead','on','MaxHeadSize',1,'LineWidth',2) %kind of ugly arrow..check the mathworks page
        quiver(CI(1),height,0,-4,0,'Color',color,'ShowArrowHead','off','LineStyle',':','LineWidth',2); 
        quiver(CI(2),height,0,-4,0,'Color',color,'ShowArrowHead','off','LineStyle',':','LineWidth',2); 
    end
end 

legend_cell = {'Scene categorization','Distraction'}; %can figure out a way to add the arrowws and CIs to the legend
plotting_parameters(plot_title,legend_cell,onset_time,12,'best','Decoding accuracy (%)'); %[0.75 0.7 0.1 0.1]

%% Save the plot
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_both_tasks',analysis,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_both_tasks.svg',analysis,subjects(1),subjects(end)))); %save as svg
close(gcf);    

end