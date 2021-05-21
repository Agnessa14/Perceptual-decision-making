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
    decoding_accuracies_all_subjects_dis = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects_cat = NaN(sorted_subjects(end),numTimepoints);
    decoding_accuracies_all_subjects_dis = NaN(sorted_subjects(end),numTimepoints);
end

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        cat_filename = 'svm_decoding_accuracy_categorization.mat';
        dis_filename = 'svm_decoding_accuracy_fixation.mat';
        all_dimensions = repmat({':'},1,3); %conditions x conditions x timepoints
    elseif strcmp(analysis,'category_decoding')
        cat_filename = 'svm_artificial_vs_natural_decoding_accuracy_categorization.mat';
        dis_filename = 'svm_artificial_vs_natural_decoding_accuracy_fixation.mat';
        all_dimensions = {':'}; % timepoints
    end
    load(fullfile(subject_results_dir,cat_filename));
    decoding_accuracies_all_subjects_cat(subject,all_dimensions{:}) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,dis_filename));
    decoding_accuracies_all_subjects_dis(subject,all_dimensions{:}) = decodingAccuracy_avg;
end   

%% Average over subjects + conditions and remove any NaN (for non-included subjects)
if strcmp(analysis,'object_decoding')
    avg_over_conditions_all_subjects_cat = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_cat,1),2),3))';
    avg_over_conditions_all_subjects_dis = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects_dis,1),2),3))';
elseif strcmp(analysis,'category_decoding')
    avg_over_conditions_all_subjects_cat = squeeze(nanmean(decoding_accuracies_all_subjects_cat,1));
    avg_over_conditions_all_subjects_dis = squeeze(nanmean(decoding_accuracies_all_subjects_dis,1));
end

%% Setup the figure
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));

%% Plot stats if needed
if with_stats
    num_perms = 10000;
    cluster_th = 0.05;
    significance_th = 0.05;
    
    for task = 1:2
        %invert the task indices so the categorization task is plotted last
        if task == 1
            task_plot = 2;
        elseif task == 2
            task_plot = 1;
        end
        if task_plot == 1
            task_name = 'categorization';
            data = avg_over_conditions_all_subjects_cat;
            plot_location = 47.5;
            color_err = [0.75 0.75 0.95]; %can't be the same one as that of the plot
            color_data = [0 0.4 0.85];     
        elseif task_plot == 2
            task_name = 'fixation';
            data = avg_over_conditions_all_subjects_dis;
            plot_location = 46;
            color_err = [0.95 0.75 0.9]; 
            color_data = [0.9 0.2 0.8];       
        end
        
        %%error bars
        filename_forstats = fullfile(results_avg_dir,sprintf('for_stats_subjects_%d_%d_%s_task_%s.mat',...
        subjects(1),subjects(end),task_name,analysis));
        if exist(filename_forstats,'file')
            load(filename_forstats);
        else
            for_stats = all_subjects_for_stats(subjects,task_plot,analysis);
        end
        stdDM = std(for_stats); 
        err = stdDM/sqrt(size(for_stats,1)); %standard deviation/sqrt of num subjects  
        
        %plot as a shaded area
        top_curve = data + err;
        bottom_curve = data - err;
        x2 = [1:numTimepoints, fliplr(1:numTimepoints)];
        shaded_area = [top_curve, fliplr(bottom_curve)];
        fill(x2, shaded_area, color_err);
        hold on;

        %plot the data
        p = plot(data,'Linewidth',3, 'Color', color_data);
        if task_plot==1
            p1 = p; %for the legend
        elseif task_plot == 2
            p2 = p;
        end
        
        %significant timepoints
        filename_sign = fullfile(results_avg_dir,sprintf('significant_timepoints_subjects_%d_%d_%s_task_%s.mat',...
        subjects(1),subjects(end),task_name,analysis));
        if exist(filename_sign,'file')
            load(filename_sign);
        else
            significant_timepoints = run_permutation_stats(subjects,task_plot,analysis,...
                for_stats,num_perms,cluster_th,significance_th);
        end
        st = (significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color_data); 
        
        %peak latency and 95% confidence interval 
        [peak_latency, CI] = bootstrap_peak_latency(subjects,task_plot,analysis);
        if strcmp(analysis,'object_decoding')
            arrow_x = 51;
            height_cat = 82.5;
            height_dis = 80;
        elseif strcmp(analysis,'category_decoding')
            arrow_x = 60;
            height_cat = 78;
            height_dis = 75.5;
        end
        if task_plot == 1
            height = height_cat;
        elseif task_plot == 2 
            height = height_dis;
        end
        quiver(arrow_x,data(peak_latency),6,0,0,'Color',color_data,'ShowArrowHead','on','MaxHeadSize',0.75,'LineWidth',2);
        quiver(CI(1),height,0,-2,0,'Color',color_data,'ShowArrowHead','off','LineStyle',':','LineWidth',2); 
        quiver(CI(2),height,0,-2,0,'Color',color_data,'ShowArrowHead','off','LineStyle',':','LineWidth',2);
        save(fullfile(results_avg_dir,sprintf('peak_latency_subjects_%d_%d_%s_task_%s.mat',...
        subjects(1),subjects(end),task_name,analysis)),'peak_latency');
        save(fullfile(results_avg_dir,sprintf('confidence_interval_peak_latency_subjects_%d_%d_%s_task_%s.mat',...
        subjects(1),subjects(end),task_name,analysis)),'CI');    
    end
end 

%Plot parameters
if strcmp(analysis,'object_decoding')
    analysis_title = 'Object';
elseif strcmp(analysis,'category_decoding')
    analysis_title = 'Category';
end
plot_title = sprintf('%s decoding over time (N=%d)',analysis_title,numel(subjects));
onset_time = 40; 
xticks(0:10:200);
legend_cell = {'Scene categorization','Distraction'}; %can figure out a way to add the arrowws and CIs to the legend
plotting_parameters(plot_title,'',onset_time,12,'best','Decoding accuracy (%)'); %[0.75 0.7 0.1 0.1]
legend([p1,p2],legend_cell);

%% Save the plot
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_both_tasks',analysis,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_both_tasks.svg',analysis,subjects(1),subjects(end)))); %save as svg
close(gcf);    

%% Plot the difference curve
diff_curve = avg_over_conditions_all_subjects_cat-avg_over_conditions_all_subjects_dis;
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
plot(diff_curve,'--','LineWidth',3,'Color','k');
diff_title = sprintf('%s decoding: difference in accuracy between the categorization and the distraction tasks over time (N=%d)',analysis_title,numel(subjects));
plotting_parameters(diff_title,'',onset_time,12,'best','%'); 
legend('off');

%Save
saveas(gcf,fullfile(results_avg_dir,sprintf('diff_curve_svm_%s_subjects_%d_%d',analysis,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('diff_curve_svm_%s_subjects_%d_%d.svg',analysis,subjects(1),subjects(end)))); %save as svg
close(gcf);    

end