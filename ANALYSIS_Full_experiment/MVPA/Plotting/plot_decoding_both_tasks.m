function plot_decoding_both_tasks(subjects,with_stats,analysis)
%PLOT_DECODING_BOTH_TASKS Plot the results from object decoding, averaged over
%all participants for both tasks (categorization and distraction).
%
%Input: subject IDs, with_stats (1 plot with stats, 0 plot without),
%analysis ('object_decoding' or 'category_decoding')
%
%Output: curve of decoding accuracies per timepoint, for two tasks
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate
numConditions = 60;
numTimepoints = 200;
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects_cat = NaN(max(subjects),numConditions,numConditions,numTimepoints);
    decoding_accuracies_all_subjects_dis = NaN(max(subjects),numConditions,numConditions,numTimepoints);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects_cat = NaN(max(subjects),numTimepoints);
    decoding_accuracies_all_subjects_dis = NaN(max(subjects),numTimepoints);
end

%% Loop: collect results from all subjects
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        cat_filename = 'svm_decoding_accuracy_categorization.mat';
        dis_filename = 'svm_decoding_accuracy_fixation.mat';
        all_dimensions = repmat({':'},1,3); %conditions x conditions x timepoints
    elseif strcmp(analysis,'category_decoding')
        cat_filename = 'cross_validated_dth_pseudotrials_svm_decodingAccuracy_categorization.mat';
        dis_filename = 'cross_validated_dth_pseudotrials_svm_decodingAccuracy_fixation.mat';
        all_dimensions = {':'}; % timepoints
    end
    load(fullfile(subject_results_dir,cat_filename),'decodingAccuracy_avg');
    decoding_accuracies_all_subjects_cat(subject,all_dimensions{:}) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,dis_filename),'decodingAccuracy_avg');
    decoding_accuracies_all_subjects_dis(subject,all_dimensions{:}) = decodingAccuracy_avg;
end

%% For stats matrix: num subjects x num timepoints (average over conditions)
if strcmp(analysis,'object_decoding')
    for_stats_cat = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects_cat,2),3));
    for_stats_dis = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects_dis,2),3));
elseif strcmp(analysis,'category_decoding')
    for_stats_cat = decoding_accuracies_all_subjects_cat;
    for_stats_dis = decoding_accuracies_all_subjects_dis;
end

%% Difference between tasks
for_stats_cat = for_stats_cat(subjects,:);
for_stats_dis = for_stats_dis(subjects,:);
for_stats_diff = for_stats_cat-for_stats_dis;
diff_curve = squeeze(nanmean(for_stats_diff,1));

%% Average over subjects 
avg_over_conditions_all_subjects_cat = squeeze(nanmean(for_stats_cat,1));
avg_over_conditions_all_subjects_dis = squeeze(nanmean(for_stats_dis,1));

%% Setup the figure
f = figure(abs(round(randn*10))); %Random figure number
f.Position(3:4) = f.Position(3:4)*1.5;
legend_plot = cell(3,1);

%% Plot
%loop over both tasks
for task = 1:3
    %invert the task indices so the categorization task is plotted last
    if task == 1
        task_plot = 2;
    elseif task == 2
        task_plot = 1;
    elseif task == 3
        task_plot = 3;
    end
    if task_plot == 1
        task_name = 'categorization';
        data = avg_over_conditions_all_subjects_cat-50;
        plot_location = -5;
        color_err = [0.75 0.75 0.95]; %error bars color: can't be the same one as that of the plot
        color_data = [0 0.4 0.85];    
        for_stats_data = for_stats_cat-50;
    elseif task_plot == 2
        task_name = 'fixation';
        data = avg_over_conditions_all_subjects_dis-50;
        plot_location = -7;
        color_err = [0.95 0.75 0.9]; 
        color_data = [0.9 0.2 0.8];   
        for_stats_data = for_stats_dis-50;
    elseif task_plot == 3
        task_name = 'difference';
        data = diff_curve;
        plot_location = -9;
        color_err = [0.5 0.5 0.5];%/255;
        color_data = 'k';
        for_stats_data = for_stats_diff;
    end

    if ~with_stats
        %plot the data
        p = plot(data,'Linewidth',3, 'Color', color_data);
        hold on;
        if task_plot==1
            p1 = p; %for the legend
            legend_plot{task_plot} = 'Scene categorization';
        elseif task_plot == 2
            p2 = p;
            legend_plot{task_plot} = 'Distraction';
        elseif task_plot == 3
            p3 = p;
            legend_plot{task_plot} = 'Scene categorization-distraction';
        end
    else
        %Stat parameters
        filename = fullfile(results_avg_dir,...
            sprintf('stats_fdr_%s_%s_subjects_%d_%d.mat',analysis,task_name,subjects(1),subjects(end)));
        if exist(filename,'file')
            load(filename,'stats_decoding');
        else
            peak_latency_true = find(data==max(data));
            stats_decoding.peak_latency_ground_ms = (peak_latency_true-40)*5;
            stats_decoding.num_perms = 1000;
            stats_decoding.tail = 'right';
            stats_decoding.qvalue = 0.01;
            [stats_decoding.significant_timepoints,stats_decoding.pvalues,...
                stats_decoding.crit_p, stats_decoding.adjusted_pvalues]...
                = fdr_permutation_cluster_1sample_alld(for_stats_data,...
                stats_decoding.num_perms,stats_decoding.tail,stats_decoding.qvalue);
            [stats_decoding.peak_latency_bs, stats_decoding.CI] = bootstrap_peak_latency(for_stats_data);
            save(filename,'stats_decoding');
        end

        
        %1) significant timepoints
        st = (stats_decoding.significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color_data); 
        hold on;
        
        %2) error bars
        stdDM = std(for_stats_data); 
        err = stdDM/sqrt(size(for_stats_data,1)); %standard deviation/sqrt of num subjects  

        %plot as a shaded area
        top_curve = data + err;
        bottom_curve = data - err;
        x2 = [1:numTimepoints, fliplr(1:numTimepoints)];
        shaded_area = [top_curve, fliplr(bottom_curve)];
        fill(x2, shaded_area, color_err,'FaceAlpha',0.5);
        hold on;
        
        %3) plot the data - on top of the error bars
        p = plot(data,'Linewidth',3, 'Color', color_data);
        hold on;
        if task_plot==1
            p1 = p; %for the legend
            legend_plot{task_plot} = 'Scene categorization';
        elseif task_plot == 2
            p2 = p;
            legend_plot{task_plot} = 'Distraction';
        elseif task_plot == 3
            p3 = p;
            legend_plot{task_plot} = 'Scene categorization-distraction';
        end

        %4) peak latency and 95% confidence interval 
        if strcmp(analysis,'object_decoding')
            height_cat = 35;
            height_dis = 32.5;
            height_diff = 30;
        elseif strcmp(analysis,'category_decoding')
            height_cat = 30;
            height_dis = 27.5;
            height_diff = 25;
        end
        if task_plot == 1
            height = height_cat;
        elseif task_plot == 2 
            height = height_dis;
        elseif task_plot == 3
            height = height_diff;
        end
        
        %plot arrows and lines for peak latency and CI
        fontsize = 20;
        peak_acc = max(data);
        str_pl = num2str(stats_decoding.peak_latency_ground_ms);
        quiver(stats_decoding.CI(1),height,0,-2,0,'Color',color_data,'ShowArrowHead','off','LineStyle',':','LineWidth',2); 
        quiver(stats_decoding.CI(2),height,0,-2,0,'Color',color_data,'ShowArrowHead','off','LineStyle',':','LineWidth',2);
        if task_plot == 2
            arrow_x = peak_latency_true+8;
            text(arrow_x,peak_acc,['\leftarrow',str_pl,' ms'],'Color',color_data,'FontSize',fontsize);
        elseif task_plot == 1
            arrow_x = peak_latency_true-50;
            text(arrow_x,peak_acc,[str_pl,' ms \rightarrow'],'Color',color_data,'FontSize',fontsize);
        elseif task_plot == 3
            arrow_x = peak_latency_true-50;
            text(arrow_x,peak_acc,[str_pl,' ms \rightarrow'],'Color',color_data,'FontSize',fontsize);    
        end
        
    end
end 

%% Plot parameters
plot_title = sprintf('%s over time (N=%d)',analysis,numel(subjects));
legend_bool = 0;
title_bool = 0;
plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,12,'best','Decoding accuracy-50 (%)');
if legend_bool==1
    legend([p1,p2,p3],legend_plot);
end
set(gca,'FontName','Arial','FontSize',18);
ylim([-10,40]);

%% Save the plot and matrices
save(fullfile(results_avg_dir,sprintf('for_stats_cat_svm_%s_subjects_%d_%d_both_tasks.mat',analysis,subjects(1),subjects(end))),'for_stats_cat'); 
save(fullfile(results_avg_dir,sprintf('for_stats_dis_svm_%s_subjects_%d_%d_both_tasks.mat',analysis,subjects(1),subjects(end))),'for_stats_dis'); 
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_both_tasks',analysis,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_both_tasks.svg',analysis,subjects(1),subjects(end)))); %save as svg
close(gcf);    

end