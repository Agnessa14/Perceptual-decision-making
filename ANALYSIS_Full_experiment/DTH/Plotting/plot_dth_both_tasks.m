function plot_dth_both_tasks(subjects,with_cross_task,with_stats) %distance art, distance nat, RT
%PLOT_DTH_BOTH_TASKS Plot the fixed effects distance-to-hyperplane results for both
%tasks. Assumes that the correlations & stats have already been computed.
%
%Input: subjects' ID (e.g., 1:13), add the cross task curves (1 yes, 0 no),
%add stats (1 yes, 0 no)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%Author: Agnessa Karapetian,2021
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%Load the correlations
load(fullfile(results_dir,sprintf...
    ('cv_all_dist_med_rt_dth_subjects_%d_%d_categorization.mat',subjects(1),subjects(end))),'dth_results');
corr_both_categorization = dth_results.corr_both_categories;
clear dth_results;
load(fullfile(results_dir,sprintf...
    ('cv_all_dist_med_rt_dth_subjects_%d_%d_fixation.mat',subjects(1),subjects(end))),'dth_results');
corr_both_distraction = dth_results.corr_both_categories;
clear dth_results;

if with_cross_task  
    load(fullfile(results_dir,sprintf...
        ('cv_all_dist_med_rt_dth_subjects_%d_%d_cross_task_categorization_distances.mat',subjects(1),subjects(end))),'dth_results');
    corr_both_cross_task_categorization_dist = dth_results.corr_both_categories;
    clear dth_results;
    load(fullfile(results_dir,sprintf...
        ('cv_all_dist_med_rt_dth_subjects_%d_%d_cross_task_fixation_distances.mat',subjects(1),subjects(end))),'dth_results');
    corr_both_cross_task_distraction_dist = dth_results.corr_both_categories;
    clear dth_results;
end

%% Plot 
cmap_1 = autumn;
cmap_2 = winter;

color_cat = cmap_1(100,:);
color_dis = cmap_2(100,:);
if with_cross_task
    color_cross_1 = cmap_2(175,:);
    color_cross_2 = cmap_1(175,:);
end
figure(abs(round(randn*10)));
categ = plot(corr_both_categorization,'LineWidth',2,'Color',color_cat);
hold on;
distr = plot(corr_both_distraction,'LineWidth',2,'Color',color_dis);

if with_cross_task
    categ_eeg = plot(corr_both_cross_task_categorization_dist,'LineWidth',2,'Color',color_cross_1);
    hold on;
    distr_eeg = plot(corr_both_cross_task_distraction_dist,'LineWidth',2,'Color',color_cross_2);
end

%% Plot stats if needed
if with_stats
    analysis = 'random_dth';
    if with_cross_task
        num_tasks = 4;
    else
        num_tasks = 2;
    end
    for task = 1:num_tasks
        %Stat plot parameters
        if task<3
            task_name = get_task_name(task);
        else
            task_name = get_task_name(task-2); %for iteration 3: categ, iter 4: fix (distance task name)
        end
        
        if task == 1
            plot_location = 0.15;
            color = color_cat;
        elseif task == 2
            plot_location = 0.18;
            color = color_dis;
        elseif task == 3
            plot_location = 0.19;
            color = color_cross_1;
        elseif task == 4
            plot_location = 0.17;
            color = color_cross_2;
        end
   
        %Load and plot the stats
        if task<3
            filename_sign = 'dth_permutation_stats';
        elseif (task<=4)&&(task>2)
            filename_sign = 'dth_permutation_stats_crosstask';
        end               
     
        filename = fullfile(results_dir,sprintf('%s_%d_%d_%s_task_%s_both.mat',filename_sign,...
            subjects(1),subjects(end),task_name,analysis));        
        load(filename,'permutation_stats'); 
        st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;              
    end
end 

%% Plotting parameters
plot_title =  sprintf('Correlation between the distance to hyperplane and reaction time in 60 scenes (N=%d)',numel(subjects));
if with_cross_task
    legend_plot = {'Scene categorization','Distraction EEG, categorization RTs',...
        'Distraction','Scene categorization EEG, distraction RTs'};
else
    legend_plot = {'Scene categorization','Distraction'};
end
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
ylim([-0.3 0.3]);
yticks(-0.3:0.1:0.3);
legend_bool = 0;
title_bool = 0;
plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Spearman''s coefficient'); 
if with_cross_task && legend_bool==1
    legend([categ,distr_eeg,distr,categ_eeg],legend_plot);
end

%% Save
if with_cross_task
    saveas(gcf,fullfile(results_dir,sprintf('cv_all_dist_med_rt_dth_subjects_%d_%d_both_tasks_with_crosstask',subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_dir,sprintf('cv_all_dist_med_rt_dth_subjects_%d_%d_both_tasks_with_crosstask.svg',subjects(1),subjects(end))));
else
    saveas(gcf,fullfile(results_dir,sprintf('cv_all_dist_med_rt_dth_subjects_%d_%d_both_tasks',subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_dir,sprintf('cv_all_dist_med_rt_dth_subjects_%d_%d_both_tasks.svg',subjects(1),subjects(end))));
end
close(gcf);
end