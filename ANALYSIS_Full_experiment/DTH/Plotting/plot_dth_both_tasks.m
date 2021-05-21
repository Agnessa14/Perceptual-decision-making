function plot_dth_both_tasks(subjects,with_cross_task,with_stats,method) %distance art, distance nat, RT
%PLOT_DTH_BOTH_TASKS Plot the fixed effects distance-to-hyperplane results for both
%tasks.
%
%Input: subjects' ID (e.g., 1:13), add the cross task curves (1 yes, 0 no),
%add stats (1 yes, 0 no), method ('fixed' or 'random' effects analysis)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%Load the correlations
if strcmp(method,'fixed')
    filename = 'pseudotrials_SVM_DTH';
elseif strcmp(method,'random')
    filename = 'random_effects_dth';
end
load(fullfile(results_dir,sprintf('%s_subjects_%d_%d_categorization.mat',filename,subjects(1),subjects(end))));
corr_both_categorization = dth_results.corr_both_categories;
clear dth_results;
load(fullfile(results_dir,sprintf('%s_subjects_%d_%d_fixation.mat',filename,subjects(1),subjects(end))));
corr_both_distraction = dth_results.corr_both_categories;
clear dth_results;

if with_cross_task  
    load(fullfile(results_dir,sprintf('%s_subjects_%d_%d_cross_task_categorization_distances.mat',filename,subjects(1),subjects(end))));
    corr_both_cross_task_categorization_dist = dth_results.corr_both_categories;
    clear dth_results;
    load(fullfile(results_dir,sprintf('%s_subjects_%d_%d_cross_task_fixation_distances.mat',filename,subjects(1),subjects(end))));
    corr_both_cross_task_distraction_dist = dth_results.corr_both_categories;
    clear dth_results;
end

%% Plot 
%Colors
color_cat = [0 0.45 0.75];
color_dis =  [0.9 0.1 0];
if with_cross_task
    color_cross_1 = [0.2 0.3 0.4];
    color_cross_2 = [0.7 0.2 0.7];
end
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
plot(corr_both_categorization,'LineWidth',2,'Color',color_cat);
hold on;
plot(corr_both_distraction,'LineWidth',2,'Color',color_dis);

if with_cross_task
    plot(corr_both_cross_task_categorization_dist,'LineWidth',2,'Color',color_cross_1)
    hold on;
    plot(corr_both_cross_task_distraction_dist,'LineWidth',2,'Color',color_cross_2)
end

%% Plot stats if needed
if with_stats
    if with_cross_task
        num_tasks = 4;
    else
        num_tasks = 2;
    end
    for task = 1:num_tasks
        if task<3
            task_name = get_task_name(task);
        else
            task_name = get_task_name(task-2); %for iteration 3: categ, iter 4: fix (distance task name)
        end
        
        if task == 1
            if strcmp(method,'fixed')
                plot_location = -0.7;
            elseif strcmp(method,'random')
                plot_location = -0.3;
            end
            color = color_cat;
        elseif task == 2
            if strcmp(method,'fixed')
                plot_location = -0.75;
            elseif strcmp(method,'random')
                plot_location = -0.32;
            end
            color = color_dis;
        elseif task == 3
            if strcmp(method,'fixed')
                plot_location = -0.8;
            elseif strcmp(method,'random')
                plot_location = -0.34;
            end
            color = color_cross_1;
        elseif task == 4
            if strcmp(method,'fixed')
                plot_location = -0.85;
            elseif strcmp(method,'random')
                plot_location = -0.36;
            end
            color = color_cross_2;
        end
   
        %Load and plot the stats
        if task<3
            filename_stats = 'dth_permutation_stats';
        elseif (task<=4)&&(task>2)
            filename_stats = 'dth_permutation_stats_crosstask';
        end
        if strcmp(method,'random')
            filename_stats = sprintf('random_%s',filename_stats);
        end
        filename_sign = fullfile(results_dir,sprintf('%s_both_%s_distance_subjects_%d_%d.mat',...
            filename_stats,task_name,subjects(1),subjects(end)));
        load(filename_sign);
        significant_timepoints = permutation_stats.SignificantMaxClusterWeight;

        %Plot the stats
        st = (significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;   
            
    end
end 

%Plotting parameters
plot_title =  sprintf('Correlation between the distance to hyperplane and reaction time in 60 scenes (N=%d)',numel(subjects));
if with_cross_task
    legend_plot = {'Scene categorization','Distraction',...
        'Cross-task: scene categorization distances with distraction RTs',...
        'Cross-task: distraction distances with categorization RTs'};
else
    legend_plot = {'Scene categorization','Distraction'};
end
xticks(0:10:200);
plotting_parameters(plot_title,legend_plot,40,12,'best','Spearman''s coefficient'); %[0.7 0.85 0.1 0.01]

%% Save
if with_cross_task
    saveas(gcf,fullfile(results_dir,sprintf('%s_subjects_%d_%d_both_tasks_with_crosstask',filename,subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_dir,sprintf('%s_subjects_%d_%d_both_tasks_with_crosstask.svg',filename,subjects(1),subjects(end))));
else
    saveas(gcf,fullfile(results_dir,sprintf('%s_subjects_%d_%d_both_tasks',filename,subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_dir,sprintf('%s_pseudotrials_SVM_DTH_subjects_%d_%d_both_tasks.svg',filename,subjects(1),subjects(end))));
end
close(gcf);
end