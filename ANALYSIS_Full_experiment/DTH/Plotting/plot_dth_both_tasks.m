function plot_dth_both_tasks(subjects,with_cross_task,with_stats,with_error_bars,varargin) %distance art, distance nat, RT
%PLOT_DTH_BOTH_TASKS Plot the fixed effects distance-to-hyperplane results for both
%tasks. Assumes that the correlations & stats have already been computed.
%
%Input: subjects' ID (e.g., 1:13), add the cross task curves (1 yes, 0 no),
%add stats (1 yes, 0 no), plot with/without error bars (1/0), varargin:
%stats_type ('cluster' or 'fdr')
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%Author: Agnessa Karapetian,2021
%

%% Load data
%Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
results_dir_sub = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';

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

%% If with error bars, calculate the correlation for each subject
%Load subject-specific data
if with_error_bars
    numTimepoints = 200;
    numConditions = 60;
    distances_categorization = NaN(max(subjects),numConditions,numTimepoints);
    RTs_categorization = NaN(max(subjects),numConditions);
    distances_distraction = NaN(max(subjects),numConditions,numTimepoints);
    RTs_distraction = NaN(max(subjects),numConditions);
    
    for subject = subjects
        subname = get_subject_name(subject);
        %categorization
        load(fullfile(results_dir_sub,subname,...
            'cross_validated_dth_pseudotrials_svm_decisionValues_categorization.mat'),'decisionValues_Avg');
        load(fullfile(results_dir_sub,subname,'RTs_correct_trials_categorization.mat'),'RT_per_condition');       
        distances_categorization(subject,:,:) = decisionValues_Avg;   
        RTs_categorization(subject,:) = RT_per_condition; 
        
        %distraction
        load(fullfile(results_dir_sub,subname,...
            'cross_validated_dth_pseudotrials_svm_decisionValues_fixation.mat'),'decisionValues_Avg');
        load(fullfile(results_dir_sub,subname,'RTs_correct_trials_fixation.mat'),'RT_per_condition');
        distances_distraction(subject,:,:) = decisionValues_Avg;   
        RTs_distraction(subject,:) = RT_per_condition; 
    end

    %Get the median RTs of all subjects for each scene
    medianRT_categorization = nanmedian(RTs_categorization,1);
    medianRT_distraction = nanmedian(RTs_distraction,1);
    
    %Correlate each subject's distances with the median RT
    t = 1:numTimepoints;
    size_corr = [max(subjects),numTimepoints];
    corr_cat_all_subjects = NaN(size_corr);
    corr_dis_all_subjects = NaN(size_corr);
    corr_crosstask_dist_cat_all_subjects = NaN(size_corr);
    corr_crosstask_dist_dis_all_subjects = NaN(size_corr);
    for subject = subjects
        distances_subject_cat = squeeze(distances_categorization(subject,:,:));
        distances_subject_dis = squeeze(distances_distraction(subject,:,:));
        
        %Find any missing scenes (only in categorization) and exclude them
        conds = 1:numConditions;
        conditions_subject = ~isnan(distances_subject_cat(:,1));

        %Correlate RT and distances    
        corr_cat_all_subjects(subject,:) = arrayfun(@(x) ...
            corr(squeeze(distances_subject_cat(conditions_subject,x)),...
            medianRT_categorization(conditions_subject)','type','Spearman'), t);
        corr_dis_all_subjects(subject,:) = arrayfun(@(x) ...
            corr(squeeze(distances_subject_dis(conds,x)),...
            medianRT_distraction(conds)','type','Spearman'), t);
        if with_cross_task
            corr_crosstask_dist_cat_all_subjects(subject,:) = arrayfun(@(x) ...
                corr(squeeze(distances_subject_cat(conditions_subject,x)),...
                medianRT_distraction(conditions_subject)','type','Spearman'), t);
            corr_crosstask_dist_dis_all_subjects(subject,:) = arrayfun(@(x) ...
                corr(squeeze(distances_subject_dis(conditions_subject,x)),...
                medianRT_categorization(conditions_subject)','type','Spearman'), t);
        end
    end
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
% categ = plot(corr_both_categorization,'LineWidth',2,'Color',color_cat);
% hold on;
% distr = plot(corr_both_distraction,'LineWidth',2,'Color',color_dis);
% 
% if with_cross_task
%     categ_eeg = plot(corr_both_cross_task_categorization_dist,'LineWidth',2,'Color',color_cross_1);
%     hold on;
%     distr_eeg = plot(corr_both_cross_task_distraction_dist,'LineWidth',2,'Color',color_cross_2);
% end

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
            data = corr_both_categorization;
            if with_error_bars
                for_stats = corr_cat_all_subjects(subjects,:);
            end
        elseif task == 2
            plot_location = 0.18;
            color = color_dis;
            data = corr_both_distraction;
            if with_error_bars
                for_stats = corr_dis_all_subjects(subjects,:);
            end
        elseif task == 3
            plot_location = 0.19;
            color = color_cross_1;
            data = corr_both_cross_task_categorization_dist;
            if with_error_bars
                for_stats = corr_crosstask_dist_cat_all_subjects(subjects,:);
            end
        elseif task == 4
            plot_location = 0.17;
            color = color_cross_2;
            data = corr_both_cross_task_distraction_dist;
            if with_error_bars
                for_stats = corr_crosstask_dist_dis_all_subjects(subjects,:);
            end
        end
   
        %Load and plot the stats
        if task<3
            filename_sign = 'dth_permutation_stats';
        elseif (task<=4)&&(task>2)
            filename_sign = 'dth_permutation_stats_crosstask';
        end               
        if isempty(varargin)
            error('Specify the type of stats')
        else
            stats_type = varargin{1};
        end
        filename = fullfile(results_dir,sprintf('%s_%d_%d_%s_task_%s_both.mat',filename_sign,...
            subjects(1),subjects(end),task_name,analysis));   
        
        if strcmp(stats_type,'cluster')
            load(filename,'permutation_stats'); 
            st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
        elseif strcmp(stats_type,'fdr')
            load(filename,'fdr_stats'); 
            st = (fdr_stats.significant_timepoints*plot_location); %depending on the stats            
        end
            
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;     
        
        %if needed: plot error bars
        if with_error_bars
            %calculate error bars
            stdDM = std(for_stats); 
            err = stdDM/sqrt(size(for_stats,1)); %standard deviation/sqrt of num subjects  

            %plot as a shaded area
            top_curve = data + err;
            bottom_curve = data - err;
            x2 = [1:numTimepoints, fliplr(1:numTimepoints)];
            shaded_area = [top_curve, fliplr(bottom_curve)];
            fill(x2, shaded_area, color,'FaceAlpha',0.5);
            hold on;
        end

        %Plot average curve
        plot(data,'LineWidth',2,'Color',color);
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
file_name = 'cv_all_dist_med_rt_dth_both_tasks';
if with_error_bars
    file_name = sprintf('error_bars_%s',file_name);
end

if ~isempty(varargin)
    if strcmp(stats_type,'cluster')
        file_name = sprintf('%s_cluster',file_name);
    elseif strcmp(stats_type,'fdr')
        file_name = sprintf('%s_fdr',file_name);
    end
end

if with_cross_task
    saveas(gcf,fullfile(results_dir,sprintf('%s_%d_%d_with_crosstask',file_name,subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_dir,sprintf('%s_subjects_%d_%d_with_crosstask.svg',file_name,subjects(1),subjects(end))));
else
    saveas(gcf,fullfile(results_dir,sprintf('%s_subjects_%d_%d',file_name,subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_dir,sprintf('%s_subjects_%d_%d.svg',file_name,subjects(1),subjects(end))));
end
close(gcf);
end