function dth_all_distances_median_RT(subjects,task_distance,task_RT,with_stats,with_error_bars) %distance art, distance nat, RT
%DTH_ALL_DISTANCES_MEDIAN_RT Performs the distance-to-hyperplane analysis using
%the distances-to-hyperplane for each subject and the median RT across subjects.
%
%Input: subjects' ID (e.g., 1:13), task for the distances (1=categorization, 2=distraction), task for RT,
%add stats (1 with, 0 without), plot with/without error bars (1/0)
%
%Correlates the decision values of each subject with reaction times averaged over
%participants of each condition (60 scenes), averaged, at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));
task_name_distance = get_task_name(task_distance);
task_name_RT = get_task_name(task_RT);

%% Get the distances from all subjects
numTimepoints = 200;
numConditions = 60;
distances = NaN(max(subjects),numConditions,numTimepoints);
RTs = NaN(max(subjects),numConditions);
RTs_art = NaN(max(subjects),numConditions/2);
RTs_nat = NaN(size(RTs_art));
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,...
        sprintf('cross_validated_dth_pseudotrials_svm_decisionValues_%s.mat',...
        task_name_distance)),'decisionValues_Avg');
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',...
        task_name_RT)),'RT_per_condition');
    distances(subject,:,:) = decisionValues_Avg;   
    RTs(subject,:) = RT_per_condition; 
    RTs_art(subject,:) = RT_per_condition(artificial_conditions);
    RTs_nat(subject,:) = RT_per_condition(natural_conditions);
end

%% Get the median RTs of all subjects for each condition 
medianRT = nanmedian(RTs,1);
medianRT_art = nanmedian(RTs_art,1);
medianRT_nat = nanmedian(RTs_nat,1);

%% Correlate each subject's distances with the median RT
t = 1:numTimepoints;
size_corr = [max(subjects),numTimepoints];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);

for subject = subjects
    %Find any missing scenes and exclude them
    distances_subject = squeeze(distances(subject,:,:));
    if any(ismember(artificial_conditions,find(isnan(distances_subject))))
        artificial_conditions_subject = artificial_conditions(~isnan(distances_subject(artificial_conditions,1)));
        natural_conditions_subject = natural_conditions;
    elseif any(ismember(natural_conditions,find(isnan(distances_subject))))
        artificial_conditions_subject = artificial_conditions;
        natural_conditions_subject = natural_conditions(~isnan(distances_subject(natural_conditions,1)));
    else
        artificial_conditions_subject = artificial_conditions;
        natural_conditions_subject = natural_conditions;
    end
    
    %Correlate RT and distances    
    all_conditions_subject = [artificial_conditions_subject natural_conditions_subject];
    correlation_art(subject,:) = arrayfun(@(x) ...
        corr(squeeze(distances_subject(artificial_conditions_subject,x)),...
        medianRT_art(artificial_conditions_subject)','type','Spearman'), t);
    correlation_nat(subject,:) = arrayfun(@(x) ...
        corr(squeeze(distances_subject(natural_conditions_subject,x)),...
        medianRT_nat(natural_conditions_subject-30)','type','Spearman'), t);
    correlation_both(subject,:) = arrayfun(@(x) ...
        corr(squeeze(distances_subject(all_conditions_subject,x)),...
        medianRT(all_conditions_subject)','type','Spearman'), t);
end

%% Average over participants
avg_corr_art = squeeze(nanmean(correlation_art,1));
avg_corr_nat = squeeze(nanmean(correlation_nat,1));
avg_corr_both = squeeze(nanmean(correlation_both,1));

%% Plot
figure(abs(round(randn*10)));
cmap_1 = cool;
cmap_2 = summer;
color_art = cmap_1(200,:); %purple
color_nat = cmap_2(100,:); %green

%% Plot stats if needed
for c = 1:3
    if c == 1
        data = avg_corr_art;
        color = color_art;
    elseif c == 2
        data = avg_corr_nat;
        color = color_nat;
    elseif c == 3
        data = avg_corr_both;
        color = 'k';
    end
    if ~with_stats
        %Plot average curve
        plot(data,'LineWidth',2,'Color',color);
        hold on;
    else
        %define some variables for the stats and the plot
        analysis = 'random_dth';
        permutation_stats.num_perms = 1000;
        permutation_stats.cluster_th = 0.05;
        permutation_stats.significance_th = 0.05;
        permutation_stats.tail = 'left';

        if c == 1
            category = 'artificial';            
            plot_location = 0.15;
            for_stats = correlation_art(subjects,:);
            if with_error_bars
                for_stats = correlation_art(subjects,:);
            end
        elseif c == 2
            category = 'natural';            
            plot_location = 0.17;
            for_stats = correlation_nat(subjects,:);
            if with_error_bars
                for_stats = correlation_nat(subjects,:);
            end
        elseif c == 3
            category = 'both';            
            for_stats = correlation_both(subjects,:);       
        end
        
        %Check if stats already exist, otherwise run the stats script
        if task_distance==task_RT
            filename_sign = 'dth_permutation_stats';
        else
            filename_sign = 'dth_permutation_stats_crosstask';
        end     
        filename = fullfile(results_avg_dir,sprintf('%s_%d_%d_%s_task_%s_%s.mat',filename_sign,...
            subjects(1),subjects(end),task_name_distance,analysis,category));
        if exist(filename,'file')
            load(filename,'permutation_stats');
        else
            [permutation_stats.SignificantMaxClusterWeight,permutation_stats.pValWeight,...
                permutation_stats.SignificantMaxClusterSize,permutation_stats.pValSize] = ...
                permutation_cluster_1sample_weight_alld(for_stats,permutation_stats.num_perms,...
                permutation_stats.cluster_th,permutation_stats.significance_th,permutation_stats.tail); 
        save(filename,'permutation_stats');
        end
        
        %Plot statistics
        if c < 3
            st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
            st(st==0) = NaN;
            plot(st,'*','Color',color); 
            hold on;
        end
        
        %if needed: plot error bars
        if c < 3
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
        end
    end
end 
    
%% Plotting parameters
if task_distance==task_RT
    if task_distance==1
        task_title = 'scene categorization';
    elseif task_distance==2
        task_title = 'distraction';
    end
    plot_title =  sprintf('Correlation between the distance to hyperplane and reaction time in a %s task (N=%d)',task_title,numel(subjects));
elseif task_distance==1 && task_RT==2
    plot_title =  sprintf('Correlation between the distance to hyperplane from a categorization task and reaction time from a distraction task (N=%d)',numel(subjects));
elseif task_distance==2 && task_RT==1
    plot_title =  sprintf('Correlation between the distance to hyperplane from a distraction task and reaction time from a categorization task (N=%d)',numel(subjects));
end
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
legend_plot = {'Artificial scenes','Natural scenes'}; 
ylim([-0.3 0.3]);
yticks(-0.3:0.1:0.3);
legend_bool = 0;
title_bool = 0;
plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Spearman''s r'); 

%% Save correlations and figures
dth_results.corr_both_categories = avg_corr_both;
dth_results.corr_artificial = avg_corr_art;
dth_results.corr_natural = avg_corr_nat;

save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
file_name = 'cv_all_dist_med_rt_dth';
if with_error_bars
    file_name = sprintf('error_bars_%s',file_name);
end
if isequal(task_distance,task_RT)
    save(fullfile(save_path,sprintf('%s_subjects_%d_%d_%s.mat',file_name,subjects(1),subjects(end),task_name_distance)),'dth_results');
    saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_%s',file_name,subjects(1),subjects(end),task_name_distance))); 
    saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_%s.svg',file_name,subjects(1),subjects(end),task_name_distance))); 
else
    save(fullfile(save_path,sprintf('%s_subjects_%d_%d_cross_task_%s_distances.mat',file_name,subjects(1),subjects(end),task_name_distance)),'dth_results');
    saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_cross_task_%s_distances',file_name,subjects(1),subjects(end),task_name_distance)));
    saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_cross_task_%s_distances.svg',file_name,subjects(1),subjects(end),task_name_distance)));
end

close(gcf);
end





