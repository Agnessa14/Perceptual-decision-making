function distances_both_tasks_dth_pseudotrials_full_experiment(subjects,with_stats) 
%DISTANCES_BOTH_TASKS_DTH_PSEUDOTRIALS_FULL_EXPERIMENT Performs the distance-to-hyperplane analysis using
%the svm classifier 60 scenes on normalized distances, using the EEG data
%from both tasks.
%
%Input: subjects' ID (e.g., 1:13), add stats (1 with, 0 without)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%% Get the distances from all subjects
numTimepoints = 200;
numConditions = 60;
distances_categ = NaN(max(subjects),numConditions,numTimepoints);
distances_distr = NaN(max(subjects),numConditions,numTimepoints);
RTs = NaN(max(subjects),numConditions);
RTs_art = NaN(max(subjects),numConditions/2);
RTs_nat = NaN(size(RTs_art));
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;
for subject = subjects
    subname = get_subject_name(subject);
    
    %categorization
    load(fullfile(results_dir,subname,'dth_pseudotrials_svm_decisionValues_categorization.mat'));
    distances_categ(subject,:,:) = decisionValues_Avg;
    clear decisionValues_Avg;
    
    %distraction
    load(fullfile(results_dir,subname,'dth_pseudotrials_svm_decisionValues_fixation.mat'));
    distances_distr(subject,:,:) = decisionValues_Avg;
    clear decisionValues_Avg;
    
    %RTs
    load(fullfile(results_dir,subname,'RTs_correct_trials_categorization.mat'));
    
    %normalize RTs
    RTs(subject,:) = normalize(RT_per_condition);
    RTs_art(subject,:) = normalize(RT_per_condition(artificial_conditions));
    RTs_nat(subject,:) = normalize(RT_per_condition(natural_conditions));
end

%% Get the median RTs and mean distances of all subjects for each condition 
medianRT = nanmedian(RTs,1);
medianRT_art = nanmedian(RTs_art,1);
medianRT_nat = nanmedian(RTs_nat,1);

% mean_distances_categ = squeeze(nanmean(distances_categ,1)); %avg over subjects
% mean_distances_distr = squeeze(nanmean(distances_distr,1)); 

mean_distances = nanmean(cat(4,distances_categ,distances_distr),4);

%% Correlate DTH and RT
% t = 1:numTimepoints;
% correlation_dth_rt_both = arrayfun(@(x) corr(mean_distances(:,x),medianRT','type','Spearman'),t);
% correlation_dth_rt_art  = arrayfun(@(x) corr(mean_distances(artificial_conditions,x),medianRT_art','type','Spearman'),t);
% correlation_dth_rt_nat  = arrayfun(@(x) corr(mean_distances(natural_conditions,x),medianRT_nat','type','Spearman'),t);
% correlation_dth_rt_avg  = mean([correlation_dth_rt_art;correlation_dth_rt_nat],1);

t = 1:numTimepoints;
size_corr = [max(subjects),numTimepoints];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);
correlation_avg = NaN(size_corr);

for subject = subjects
    correlation_art(subject,:) = arrayfun(@(x) corr(squeeze(mean_distances(subject,artificial_conditions,x))',medianRT_art','type','Spearman'), t);
    correlation_nat(subject,:) = arrayfun(@(x) corr(squeeze(mean_distances(subject,natural_conditions,x))',medianRT_nat','type','Spearman'), t);
    correlation_both(subject,:) = arrayfun(@(x) corr(squeeze(mean_distances(subject,:,x))',medianRT','type','Spearman'), t);
    correlation_avg(subject,:) = mean([squeeze(correlation_art(subject,:));squeeze(correlation_nat(subject,:))],1);
end

%% Average over participants
avg_corr_art = squeeze(nanmean(correlation_art,1));
avg_corr_nat = squeeze(nanmean(correlation_nat,1));
avg_corr_both = squeeze(nanmean(correlation_both,1));
avg_corr_avg = squeeze(nanmean(correlation_avg,1));

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
color_art = [0 0.45 0.75];
color_nat =  [0.9 0.1 0];
color_both = [0 0 0];
color_avg = [0.7 0.2 0.7];

plot(avg_corr_art,'LineWidth',2,'Color',color_art);
hold on;
plot(avg_corr_nat,'LineWidth',2,'Color',color_nat);
plot(avg_corr_both,'LineWidth',2,'Color',color_both);
plot(avg_corr_avg,'LineWidth',2,'Color',color_avg);

%% Plot stats if needed
if with_stats
    task_distance = 4;
    task_RT = 1;
    num_perms = 1000;
    for c = 1:4
        %Setup some plot variables
        if c == 1
            category = 'artificial';            
            plot_location = -0.75;
            color = color_art;
            medianRT_stats = medianRT_art;
            distances_stats = mean_distances(:,artificial_conditions,:);  
            true_correlation = avg_corr_art;
        elseif c == 2
            category = 'natural';            
            plot_location = -0.8;
            color = color_nat;
            medianRT_stats = medianRT_nat;
            distances_stats = mean_distances(:,natural_conditions,:);    
            true_correlation = avg_corr_nat;
        elseif c == 3
            category = 'both';            
            plot_location = -0.85;
            color = color_both;
            medianRT_stats = medianRT;      
            true_correlation = avg_corr_both;
            distances_stats = mean_distances;
        elseif c == 4
            category = 'average';            
            plot_location = -0.9;
            color = color_avg;
            medianRT_stats = medianRT;
            true_correlation = avg_corr_avg;
            distances_stats = mean_distances;
        end
        
        %Check if stats already exist, otherwise run the stats script
        filename = 'random_both_tasks_distances_dth_permutation_stats';
        filename_sign = fullfile(results_avg_dir,sprintf('%s_%s_distance_subjects_%d_%d.mat',...
            filename,category,subjects(1),subjects(end)));
        if exist(filename_sign,'file')
            load(filename_sign);
            significant_timepoints = permutation_stats.SignificantMaxClusterWeight;
        else
            significant_timepoints = weighted_cluster_perm_stats(subjects,medianRT_stats,distances_stats,true_correlation,task_distance,task_RT,category,'left',0,num_perms,'random');
%             significant_timepoints = distances_both_tasks_weighted_cluster_perm_stats(subjects,category,1,num_perms,'fixed');
        end
        
        %Plot the stats
        st = (significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;        
    end
end 

%Plotting parameters
plot_title =  sprintf('Correlation between the distance to hyperplane from the categorization and distraction tasks and reaction time from the categorization task (N=%d)',numel(subjects));
legend_plot = {'Artificial scenes','Natural scenes','All scenes',...
    'Average of artificial and natural scenes'}; %add stimulus onset?
xticks(0:10:200);
plotting_parameters(plot_title,legend_plot,40,12,'best','Spearman''s coefficient'); %[0.4 0.8 0.1 0.1

%% Save
%correlations
dth_results.corr_both_categories = correlation_dth_rt_both;
dth_results.corr_artificial = correlation_dth_rt_art;
dth_results.corr_natural = correlation_dth_rt_nat;
dth_results.corr_avg_categories = correlation_dth_rt_avg;
save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
save(fullfile(save_path,sprintf('random_distances_both_tasks_SVM_DTH_subjects_%d_%d_distances.mat',subjects(1),subjects(end))),'dth_results');
saveas(gcf,fullfile(save_path,sprintf('random_distances_both_tasks_SVM_DTH_subjects_%d_%d_distances',subjects(1),subjects(end))));
saveas(gcf,fullfile(save_path,sprintf('random_distances_both_tasks_SVM_DTH_subjects_%d_%d_distances.svg',subjects(1),subjects(end))));

close(gcf);
end

