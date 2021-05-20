function dth_all_distances_median_RT(subjects,task_distance,task_RT,with_stats) %distance art, distance nat, RT
%DTH_ALL_DISTANCES_MEDIAN_RT Performs the distance-to-hyperplane analysis using
%the distances-to-hyperplane for each subject and the median RT across subjects.
%
%Input: subjects' ID (e.g., 1:13), task for the distances (1=categorization, 2=distraction), task for RT,
%add stats (1 with, 0 without)
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
    load(fullfile(results_dir,subname,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_name_distance)));
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name_RT)));
    distances(subject,:,:) = decisionValues_Avg;   
    RTs(subject,:) = normalize(RT_per_condition);
    RTs_art(subject,:) = normalize(RT_per_condition(artificial_conditions));
    RTs_nat(subject,:) = normalize(RT_per_condition(natural_conditions));
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
correlation_avg = NaN(size_corr);

for subject = subjects
    correlation_art(subject,:) = arrayfun(@(x) corr(squeeze(distances(subject,artificial_conditions,x))',medianRT_art','type','Spearman'), t);
    correlation_nat(subject,:) = arrayfun(@(x) corr(squeeze(distances(subject,natural_conditions,x))',medianRT_nat','type','Spearman'), t);
    correlation_both(subject,:) = arrayfun(@(x) corr(squeeze(distances(subject,:,x))',medianRT','type','Spearman'), t);
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
    num_perms = 1000;
    for c = 1:4
        if c == 1
            category = 'artificial';            
            plot_location = -0.34;
            color = color_art;
            medianRT_stats = medianRT_art;
            distances_stats = distances(:,artificial_conditions,:);  
            true_correlation = avg_corr_art;
        elseif c == 2
            category = 'natural';            
            plot_location = -0.37;
            color = color_nat;
            medianRT_stats = medianRT_nat;
            distances_stats = distances(:,natural_conditions,:);    
            true_correlation = avg_corr_nat;
        elseif c == 3
            category = 'both';            
            plot_location = -0.4;
            color = color_both;
            medianRT_stats = medianRT;      
            true_correlation = avg_corr_both;
            distances_stats = distances;
        elseif c == 4
            category = 'average';            
            plot_location = -0.43;
            color = color_avg;
            medianRT_stats = medianRT;
            true_correlation = avg_corr_avg;
            distances_stats = distances;
        end

        %Check if stats already exist, otherwise run the stats script
        filename_sign = fullfile(results_avg_dir,sprintf('random_dth_permutation_stats_%s_%s_distance_subjects_%d_%d.mat',...
            category,task_name_distance,subjects(1),subjects(end)));
        if exist(filename_sign,'file')
            load(filename_sign);
            significant_timepoints = permutation_stats.SignificantMaxClusterWeight;
        else
            significant_timepoints = weighted_cluster_perm_stats(subjects,medianRT_stats,distances_stats,true_correlation,task_distance,task_RT,category,'left',1,num_perms,'random');
        end
        
        %Plot
        st = (significant_timepoints*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;        
    end
end 

%Plotting parameters
if isequal(task_distance,task_RT)
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
legend_plot = {'Artificial scenes','Natural scenes','All scenes',...
    'Average of artificial and natural scenes'}; %add stimulus onset?
xticks(0:10:200);
plotting_parameters(plot_title,legend_plot,40,12,'best','Spearman''s coefficient'); %[0.4 0.8 0.1 0.1
ylim([-0.5 0.4]);
%% Save figures
save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

if isequal(task_distance,task_RT)
    saveas(gcf,fullfile(save_path,sprintf('random_effects_dth_subjects_%d_%d_%s',subjects(1),subjects(end),task_name_distance))); 
    saveas(gcf,fullfile(save_path,sprintf('random_effects_dth_subjects_%d_%d_%s.svg',subjects(1),subjects(end),task_name_distance))); 
else
    saveas(gcf,fullfile(save_path,sprintf('random_effects_dth_subjects_%d_%d_cross_task_%s_distances',subjects(1),subjects(end),task_name_distance)));
    saveas(gcf,fullfile(save_path,sprintf('random_effects_dth_subjects_%d_%d_cross_task_%s_distances.svg',subjects(1),subjects(end),task_name_distance)));
end

close(gcf);
end





