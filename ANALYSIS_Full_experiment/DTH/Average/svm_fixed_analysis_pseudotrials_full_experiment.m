function svm_fixed_analysis_pseudotrials_full_experiment(subjects,task_distance,task_RT,with_stats) 
%SVM_FIXED_ANALYSIS_PSEUDOTRIALS_FULL_EXPERIMENT Performs the distance-to-hyperplane analysis using
%the svm classifier 60 scenes on normalized distances.
%
%Input: subjects' ID (e.g., 1:13), task of the EEG data (1=categorization, 2=distraction), 
%task of the RT data (1=categorization, 2=distraction), add stats (1 with, 0 without)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
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
    if size(decisionValues_Avg,1)==58
        keyboard;
    end
    distances(subject,:,:) = decisionValues_Avg;   
    
    %normalize RTs
    RTs(subject,:) = normalize(RT_per_condition);
    RTs_art(subject,:) = normalize(RT_per_condition(artificial_conditions));
    RTs_nat(subject,:) = normalize(RT_per_condition(natural_conditions));
end

%% Get the median RTs and mean distances of all subjects for each condition 
medianRT = nanmedian(RTs,1);
medianRT_art = nanmedian(RTs_art,1);
medianRT_nat = nanmedian(RTs_nat,1);
mean_distances = squeeze(nanmean(distances,1)); %avg over subjects

%% Correlate DTH and RT
t = 1:numTimepoints;
correlation_dth_rt_both = arrayfun(@(x) corr(mean_distances(:,x),medianRT','type','Spearman'),t);
correlation_dth_rt_art  = arrayfun(@(x) corr(mean_distances(artificial_conditions,x),medianRT_art','type','Spearman'),t);
correlation_dth_rt_nat  = arrayfun(@(x) corr(mean_distances(natural_conditions,x),medianRT_nat','type','Spearman'),t);
correlation_dth_rt_avg = mean([correlation_dth_rt_art;correlation_dth_rt_nat],1);

%% Plot 
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
color_art = [0 0.45 0.75];
color_nat =  [0.9 0.1 0];
color_both = [0 0 0];
color_avg = [0.7 0.2 0.7];

plot(correlation_dth_rt_art,'LineWidth',2,'Color',color_art);
hold on;
plot(correlation_dth_rt_nat,'LineWidth',2,'Color',color_nat);
hold on;
plot(correlation_dth_rt_both,'LineWidth',2,'Color',color_both);
hold on;
plot(correlation_dth_rt_avg,'LineWidth',2,'Color',color_avg);
hold on;

%% Plot stats if needed
if with_stats
    for c = 1:4
        if c == 1
            category = 'artificial';            
            plot_location = -0.75;
            color = color_art;
        elseif c == 2
            category = 'natural';            
            plot_location = -0.8;
            color = color_nat;
        elseif c == 3
            category = 'both';            
            plot_location = -0.85;
            color = color_both;
        elseif c == 4
            category = 'average';            
            plot_location = -0.9;
            color = color_avg;
        end

        significant_timepoints = weighted_cluster_perm_stats(subjects,task,category,0);
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

%% Save
%correlations
dth_results.corr_both_categories = correlation_dth_rt_both;
dth_results.corr_artificial = correlation_dth_rt_art;
dth_results.corr_natural = correlation_dth_rt_nat;
dth_results.corr_avg_categories = correlation_dth_rt_avg;
save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
if isequal(task_distance,task_RT)
    save(fullfile(save_path,sprintf('pseudotrials_SVM_DTH_%d_subjects_%s.mat',numel(subjects),task_name_distance)),'dth_results');
else
    save(fullfile(save_path,sprintf('pseudotrials_SVM_DTH_%d_subjects_cross_task_%s_distances.mat',numel(subjects),task_name_distance)),'dth_results');
end

%figures
saveas(gcf,fullfile(save_path,sprintf('pseudotrials_SVM_DTH_%d_subjects_%s',numel(subjects),task_name))); 
saveas(gcf,fullfile(save_path,sprintf('pseudotrials_SVM_DTH_%d_subjects_%s.svg',numel(subjects),task_name)));

end