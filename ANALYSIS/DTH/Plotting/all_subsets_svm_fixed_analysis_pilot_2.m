function all_subsets_svm_fixed_analysis_pilot_2(subjects,subsets) %distance art, distance nat, RT
%ALL_SUBSETS_SVM_FIXED_ANALYSIS_PILOT_2 Performs the distance-to-hyperplane analysis using
%the svm classifier and 60 scenes. Combine the plots from different subsets
%of data.
%
%Input: subjects' ID (e.g., 1:13)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS';
addpath(genpath(results_dir));

%% Get the distances from all subjects
%Setup the subset names
subset_str = char(numel(subsets));
subset_str_all = char(numel(subsets));
for s = 1:numel(subsets)
    subset_str(s) = num2str(subsets(s));
    subset_str_all = cat(subset_str_all,num2str(subsets(s)));
end

%Preallocate
numTimepoints = 200;
numConditions = 60;
distances = NaN(numel(subjects),numConditions,numTimepoints);
RTs = NaN(numel(subjects),numConditions);
legend_cell = cell(numel(subsets),1);

%Setup figure
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen

%Subject loop
for subject = subjects
    subname = get_subject_name(subject);
    
    for s = 1:numel(subsets)
        load(fullfile(results_dir,subname,sprintf('subset_%s_decisionValues.mat',subset_str(s))));
        distances(subject,:,:) = decisionValues_Avg;   
        load(fullfile(results_dir,subname,sprintf('subset_%s_RTs_correct_trials.mat',subset_str(s))));
        RTs(subject,:) = normalize(RT_per_condition);
    end

    %% Get the median RTs of all subjects for each condition
    %Normalize and get median
    medianRT = nanmedian(RTs,1);

    %% Correlate DTH and RT
    mean_distances = squeeze(nanmean(distances,1)); %avg over subjects
    correlation_dth_rt_art = NaN(1,numTimepoints);
    last_category1_sample = 30;
    
    for t = 1:numTimepoints
        correlation_dth_rt_art(t) = corr(mean_distances(1:last_category1_sample,t),...
            medianRT(1:last_category1_sample)','type','Spearman');
    end

    %% Plot 
    %Artificial scenes
    plot(correlation_dth_rt_art,'LineWidth',2);
    hold on;
    legend_cell{s} = sprintf('%s_percent of the data',num2str(s*100));

end
    
%Plot parameters
title(sprintf('Correlation between the distance to hyperplane and reaction time for 30 artificial in a categorization task (N=%d)', numel(subjects)));
legend(legend_cell);
axis([0,200,-0.7,0.6]);

%% Save figure
save_path = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';
saveas(gcf,fullfile(save_path,sprintf('subset_%s_SVM_DTH_artificial_natural_%d_subjects',num2str(subset*100),numel(subjects)))); %turn subset into percentage: otherwise thinks it's a file extension
saveas(gcf,fullfile(save_path,sprintf('subset_%s_SVM_DTH_artificial_natural_%d_subjects.svg',num2str(subset*100),numel(subjects))));

end