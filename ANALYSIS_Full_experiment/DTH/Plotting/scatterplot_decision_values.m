function scatterplot_decision_values(subjects,task) 
%SCATTERPLOT_DECISION_VALUES Creates a scatterplot of decision values for
%each timepoint.
%
%Input: subjects' ID (e.g., 1:13), task (1=categorization, 2=distraction)
%
%Output: P plots of decision values per scene, where P is the number of
%timepoints.
%
%
%Author: Agnessa Karapetian, 2021
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
addpath(genpath(results_dir));
task_name = get_task_name(task);

%% Get the distances from all subjects
numTimepoints = 200;
numConditions = 60;
decision_values = NaN(max(subjects),numConditions,numTimepoints);
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('classification_results_dth_%s.mat',task_name)));
    decision_values(subject,:,:) = classification_results.decisionValues;      
end

mean_dvs = squeeze(nanmean(decision_values,1)); %avg over subjects

%Plot
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
title_part = 'Decision values for scenes at';

for t = 1:numTimepoints
    scatter(1:numConditions,mean_dvs(:,t),'filled');
    xlabel('Scene ID');
    ylabel('Decision value (arbitrary unit)');
    title(sprintf('%s %d ms',title_part,(t-40)*5));
    saveas(gcf,fullfile(results_dir,subname,sprintf('scatterplot_decision_values_timepoint_%d',t)));
end

keyboard;

end

