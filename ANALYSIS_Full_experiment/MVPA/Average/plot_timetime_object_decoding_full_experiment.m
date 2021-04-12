function plot_timetime_object_decoding_full_experiment(subjects,task)
%PLOT_TIMETIME_OBJECT_DECODING_FULL_EXPERIMENT Plot the results from time-generalized object decoding, averaged over
%all participants and subject-specific.
%
%Input: subject IDs, task (1=categorization,2=distraction)
%
%Output: curve of decoding accuracies per timepoint

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate - too large matrix, takes too long
% numConditions = 60;
% numTimepoints = 50;

% sorted_subjects = sort(subjects); %order by ID
  
% legend_cell = cell(1,sorted_subjects(end));
% legend_cell(:) = {NaN};

% cmap = jet(max(subjects));

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'time_generalized_svm_decoding_accuracy_crosstask.mat'));
    decoding_accuracies_all_subjects(subject,:,:,:,:) = timeg_decodingAccuracy_avg;
    
%     %plot the object decoding curve for the participant
%     decodingAccuracy_avg = squeeze(nanmean(nanmean(decodingAccuracy_avg,1),2));   %avg over conditions
%     plot(decodingAccuracy_avg, 'Linewidth',2, 'Color', cmap(subject, :));
%     hold on;
%     legend_cell{subject} = sprintf('Subject %d',subject);   
end   

%% Average over subjects and conditions and remove any NaN (for non-included subjects and conditions)
avg_over_conditions_all_subjects = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects,1),2),3)); 
% legend_cell(cellfun(@(x) any(isnan(x)),legend_cell)) = [];

%% Plot the average of all subjects
hold on;
if task==1
    task_title = 'scene categorization';
elseif task==2
    task_title = 'distraction';
end

plot_title = sprintf('Time-generalized object decoding for %d subjects in a %s task',numel(subjects),task_title);
title(plot_title);
ylabel('Trained on timepoints: Categorization');
xlabel('Tested on timepoints: Distraction');
% legend_cell{numel(subjects)+1} = 'Average over all subjects';
% legend_cell{numel(subjects)+2} = 'Stimulus onset'; %last legend element
xticks(0:10:200);


% Set up the figure for plotting
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
h = pcolor(avg_over_conditions_all_subjects); 
set(h, 'EdgeColor', 'none');
axis square;
% plotting_parameters(plot_title,legend_cell,onset_time,8,'best','Decoding accuracy (%)'); %'best' location or specify position [0.75 0.7 0.1 0.1]  

%% Save the plot
saveas(gcf,fullfile(results_avg_dir,sprintf('timegen_svm_object_decoding_%d_subjects_crosstask',numel(subjects)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('timegen_svm_object_decoding_%d_subjects_crosstask.svg',numel(subjects)))); %save as svg
close(gcf)
end