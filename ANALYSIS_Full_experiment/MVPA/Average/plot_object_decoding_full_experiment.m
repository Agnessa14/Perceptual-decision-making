function plot_object_decoding_full_experiment(subjects,task)
%PLOT_OBJECT_DECODING_FULL_EXPERIMENT Plot the results from object decoding, averaged over
%all participants and subject-specific.
%
%Input: subject IDs, task (1=categorization,2=distraction)
%
%Output: curve of decoding accuracies per timepoint

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
task_name = get_task_name(task);

%% Preallocate
numConditions = 60;
numTimepoints = 200;

%% Set up the figure for plotting
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
sorted_subjects = sort(subjects); %order by ID
decoding_accuracies_all_subjects = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
legend_cell = cell(1,sorted_subjects(end));
legend_cell(:) = {NaN};

cmap = jet(max(subjects));
%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)));
    decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;
   
    %plot the object decoding curve for the participant    
    avg_over_conditions = squeeze(nanmean(nanmean(decodingAccuracy_avg,1),2));
    plot(avg_over_conditions, 'Linewidth',2, 'Color', cmap(subject, :));
    hold on;
    legend_cell{subject} = sprintf('Subject %d',subject);   
end   

%% Remove any NaN (for non-included subjects)
avg_over_conditions_all_subjects = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects,1),2),3));
legend_cell(cellfun(@(x) any(isnan(x)),legend_cell)) = [];

%% Plot the average of all subjects
plot(avg_over_conditions_all_subjects,'--','Color','k','Linewidth',3);
hold on;
if task==1
    task_title = 'scene categorization';
elseif task==2
    task_title = 'distraction';
end
title = sprintf('Object decoding per timepoint for %d subjects in a %s task',numel(subjects),task_title);
onset_time = 40; 
legend_cell{numel(subjects)+1} = 'Average over all subjects';
legend_cell{numel(subjects)+2} = 'Stimulus onset'; %last legend element
xticks(0:10:200);
plotting_parameters(title,legend_cell,onset_time,8,'best','Decoding accuracy (%)'); %'best' location or specify position [0.75 0.7 0.1 0.1]

%save the plot
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_object_decoding_%d_subjects_%s',numel(subjects),task_name))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_object_decoding_%d_subjects_%s.svg',numel(subjects),task_name))); %save as svg
close(gcf);    

end