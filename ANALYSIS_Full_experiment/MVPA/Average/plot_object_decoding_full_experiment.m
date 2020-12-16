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
decoding_accuracies_all_subjects = NaN(numel(subjects),numConditions,numConditions,numTimepoints);

%% Set up the figure for plotting
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
sorted_subjects = sort(subjects); %order by ID
legend_cell = cell(1,sorted_subjects(end));
legend_cell(:) = {NaN};

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)));
    decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;
   
    %plot the object decoding curve for the participant    
    avg_over_conditions = squeeze(nanmean(nanmean(decodingAccuracy_avg,1),2));
    plot(avg_over_conditions, 'Linewidth',2);
    hold on;
    legend_cell{subject} = sprintf('Subject %d',subject);   
end   

%% Remove any NaN (for non-included subjects)
avg_over_conditions_all_subjects = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects,1),2),3));
legend_cell(cellfun(@(x) any(isnan(x)),legend_cell)) = [];

%% Plot the average of all subjects
plot(avg_over_conditions_all_subjects,'--','Color','k','Linewidth',3);
hold on;
title = sprintf('Object decoding per timepoint for %d subjects',numel(subjects));
onset_time = 40; 
legend_cell{numel(subjects)+1} = 'Average over all subjects';
legend_cell{numel(subjects)+2} = 'Stimulus onset'; %last legend element
plotting_parameters(title,legend_cell,onset_time,10);

%save the plot
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_object_decoding_%d_subjects_%s',numel(subjects),task_name))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('svm_object_decoding_%d_subjects_%s.svg',numel(subjects),task_name))); %save as svg
close(gcf);    

end