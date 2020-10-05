function plot_object_decoding(subjects)
%PLOT_OBJECT_DECODING Plot the results from object decoding, averaged over
%all participants and subject-specific.
%
%Input: subject IDs
%
%Output: curve of decoding accuracies per timepoint

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS/'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS/';

%% Preallocate
numConditions = 60;
numTimepoints = 200;
decoding_accuracies_all_subjects = NaN(numel(subjects),numConditions,numConditions,numTimepoints);

%% Set up the figure for plotting
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
legend_cell = cell(1,4);

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy.mat')));
    
    %matrix of all subjects
    decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;
   
    %plot the object decoding curve for the participant
    avg_over_conditions = squeeze(nanmean(decodingAccuracy_avg,[1,2]));
    plot(avg_over_conditions, 'Linewidth',2);
    hold on;
    
    %legend
    legend_cell{subject} = sprintf('Subject %d',subject);
    
end   
title = sprintf('Classification accuracy of individual scenes per timepoint in a categorization task (N=%d)',...
        numel(subjects));
onset_time = 40; 
legend_cell{5} = 'Stimulus onset';
plotting_parameters(title,legend_cell,onset_time);

%save the plot
saveas(gcf,fullfile(subject_results_dir,'svm_object_decoding')); %save as matlab figure
saveas(gcf,fullfile(subject_results_dir,'svm_object_decoding.svg')); %save as svg
close(gcf);    

%% Plot the average of all subjects
avg_over_conditions_all_subjects = squeeze(nanmean(decoding_accuracies_all_subjects,1:3));
plot(avg_over_conditions_all_subjects,'Linewidth',2) 

end