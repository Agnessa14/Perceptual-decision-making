function encoding_all_subjects(subjects,task) %distance art, distance nat, RT
%ENCODING_ALL_SUBJECTS Average encoding accuracy across subjects.
%
%Input: subjects' ID (e.g., 1:13), task(1=categorization, 2=distraction) 
%
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));
task_name = get_task_name(task);


%% Get the distances from all subjects
numTimepoints = 200;
encoding_accuracies = NaN(max(subjects),numTimepoints);
filename = sprintf('cross_validated_regression_RTs_encodingAccuracy_%s',task_name);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,filename),'encodingAccuracy_avg');
    encoding_accuracies(subject,:) = encodingAccuracy_avg;   
end

%% Average over participants
encoding_accuracies_subjects = encoding_accuracies(subjects,:);
avg_encoding_acc = squeeze(mean(encoding_accuracies_subjects,1));

%% Plot
plot(avg_encoding_acc,'LineWidth',2);
plot_title = 'Regression of RTs on EEG patterns';
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
legend_plot = {'Artificial scenes','Natural scenes'}; 
legend_bool = 0;
title_bool = 1;
plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Encoding accuracy'); 

%% Save
filename_save = sprintf('average_%s',filename);
save(fullfile(results_avg_dir,sprintf('subject_level_%s',filename)),'encoding_accuracies_subjects');
save(fullfile(results_avg_dir,filename_save),'avg_encoding_acc');
saveas(gcf,fullfile(results_avg_dir,filename_save));
saveas(gcf,fullfile(results_avg_dir,sprintf('%s.svg',filename_save)));

close(gcf);

