function plot_decoding_full_experiment(subjects,task,analysis)
%PLOT_DECODING_FULL_EXPERIMENT Plot the results from decoding, averaged over
%all participants and subject-specific.
%
%Input: subject IDs, task (1=categorization,2=distraction), analysis
%('object_decoding' or 'category_decoding')
%
%Output: curve of decoding accuracies per timepoint
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
task_name = get_task_name(task);

%% Preallocate
numConditions = 60;
numTimepoints = 200;

% Set up the figure for plotting
figure(abs(round(randn*10))); %Random figure number
%for full screen: % set(gcf, 'Position', get(0, 'Screensize'));
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects = NaN(max(subjects),numTimepoints);
end
   
legend_cell = cell(1,max(subjects));
legend_cell(:) = {NaN};

cmap = jet(max(subjects));

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)),...
            'decodingAccuracy_avg');
        decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;
    elseif strcmp(analysis,'category_decoding')
        load(fullfile(subject_results_dir,...
            sprintf('cross_validated_dth_pseudotrials_svm_decodingAccuracy_%s.mat',task_name)),...
            'decodingAccuracy_avg');
        decoding_accuracies_all_subjects(subject,:) = decodingAccuracy_avg;
    end
    
    %plot the object decoding curve for the participant
    if strcmp(analysis,'object_decoding')
        decodingAccuracy_avg = squeeze(nanmean(nanmean(decodingAccuracy_avg,1),2)); %avg over conditions
    end
    plot(decodingAccuracy_avg-50, 'Linewidth',2, 'Color', cmap(subject, :));
    hold on;
    legend_cell{subject} = sprintf('Subject %d',subject);   
end   

%% Remove any NaN (for non-included subjects and conditions)
if strcmp(analysis,'object_decoding')
    avg_over_conditions_all_subjects = squeeze(nanmean(nanmean(nanmean(decoding_accuracies_all_subjects,1),2),3));
elseif strcmp(analysis,'category_decoding')
    avg_over_conditions_all_subjects = squeeze(nanmean(decoding_accuracies_all_subjects,1));
end

legend_cell(cellfun(@(x) any(isnan(x)),legend_cell)) = [];

%% Plot the average of all subjects
plot(avg_over_conditions_all_subjects-50,'--','Color','k','Linewidth',3);
hold on;
if task==1
    task_title = 'scene categorization';
elseif task==2
    task_title = 'distraction';
end

if strcmp(analysis,'object_decoding')
    analysis_title = 'Object';
elseif strcmp(analysis,'category_decoding')
    analysis_title = 'Category';
end

plot_title = sprintf('%s decoding per timepoint for %d subjects in a %s task',...
    analysis_title,numel(subjects),task_title);
title_boolean = 0;
legend_cell{numel(subjects)+1} = 'Average over all subjects';
legend_cell{numel(subjects)+2} = 'Stimulus onset'; %last legend element
legend_boolean = 0;
legend_location = 'best';
font_size = 8;
set(gca,'FontName','Arial','FontSize',18);
ylim([-10 50]);
y_label = 'Decoding accuracy-50 (%)';
plotting_parameters(plot_title,title_boolean,legend_cell,...
    legend_boolean,font_size,legend_location,y_label); 

%% Save the plot
if strcmp(analysis,'object_decoding')
    filename = 'svm_object_decoding_subjects';
elseif strcmp(analysis,'category_decoding')
    filename = 'svm_artificial_vs_natural_subjects';
end

%Save the plot with or without legend
if legend_boolean == 1
    saveas(gcf,fullfile(results_avg_dir,sprintf('legend_%s_%d_%d_%s',filename,subjects(1),subjects(end),task_name))); 
    saveas(gcf,fullfile(results_avg_dir,sprintf('legend_%s_%d_%d_%s.svg',filename,subjects(1),subjects(end),task_name))); 
else
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s',filename,subjects(1),subjects(end),task_name))); %save as matlab figure
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s.svg',filename,subjects(1),subjects(end),task_name))); %save as svg    
end
close(gcf);

end