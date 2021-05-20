function plot_decoding_full_experiment(subjects,task,analysis)
%PLOT_OBJECT_DECODING_FULL_EXPERIMENT Plot the results from object decoding, averaged over
%all participants and subject-specific.
%
%Input: subject IDs, task (1=categorization,2=distraction), analysis
%('object_decoding' or 'category_decoding')
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

% Set up the figure for plotting
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
sorted_subjects = sort(subjects); %order by ID
if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects = NaN(sorted_subjects(end),numTimepoints);
end
   
legend_cell = cell(1,sorted_subjects(end));
legend_cell(:) = {NaN};

cmap = jet(max(subjects));

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)));
        decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;
    elseif strcmp(analysis,'category_decoding')
        load(fullfile(subject_results_dir,sprintf('svm_artificial_vs_natural_decoding_accuracy_%s.mat',task_name)));
        decoding_accuracies_all_subjects(subject,:) = decodingAccuracy_avg;
    end
    
    %plot the object decoding curve for the participant
    if strcmp(analysis,'object_decoding')
        decodingAccuracy_avg = squeeze(nanmean(nanmean(decodingAccuracy_avg,1),2)); %avg over conditions
    end
    plot(decodingAccuracy_avg, 'Linewidth',2, 'Color', cmap(subject, :));
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
plot(avg_over_conditions_all_subjects,'--','Color','k','Linewidth',3);
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

plot_title = sprintf('%s decoding per timepoint for %d subjects in a %s task',analysis_title,numel(subjects),task_title);
onset_time = 40; 
legend_cell{numel(subjects)+1} = 'Average over all subjects';
legend_cell{numel(subjects)+2} = 'Stimulus onset'; %last legend element
xticks(0:10:200);
plotting_parameters(plot_title,legend_cell,onset_time,8,'best','Decoding accuracy (%)'); %'best' location or specify position [0.75 0.7 0.1 0.1]  

%% Save the plot
if strcmp(analysis,'object_decoding')
    filename = 'svm_object_decoding_subjects';
elseif strcmp(analysis,'category_decoding')
    filename = 'svm_artificial_vs_natural_subjects';
end

%Save the plot with legend
saveas(gcf,fullfile(results_avg_dir,sprintf('legend_%s_%d_%d_%s.svg',filename,subjects(1),subjects(end),task_name))); %save as svg

%Save the plot without legend
leg = legend;
set(leg,'visible','off');
saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s',filename,subjects(1),subjects(end),task_name))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s.svg',filename,subjects(1),subjects(end),task_name))); %save as svg
close(gcf);
end