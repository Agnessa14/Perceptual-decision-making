function all_subjects_for_stats(subjects,task)
%ALL_SUBJECTS_FOR_STATS Gather the subject x time matrices needed to perform cluster-based permutation tests.
%
%Input: subject IDs, task (1=categorization,2=distraction)
%
%Output: SxP matrix containing the correlations for 

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
