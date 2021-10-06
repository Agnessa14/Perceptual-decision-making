function accuracy_subjects(subjects,task)
%ACCURACY_SUBJECTS Collect, for each participant, their overall task
%accuracy.
%
%Input: subject IDs, task (1=categorization, 2=distraction/fixation) 
%
%Returns 
%- a Sx1 vector of accuracies, where S is the # of participants.


%% Setup the paths and preallocate variables
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
data_dir = '/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT';
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
addpath(genpath(results_dir));
task_name = get_task_name(task);
accuracy = NaN(numel(subjects),1);

%% Get the trial numbers
for s = subjects
    subname = get_subject_name(s);
    load(fullfile(data_dir,subname,sprintf('preprocessed_behavioural_data_%s',task_name)),'behav');
    accuracy(s) = mean(behav.points); 
end

%save
save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
save(fullfile(save_path,sprintf('accuracy_%d_subjects_%s.mat',numel(subjects),task_name)),...
    'accuracy');
end

