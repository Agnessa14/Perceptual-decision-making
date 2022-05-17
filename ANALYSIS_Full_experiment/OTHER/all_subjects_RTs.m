function all_subjects_RTs(subjects,task)
%ALL_SUBJECTS_RTs Collects the RTs of each subject and saves in one matrix.
%
%Input: subjects' ID (e.g., 1:13), task (1=categorization, 2=distraction)
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
numConditions = 60;
RTs = NaN(max(subjects),numConditions);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name)),'RT_per_condition');
    RTs(subject,:) = RT_per_condition; 
end

%Remove the NaN subjects
RTs = RTs(subjects,:);

%% Save
save(fullfile(results_avg_dir,sprintf('RT_all_subjects_%d_%d_%s.mat',subjects(1),subjects(end),task_name)),'RTs');

end





