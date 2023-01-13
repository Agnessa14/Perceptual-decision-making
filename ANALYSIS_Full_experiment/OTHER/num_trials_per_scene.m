function num_trials_per_scene(subjects,task)
%NUM_TRIALS_PER_SCENE Collect, for each participant, the number of correct
%trials for each scene used in the analysis (based on the minimum number of correct trials for a scene), and save it in mat file. 
%In another file, save the number of trials used per participant.
%
%Input: subject IDs, task (1=categorization, 2=distraction/fixation) 
%
%Returns 
%- a SxN vector of the number of trials for each scence, where S is the
%number of participants and N is the number of scenes.
%- a Sx1 vector of the number of trials used for each participant, S
%being the number of participants.
%

%% Setup the paths and preallocate variables
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
data_dir = '/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT';
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
addpath(genpath(results_dir));
task_name = get_task_name(task);
num_conditions_max = 60;
num_trials_used = NaN(numel(subjects),1);
num_trials_condition = NaN(numel(subjects),num_conditions_max);

%% Get the trial numbers
for s = subjects
    subname = get_subject_name(s);
    load(fullfile(data_dir,subname,sprintf('timelock_%s',task_name)));
    load(fullfile(data_dir,subname,sprintf('preprocessed_behavioural_data_%s',task_name)));
    timelock_triggers = timelock.trialinfo(behav.RT>0 & behav.points==1); 
    [num_trials_used(s),num_trials_condition(s,:)] = min_number_trials(timelock_triggers, num_conditions_max); %minimum number of trials per scene   
end

%save
save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
save(fullfile(save_path,sprintf('num_trials_per_scene_%d_subjects_%s.mat',numel(subjects),task_name)),...
    'num_trials_condition');
save(fullfile(save_path,sprintf('num_trials_used_%d_subjects_%s.mat',numel(subjects),task_name)),...
    'num_trials_used');
end

