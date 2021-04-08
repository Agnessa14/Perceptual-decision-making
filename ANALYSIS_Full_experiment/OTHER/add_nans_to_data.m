function add_nans_to_data(subjects,task)
%ADD_NANS_TO_DATA For the subject that are missing data for certain
%conditions, add NaNs so that their data still has 60 conditions.
%
%Input: subject IDs, task (1=categorization,2=distraction), conditions to
%replace
%
%Returns a data matrix of 60 conditions x dimensions of the data; or 60
%conditions x 1 for the RT vector.

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
task_name = get_task_name(task);

%% Preallocate
numTimepoints = 200;

% 22 
added_matrix_22_row = NaN(1,58,numTimepoints);
added_matrix_22_column = NaN(59,1,numTimepoints);

% 56
added_matrix_56_row = NaN(1,59,numTimepoints);
added_matrix_56_column = NaN(60,1,numTimepoints);

% 47
added_matrix_47_row = NaN(1,59,numTimepoints);
added_matrix_47_column = NaN(60,1,numTimepoints);

%for dth
added_2d_array = NaN(1,numTimepoints);
added_NaN = NaN(1,1);

%% Load, adjust and save the data for each subject
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)));
    load(fullfile(results_dir,subname,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_name)));
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name)));

    if ismember(subject,1:4) && task==2
        %% Object decoding
        if isequal(size(decodingAccuracy_avg),[58,58,200])
            decodingAccuracy_with_22 = cat(1,decodingAccuracy_avg(1:21,:,:),added_matrix_22_row,...
                decodingAccuracy_avg(22:end,:,:));
            decodingAccuracy_with_22 = cat(2,decodingAccuracy_with_22(:,1:21,:),added_matrix_22_column,...
                decodingAccuracy_with_22(:,22:end,:));
            decodingAccuracy_with_56 = cat(1,decodingAccuracy_with_22(1:55,:,:),added_matrix_56_row,...
                decodingAccuracy_with_22(56:end,:,:));
            decodingAccuracy_with_56 = cat(2,decodingAccuracy_with_56(:,1:55,:),added_matrix_56_column,...
                decodingAccuracy_with_56(:,56:end,:));
            %rename
            decodingAccuracy_avg = decodingAccuracy_with_56;
        end
        %% DTH
        if isequal(size(decisionValues_Avg),[58,200])
            decisionValues_Avg_with_22 = cat(1,decisionValues_Avg(1:21,:),added_2d_array,decisionValues_Avg(22:end,:)); 
            decisionValues_Avg_with_56 = cat(1,decisionValues_Avg_with_22(1:55,:),added_2d_array,decisionValues_Avg_with_22(56:end,:));
            RT_per_condition_with_22 = [RT_per_condition(1:21);added_NaN;RT_per_condition(22:end)];
            RT_per_condition_with_56 = [RT_per_condition_with_22(1:55);added_NaN;RT_per_condition_with_22(56:end)];
            decisionValues_Avg = decisionValues_Avg_with_56;
            RT_per_condition = RT_per_condition_with_56;
        end
    elseif subject == 10 && task == 1
        %% OD
        decodingAccuracy_with_47 = cat(1,decodingAccuracy_avg(1:46,:,:),added_matrix_47_row,...
            decodingAccuracy_avg(47:end,:,:));
        decodingAccuracy_with_47 = cat(2,decodingAccuracy_with_47(:,1:46,:),added_matrix_47_column,...
            decodingAccuracy_with_47(:,47:end,:));
        %rename
        decodingAccuracy_avg = decodingAccuracy_with_47;
        
        %% DTH
        decisionValues_Avg = cat(1,decisionValues_Avg(1:46,:),added_2d_array,decisionValues_Avg(47:end,:));
        RT_per_condition = cat(1,RT_per_condition(1:46),added_NaN,RT_per_condition(47:end));
    end

    save(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)),'decodingAccuracy_avg');
    save(fullfile(subject_results_dir,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_name)),'decisionValues_Avg');
    save(fullfile(subject_results_dir,sprintf('RTs_correct_trials_%s.mat',task_name)),'RT_per_condition');
end   

end