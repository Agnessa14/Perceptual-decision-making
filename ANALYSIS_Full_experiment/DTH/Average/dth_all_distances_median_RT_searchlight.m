function dth_all_distances_median_RT_searchlight(subjects,task_distance,task_RT) %distance art, distance nat, RT
%DTH_ALL_DISTANCES_MEDIAN_RT_SEARCHLIGHT Performs the distance-to-hyperplane searchlight analysis using
%the distances from each subject and the median RT across
%subjects.
%
%Input: subjects' ID (e.g., 1:13), task for the distances (1=categorization, 2=distraction), task for RT,
%add stats (1 with, 0 without)
%
%Correlates the decision values of each subject with reaction times averaged over
%participants of each condition (60 scenes), averaged, at each timepoint, resulting in a plot of Spearman's correlation vs time.
%
%Author: Agnessa Karapetian, 2021 and modified by Muthukumar Pandaram, 2022

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
addpath(genpath(results_dir));
task_name_distance = get_task_name(task_distance);
task_name_RT = get_task_name(task_RT);

%% Get the distances from all subjects
%Define peak times
if task_distance == 1
    times = [110,130,160,165]; %peaks from different analyses
%     all_dimensions = repmat({':'},1,3); %3D: conditions x times x channels
elseif task_distance == 2
    times = [85,105,110,155,165];
end

%Preallocate
numConditions = 60;
numChannels = 63;
chanIdx = 1:numChannels;
distances = NaN(max(subjects),numConditions,numel(times),numChannels);
RTs = NaN(max(subjects),numConditions);
RTs_art = NaN(max(subjects),numConditions/2);
RTs_nat = NaN(size(RTs_art));
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;

%Loop over subjects
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,...
        sprintf('cross_validated_dth_pseudotrials_svm_decisionValues_searchlight_peak_%s.mat',...
        task_name_distance)),'decisionValues_Avg');
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',...
        task_name_RT)),'RT_per_condition');
    distances(subject,:,:,:) = decisionValues_Avg;   
    RTs(subject,:) = RT_per_condition; 
    RTs_art(subject,:) = RT_per_condition(artificial_conditions);
    RTs_nat(subject,:) = RT_per_condition(natural_conditions);
end

%% Get the median RTs of all subjects for each condition
medianRT = nanmedian(RTs,1);
medianRT_art = nanmedian(RTs_art,1);
medianRT_nat = nanmedian(RTs_nat,1);

%% Correlate each subject's distances with the median RT
size_corr = [numel(subjects),numChannels]; 
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);

for subject = subjects
    %Find any missing channels and exclude them
    distances_subject = squeeze(distances(subject,:,:,:));
    missing_channel_ids = find(isnan(distances_subject(1,1,:)));
    for iChan = chanIdx
        if ~ismember(iChan,missing_channel_ids) || isempty(missing_channel_ids)
            
            %Remove any missing conditions
            if any(ismember(artificial_conditions,find(isnan(distances_subject(:,:,:)))))
                artificial_conditions_subject = artificial_conditions(~isnan(distances_subject(artificial_conditions,1,2)));%2nd channel is never missing
                natural_conditions_subject = natural_conditions;
            elseif any(ismember(natural_conditions,find(isnan(distances_subject))))
                artificial_conditions_subject = artificial_conditions;
                natural_conditions_subject = natural_conditions(~isnan(distances_subject(natural_conditions,1,1)));
            else
                artificial_conditions_subject = artificial_conditions;
                natural_conditions_subject = natural_conditions;
            end
            
            %define the peak time
            %within-categorization
            if isequal(task_distance,task_RT) && task_distance==1
                idx_peak_time_artificial = 4; 
                idx_peak_time_natural = 3;
                idx_peak_time_both = 3; 
            %within-distraction
            elseif isequal(task_distance,task_RT) && task_distance==2
                idx_peak_time_artificial = 1; 
                idx_peak_time_natural = 3; 
                idx_peak_time_both = 3; 
            %cross-task, categorization EEG    
            elseif ~isequal(task_distance,task_RT) && task_distance==1
                idx_peak_time_artificial = 1; 
                idx_peak_time_natural = 2;
                idx_peak_time_both = 1;
            %cross-task, distraction EEG
            elseif ~isequal(task_distance,task_RT) && task_distance==2                
                idx_peak_time_artificial = 2; 
                idx_peak_time_natural = 5;
                idx_peak_time_both = 5; 
            end
            
            %Correlate RT and distances
            all_conditions_subject = [artificial_conditions_subject natural_conditions_subject];
            correlation_art(subject,iChan) = corr(squeeze(distances_subject(artificial_conditions_subject,idx_peak_time_artificial,iChan)),...
                medianRT_art(artificial_conditions_subject)','type','Spearman');
            correlation_nat(subject,iChan) = corr(squeeze(distances_subject(natural_conditions_subject,idx_peak_time_natural,iChan)),...
                medianRT_nat(natural_conditions_subject-30)','type','Spearman');
            correlation_both(subject,iChan) = corr(squeeze(distances_subject(all_conditions_subject,idx_peak_time_both,iChan)),...
                medianRT(all_conditions_subject)','type','Spearman');           
        end
    end
end

%% Average over participants
avg_corr_art = squeeze(nanmean(correlation_art,1));
avg_corr_nat = squeeze(nanmean(correlation_nat,1));
avg_corr_both = squeeze(nanmean(correlation_both,1));

%% Save correlations (avg over subjects and all)
dth_results.corr_both_categories = avg_corr_both;
dth_results.corr_artificial = avg_corr_art;
dth_results.corr_natural = avg_corr_nat;
dth_results.for_stats_corr_both_categories = correlation_both;
dth_results.for_stats_corr_artificial= correlation_art;
dth_results.for_stats_corr_natural = correlation_nat;

save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
file_name = 'dth_searchlight_peak';

if ~isequal(task_distance,task_RT)
    file_name = sprintf('%s_cross_task',file_name);
end
save(fullfile(save_path,sprintf('%s_subjects_%d_%d_%s.mat',file_name,subjects(1),subjects(end),task_name_distance)),'dth_results');

end





