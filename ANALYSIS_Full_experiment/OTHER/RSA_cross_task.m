function RSA_cross_task(subjects,with_stats) 
%RSA_CROSS_TASK Perform representational similarity analysis on SVM object decoding from both tasks. 
%
%Input: subject IDs, with stats (1) or without (0)
%
%Output: 
% 1)NxNxP representational dissimilarity matrices (1-Pearson's
% coefficient), one for each task
% 2)PxP matrix of 1-Spearman's correlations, where P is the number of timepoints. 
% 3)RSA time-time plot based on 2)
% 4)If with stats, NxN plot of statistical significance
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Preallocate
numConditions = 60;
numTimepoints = 200;
sorted_subjects = sort(subjects);
decoding_accuracies_all_subjects_cat = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);
decoding_accuracies_all_subjects_dis = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);

%% Loop: collect results from all subjects 
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    cat_filename = 'rdm_pearson_categorization.mat'; %svm_decoding_accuracy
    dis_filename = 'rdm_pearson_fixation.mat'; %svm_decoding_accuracy

    load(fullfile(subject_results_dir,cat_filename));
    decoding_accuracies_all_subjects_cat(subject,:,:,:) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,dis_filename));
    decoding_accuracies_all_subjects_dis(subject,:,:,:) = decodingAccuracy_avg;
end   

%% RSA: correlate the RDMs from both tasks at every timepoint 
rdm_cat = squeeze(nanmean(decoding_accuracies_all_subjects_cat,1)); %avg over subjects
rdm_dis = squeeze(nanmean(decoding_accuracies_all_subjects_dis,1)); 

% Fill up the RDM so it's symmetrical: replace the NaNs by 0s & add the transpose of the upper diagonal
rdm_cat(isnan(rdm_cat)) = 0;
rdm_dis(isnan(rdm_dis)) = 0;

rdm_cat_flattened_cell = arrayfun(@(x) squareform(rdm_cat(:,:,x)+(rdm_cat(:,:,x))'),...
    1:numTimepoints,'UniformOutput',false);
rdm_dis_flattened_cell = arrayfun(@(x) squareform(rdm_dis(:,:,x)+(rdm_dis(:,:,x))'),...
    1:numTimepoints,'UniformOutput',false);

rdm_cat_flattened = reshape(cell2mat(rdm_cat_flattened_cell),[],numTimepoints);
rdm_dis_flattened = reshape(cell2mat(rdm_dis_flattened_cell),[],numTimepoints);

% Perfom RSA at each combination of timepoints
rsa = NaN(numTimepoints,numTimepoints);
for tp1 = 1:numTimepoints
    for tp2 = 1:numTimepoints
        rsa(tp1,tp2) = 1-corr(rdm_cat_flattened(:,tp1),rdm_dis_flattened(:,tp2),'type','Spearman');
    end
end

%% Plot
h = pcolor(rsa); 
set(gcf, 'Position', get(0, 'Screensize'));
set(h, 'EdgeColor', 'none');
axis square;
cbar = colorbar;
ylabel(cbar,'1-Spearman''s coefficient');
xlabel('Timepoints: Distraction task');
ylabel('Timepoints: Categorization task');
title(sprintf('Time-generalized RSA of scene processing in categorization and distraction tasks (N=%d)',numel(subjects)));

%% Save RDMs, RSA results and figures
%Matrices
save(fullfile(results_avg_dir,sprintf('average_rdm_categorization_subjects_%d_%d',subjects(1),subjects(end))),'rdm_cat');
save(fullfile(results_avg_dir,sprintf('average_rdm_fixation_subjects_%d_%d',subjects(1),subjects(end))),'rdm_dis');
save(fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end))),'rsa');
%Figures
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d.png',subjects(1),subjects(end)))); %save as png
close(gcf);

%% Plot stats if needed
if with_stats
    %load the subjects x timepoints matrix
    task = 3; %for cross task
    filename_forstats = fullfile(results_avg_dir,sprintf('for_stats_subjects_%d_%d_cross_task_time_object_decoding.mat',...
    subjects(1),subjects(end)));
    if exist(filename_forstats,'file')
        load(filename_forstats);
    else
        for_stats = all_subjects_for_stats(subjects,task,'rsa_time_object');
    end
    
    %load the stats if they have been run already
    filename_sign = fullfile(results_avg_dir,sprintf('significant_timepoints_subjects_%d_%d_%s_cross_task_time_object_decoding.mat',...
    subjects(1),subjects(end)));
    if exist(filename_sign,'file')
        load(filename_sign);
    else
        significant_timepoints = run_permutation_stats(subjects,task,'rsa_time_object',for_stats);
    end
    
    %plot the stats - different plot from the RSA
    significant_timepoints(isnan(significant_timepoints))=0;
    h = pcolor(significant_timepoints); 
    set(gcf, 'Position', get(0, 'Screensize'));
    set(h, 'EdgeColor', 'none');
    axis square;
    cbar = colorbar;
    ylabel(cbar,'Significance');
    xlabel('Timepoints: Distraction task');
    ylabel('Timepoints: Categorization task');
    title(sprintf('Statistical analysis of time-generalized RSA of scene processing in categorization and distraction tasks (N=%d)',numel(subjects)));  
    
    %save the plot
    saveas(gcf,fullfile(results_avg_dir,sprintf('stats_rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end)))); 
    saveas(gcf,fullfile(results_avg_dir,sprintf('stats_rsa_cross_task_subjects_%d_%d.png',subjects(1),subjects(end)))); 
    close(gcf);
end 

end
    
