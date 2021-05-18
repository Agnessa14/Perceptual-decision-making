function RSA_cross_task(subjects) 
%RSA_CROSS_TASK Perform representational similarity analysis on SVM object decoding from both tasks. 
%
%Input: subject IDs
%
%Output: 
% -NxNxP representational dissimilarity matrices (1-Pearson's
% coefficient), one for each task
% -PxP matrix of 1-Spearman's correlations, where P is the number of timepoints. 
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
save(fullfile(results_avg_dir,sprintf('average_rdm_categorization_subjects_%d_%d',subjects(1),subjects(end))),'rdm_cat');
save(fullfile(results_avg_dir,sprintf('average_rdm_fixation_subjects_%d_%d',subjects(1),subjects(end))),'rdm_dis');
save(fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end))),'rsa');
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d.png',subjects(1),subjects(end)))); %save as png
close(gcf);
end
    
