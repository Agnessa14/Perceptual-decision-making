function RSA_cross_task(subjects) 
%RSA_CROSS_TASK Perform representational similarity analysis on the data from both tasks. 
%
%Input: subject ID
%
%Output: 1xP vector of correlations in, where P is the number of timepoints. 
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
decoding_accuracies_all_subjects_fix = NaN(sorted_subjects(end),numConditions,numConditions,numTimepoints);


%% Loop: collect results from all subjects 
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    cat_filename = 'svm_decoding_accuracy_categorization.mat';
    fix_filename = 'svm_decoding_accuracy_fixation.mat';

    load(fullfile(subject_results_dir,cat_filename));
    decoding_accuracies_all_subjects_cat(subject,:,:,:) = decodingAccuracy_avg;
    load(fullfile(subject_results_dir,fix_filename));
    decoding_accuracies_all_subjects_fix(subject,:,:,:) = decodingAccuracy_avg;
    
end   

%% RSA: correlate the RDMs from both tasks at every timepoint 
average_DA_cat = squeeze(nanmean(decoding_accuracies_all_subjects_cat,1)); %avg over subjects
average_DA_dis = squeeze(nanmean(decoding_accuracies_all_subjects_cat,1)); %avg over subjects

%fill up the RDM so it's symmetrical
t = 1:numTimepoints;
% rdm_cat = NaN(size(average_DA_cat));
% rdm_dis = NaN(size(average_DA_dis));
rdm_cat = average_DA_cat;
rdm_dis = average_DA_dis;

rdm_cat = reshape(rdm_cat(isnan(rdm_cat))==0,size(rdm_cat,1),size(rdm_cat,2),size(rdm_cat,3));
rdm_dis = reshape(isnan(rdm_dis)==0,size(rdm_dis));

rdm_cat = arrayfun(@(x) squareform(rdm_cat(:,:,x)+(rdm_cat(:,:,x))'), t);
rdm_dis = arrayfun(@(x) squareform(rdm_dis(:,:,x)+(rdm_dis(:,:,x))'), t);

% rdm_cat = arrayfun(@(x) squareform(average_DA_cat(:,:,x)),t);
% rdm_dis = arrayfun(@(x) squareform(average_DA_dis(:,:,x)),t);
rsa = arrayfun(@(x) corr(rdm_cat(t),rdm_dis(t),'type','Spearman'),t);

%% Plot
figure(abs(round(randn*10))); %Random figure number
set(gcf, 'Position', get(0, 'Screensize'));
plot(rsa,'LineWidth',3);
xlabel('Timepoint');
ylabel('Spearman''s coefficient');
title(sprintf('Representational similarity analysis of scene processing in categorization and distraction tasks (N=%d)',numel(subjects)));

%% Save
save(fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end))),'rsa');
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d.svg',subjects(1),subjects(end)))); %save as svg
close(gcf);
end
    
