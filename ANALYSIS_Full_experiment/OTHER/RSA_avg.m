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
filename_rsa = fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d.mat',subjects(1),subjects(end)));
filename_rdm_cat = fullfile(results_avg_dir,sprintf('average_rdm_categorization_subjects_%d_%d',subjects(1),subjects(end)));
filename_rdm_dis = fullfile(results_avg_dir,sprintf('average_rdm_fixation_subjects_%d_%d',subjects(1),subjects(end)));

if exist(filename_rsa,'file')
    load(filename_rsa,'rdm_rsa');
    load(filename_rdm_cat,'rdm_cat');
    load(filename_rdm_dis,'rdm_dis');
else

    rdm_all_subjects_cat = NaN(max(subjects),numConditions,numConditions,numTimepoints);
    rdm_all_subjects_dis = NaN(max(subjects),numConditions,numConditions,numTimepoints);

    %% Loop: collect results from all subjects 
    for subject = subjects
        subname = get_subject_name(subject);
        subject_results_dir = fullfile(results_dir,subname);
        cat_filename = 'rdm_pearson_categorization.mat'; %svm_decoding_accuracy
        dis_filename = 'rdm_pearson_fixation.mat'; %svm_decoding_accuracy

        load(fullfile(subject_results_dir,cat_filename),'rdm_avg');
        rdm_all_subjects_cat(subject,:,:,:) = rdm_avg;
        load(fullfile(subject_results_dir,dis_filename),'rdm_avg');
        rdm_all_subjects_dis(subject,:,:,:) = rdm_avg;
    end   

    %% Average over subjects and if method 2, perform RSA
    rdm_cat = squeeze(nanmean(rdm_all_subjects_cat,1)); %avg over subjects
    rdm_dis = squeeze(nanmean(rdm_all_subjects_dis,1)); 
    [rdm_rsa,~,~] = representational_SA(rdm_cat,rdm_dis,numTimepoints);
       
    % Save matrices
    save(filename_rdm_cat,'rdm_cat');
    save(filename_rdm_dis,'rdm_dis');
    save(filename_rdm_rsa,'rdm_rsa');
end

%% Plot
h = pcolor(rdm_rsa); 
hold on;
set(h, 'EdgeColor', 'none');
axis square;
cbar = colorbar;
ylabel(cbar,'Spearman''s coefficient');
xlabel('Timepoints: Distraction task');
ylabel('Timepoints: Categorization task');
plot_title = sprintf('Time-generalized RSA of scene processing in categorization and distraction tasks (N=%d)',numel(subjects));
title_bool = 0;
if title_bool==1
    title(plot_title);
end
% caxis([-5 20]);
xticks(0:20:200);
yticks(0:20:200);
xticklabels(-200:100:800);
yticklabels(-200:100:800);
xline(40,'--','Color','w');
yline(40,'--','Color','w');
plot(1:numTimepoints,1:numTimepoints,'LineWidth',2.5,'Color','k');

%% Plot stats if needed
if with_stats

    stats_decoding.num_perms = 1000;
    stats_decoding.qvalue = 0.01;
    stats_decoding.tail = 'right';
    filename = fullfile(results_avg_dir,...
        sprintf('stats_fdr_rsa_cross_task_subjects_%d_%d.mat',subjects(1),subjects(end)));
    if exist(filename,'file')
        load(filename,'stats_decoding');
    else
        [stats_decoding.significant_timepoints,stats_decoding.pvalues]...
            = fdr_rsa_stats(rdm_cat,rdm_dis,rdm_rsa,...
            stats_decoding.num_perms,stats_decoding.tail,stats_decoding.qvalue);
        save(filename,'stats_decoding');
    end

    %plot the stats
    contour(stats_decoding.significant_timepoints,1,'LineColor','w','LineWidth',2);

end

%% Save figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d',subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_cross_task_subjects_%d_%d.png',subjects(1),subjects(end)))); %save as png
close(gcf);

end
    
