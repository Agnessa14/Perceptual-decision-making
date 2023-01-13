function RSA_avg(subjects,task,with_stats) 
%RSA_AVG Perform representational similarity analysis on SVM object decoding. 
%
%Input: subject IDs, with stats (1) or without (0), task
%(1=categorization, 2=distraction, 3=cross task), with or without stats
%(1/0)
%
%Output: 
% 1)NxNxP representational dissimilarity matrices (1-Pearson's
% coefficient), one for each task (N = # conditions, P = # timepoints)
% 2)PxP matrix of 1-Spearman's correlations 
% 3)RSA time-time plot based on 2)
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load existing RDMs or construct them from individual ones
numConditions = 60;
numTimepoints = 200;
if task < 3
    task_name = get_task_name(task);
elseif task == 3
    task_name = 'cross_task';
end

filename_rsa = fullfile(results_avg_dir,sprintf('rsa_%s_subjects_%d_%d.mat',task_name,subjects(1),subjects(end)));
filename_rdms = fullfile(results_avg_dir,sprintf('average_rdms_%s_subjects_%d_%d.mat',task_name,subjects(1),subjects(end)));
    
if exist(filename_rsa,'file') && exist(filename_rdms,'file')
    load(filename_rsa,'rdm_rsa');
    load(filename_rdms,'rdm_avg_subjects');
else
    if task == 3
        rdm_all_subjects = NaN(max(subjects),2,2,numConditions,numConditions,numTimepoints); %2 tasks, 2 halves each
    else
        rdm_all_subjects = NaN(max(subjects),2,numConditions,numConditions,numTimepoints); 
    end
        
    %% Loop: collect results from all subjects 
    for subject = subjects
        subname = get_subject_name(subject);
        subject_results_dir = fullfile(results_dir,subname);
        if task == 3
            cat_filename = 'split_within_task_rdm_pearson_categorization.mat';
            dis_filename = 'split_within_task_rdm_pearson_fixation.mat';
            load(fullfile(subject_results_dir,cat_filename),'rdm_avg');
            rdm_all_subjects(subject,1,:,:,:,:) = rdm_avg;
            load(fullfile(subject_results_dir,dis_filename),'rdm_avg');
            rdm_all_subjects(subject,2,:,:,:,:) = rdm_avg;
         else
            mat_filename = 'split_within_task_rdm_pearson';
            load(fullfile(subject_results_dir,sprintf('%s_%s.mat',mat_filename,task_name)),'rdm_avg');           
            rdm_all_subjects(subject,:,:,:,:) = rdm_avg;
        end      
    end   

    %% Average over subjects and perform RSA
    rdm_avg_subjects = squeeze(nanmean(rdm_all_subjects,1)); %avg over subjects
    if task < 3
        rdm_1 = squeeze(rdm_avg_subjects(1,:,:,:));
        rdm_2 = squeeze(rdm_avg_subjects(2,:,:,:));
        rdm_rsa = representational_SA(rdm_1,rdm_2,numTimepoints);
    else 
        rdm_rsa = NaN(2,2,numTimepoints,numTimepoints);
        for h = 1:2
            for j = 1:2
                rdm_1 = squeeze(rdm_avg_subjects(1,h,:,:,:));
                rdm_2 = squeeze(rdm_avg_subjects(2,j,:,:,:));
                rdm_rsa(h,j,:,:) = representational_SA(rdm_1,rdm_2,numTimepoints);
            end
        end
        rdm_rsa = squeeze(mean(rdm_rsa,1:2));
    end
       
    % Save matrices
    save(filename_rdms,'rdm_avg_subjects');
    save(filename_rsa,'rdm_rsa');
end

%% Plot
h = pcolor(rdm_rsa); 
hold on;
set(h, 'EdgeColor', 'none');
axis square;
cbar = colorbar;
ylabel(cbar,'Spearman''s coefficient');
if task == 3
    xlabel('Timepoints: Distraction task');
    ylabel('Timepoints: Categorization task');
else
    xlabel('Timepoints tested on');
    ylabel('Timepoints trained on');
end
plot_title = sprintf('Time-generalized RSA of scene processing, %s (N=%d)',task_name,numel(subjects));
title_bool = 0;
if title_bool==1
    title(plot_title);
end
caxis([-0.2 1]);
xticks(0:20:200);
yticks(0:20:200);
xticklabels(-200:100:800);
yticklabels(-200:100:800);
xline(40,'--','Color','w');
yline(40,'--','Color','w');
set(gca,'FontName','Arial','FontSize',18);

if task == 3
    plot(1:numTimepoints,1:numTimepoints,'--','LineWidth',2.5,'Color','k');
end

%% Plot stats if needed
if with_stats

    if task < 3
        rdm_1 = squeeze(rdm_avg_subjects(1,:,:,:));
        rdm_2 = squeeze(rdm_avg_subjects(2,:,:,:));
    else
        rdm_1 = squeeze(mean(rdm_avg_subjects(1,:,:,:,:),2)); %avg over halves
        rdm_2 = squeeze(mean(rdm_avg_subjects(2,:,:,:,:),2));
    end


 

    stats_decoding.num_perms = 1000;
    stats_decoding.qvalue = 0.01;
    stats_decoding.tail = 'right';
    filename = fullfile(results_avg_dir,...
        sprintf('stats_fdr_rsa_%s_subjects_%d_%d.mat',task_name,subjects(1),subjects(end)));
    if exist(filename,'file')
        load(filename,'stats_decoding');
    else
        [stats_decoding.significant_timepoints,stats_decoding.pvalues]...
            = fdr_rsa_stats(rdm_1,rdm_2,rdm_rsa,...
            stats_decoding.num_perms,stats_decoding.tail,stats_decoding.qvalue);
        save(filename,'stats_decoding');
    end

    %plot the stats
    contour(stats_decoding.significant_timepoints,1,'LineColor','w','LineWidth',2);
end

%% Save figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_%s_subjects_%d_%d',task_name,subjects(1),subjects(end)))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('rsa_%s_subjects_%d_%d.png',task_name,subjects(1),subjects(end)))); %save as png
close(gcf);

end
    
