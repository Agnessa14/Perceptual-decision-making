function encoding_all_subjects(subjects,task_eeg,task_RT,with_stats,with_error_bars,varargin) %distance art, distance nat, RT
%ENCODING_ALL_SUBJECTS Average encoding accuracy across subjects.
%
%Input: subjects' ID (e.g., 1:13), task for the distances (1=categorization, 2=distraction), task for RT,
%add stats (1 with, 0 without), plot with/without error bars (1/0),
%varargin: stats_type ('cluster' or 'fdr')
%
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));
task_name_eeg = get_task_name(task_eeg);
% task_name_RT = get_task_name(task_RT);

%% Get the distances from all subjects
numTimepoints = 200;
numConditions = 60;
encoding_accuracies = NaN(max(subjects),numTimepoints);
% artificial_conditions = 1:numConditions/2;
% natural_conditions = (numConditions/2)+1:numConditions;
filename = sprintf('cross_validated_regression_RTs_encodingAccuracy_%s',task_name_eeg);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,filename),'encodingAccuracy_avg');
    encoding_accuracies(subject,:) = encodingAccuracy_avg;   
end

%% Average over participants
avg_encoding_acc = squeeze(nanmean(encoding_accuracies,1));

%% Plot
plot(avg_encoding_acc);


% %% Plot
% figure(abs(round(randn*10)));
% cmap_1 = cool;
% cmap_2 = summer;
% color_art = cmap_1(200,:); %purple
% color_nat = cmap_2(100,:); %green
% 
% %% Plot stats if needed
% for c = 1:3
%     if c == 1
%         data = avg_corr_art;
%         color = color_art;
%     elseif c == 2
%         data = avg_corr_nat;
%         color = color_nat;
%     elseif c == 3
%         data = avg_corr_both;
%         color = 'k';
%     end
%     if ~with_stats
%         %Plot average curve
%         plot(data,'LineWidth',2,'Color',color);
%         hold on;
%     else
%         analysis = 'random_dth';
%         if isempty(varargin)
%             error('Specify the type of stats')
%         else
%             stats_type = varargin{1}; 
% 
%             if strcmp(stats_type,'cluster')
%                 %define some variables for the stats and the plot
%                 permutation_stats.num_perms = 10000;
%                 permutation_stats.cluster_th = 0.05;
%                 permutation_stats.significance_th = 0.05;
%                 permutation_stats.tail = 'left';
%             elseif strcmp(stats_type,'fdr')
%                 fdr_stats.num_perms = 10000;
%                 fdr_stats.tail = 'left';
%                 fdr_stats.qvalue = 0.05;
%             end
%         end
% 
%         if c == 1
%             category = 'artificial';            
%             plot_location = 0.15;
%             for_stats = correlation_art(subjects,:);
%             if with_error_bars
%                 for_stats = correlation_art(subjects,:);
%             end
%         elseif c == 2
%             category = 'natural';            
%             plot_location = 0.17;
%             for_stats = correlation_nat(subjects,:);
%             if with_error_bars
%                 for_stats = correlation_nat(subjects,:);
%             end
%         elseif c == 3
%             plot_location = 0.19;
%             category = 'both';            
%             for_stats = correlation_both(subjects,:);       
%         end
%         
%         %Check if stats already exist, otherwise run the stats script
%         if task_distance==task_RT
%             filename_sign = 'dth_permutation_stats';
%         else
%             filename_sign = 'dth_permutation_stats_crosstask';
%         end     
% 
%         filename = fullfile(results_avg_dir,sprintf('%s_%s_%d_%d_%s_task_%s_%s.mat',stats_type,filename_sign,...
%             subjects(1),subjects(end),task_name_distance,analysis,category));
%         if exist(filename,'file')
%             if strcmp(stats_type,'cluster')
%                 load(filename,'permutation_stats');
%             elseif strcmp(stats_type,'fdr')
%                 load(filename,'fdr_stats');
%             end
%         else
%             if strcmp(stats_type,'cluster')
%                 [permutation_stats.SignificantMaxClusterWeight,permutation_stats.pValWeight,...
%                     permutation_stats.SignificantMaxClusterSize,permutation_stats.pValSize] = ...
%                     permutation_cluster_1sample_weight_alld(for_stats,permutation_stats.num_perms,...
%                     permutation_stats.cluster_th,permutation_stats.significance_th,permutation_stats.tail); 
%                 save(filename,'permutation_stats');
%             elseif strcmp(stats_type,'fdr')
%                 [fdr_stats.significant_timepoints,fdr_stats.pvalues,...
%                     fdr_stats.crit_p, fdr_stats.adjusted_pvalues]...
%                     = fdr_permutation_cluster_1sample_alld(for_stats,...
%                     fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);
%                 save(filename,'fdr_stats');
%             end
%         end
%         %Plot statistics
%         if c < 3
%             if strcmp(stats_type,'cluster')
%                 st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
%             elseif strcmp(stats_type,'fdr')
%                 st = (fdr_stats.significant_timepoints*plot_location); %depending on the stats
%             end
%             st(st==0) = NaN;
%             plot(st,'*','Color',color); 
%             hold on;
%         end
%         
%         %if needed: plot error bars
%         if c < 3
%             if with_error_bars
%                 %calculate error bars
%                 stdDM = std(for_stats); 
%                 err = stdDM/sqrt(size(for_stats,1)); %standard deviation/sqrt of num subjects  
% 
%                 %plot as a shaded area
%                 top_curve = data + err;
%                 bottom_curve = data - err;
%                 x2 = [1:numTimepoints, fliplr(1:numTimepoints)];
%                 shaded_area = [top_curve, fliplr(bottom_curve)];
%                 fill(x2, shaded_area, color,'FaceAlpha',0.5);
%                 hold on;
%             end
%             %Plot average curve
%             plot(data,'LineWidth',2,'Color',color);
%             hold on;
%         end
%     end
% end 
%     
%% Plotting parameters
% if task_distance==task_RT
%     if task_distance==1
%         task_title = 'scene categorization';
%     elseif task_distance==2
%         task_title = 'distraction';
%     end
%     plot_title =  sprintf('Correlation between the distance to hyperplane and reaction time in a %s task (N=%d)',task_title,numel(subjects));
% elseif task_distance==1 && task_RT==2
%     plot_title =  sprintf('Correlation between the distance to hyperplane from a categorization task and reaction time from a distraction task (N=%d)',numel(subjects));
% elseif task_distance==2 && task_RT==1
%     plot_title =  sprintf('Correlation between the distance to hyperplane from a distraction task and reaction time from a categorization task (N=%d)',numel(subjects));
% end
plot_title = 'Prediction of RTs based on EEG patterns';
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
legend_plot = {'Artificial scenes','Natural scenes'}; 
% ylim([-0.3 0.3]);
% yticks(-0.3:0.1:0.3);
legend_bool = 0;
title_bool = 1;
plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Encoding accuracy'); 


%% Save
filename_save = sprintf('average_%s',filename);
save(fullfile(filename_save,'avg_encoding_acc'));
saveas(gcf,fullfile(results_avg_dir,filename_save));
saveas(gcf,fullfile(results_avg_dir,sprintf('%s.svg',filename_save)));

close(gcf);

% %% Save correlations and figures
% dth_results.corr_both_categories = avg_corr_both;
% dth_results.corr_artificial = avg_corr_art;
% dth_results.corr_natural = avg_corr_nat;
% 
% save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
% file_name = 'cv_all_dist_med_rt_dth';
% if with_error_bars
%     file_name = sprintf('error_bars_%s',file_name);
% end
% if ~isempty(varargin)
%     if strcmp(stats_type,'cluster')
%         file_name = sprintf('%s_cluster',file_name);
%     elseif strcmp(stats_type,'fdr')
%         file_name = sprintf('%s_fdr',file_name);
%     end
% end
% if isequal(task_distance,task_RT)
%     save(fullfile(save_path,sprintf('%s_subjects_%d_%d_%s.mat',file_name,subjects(1),subjects(end),task_name_distance)),'dth_results');
%     saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_%s',file_name,subjects(1),subjects(end),task_name_distance))); 
%     saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_%s.svg',file_name,subjects(1),subjects(end),task_name_distance))); 
% else
%     save(fullfile(save_path,sprintf('%s_subjects_%d_%d_cross_task_%s_distances.mat',file_name,subjects(1),subjects(end),task_name_distance)),'dth_results');
%     saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_cross_task_%s_distances',file_name,subjects(1),subjects(end),task_name_distance)));
%     saveas(gcf,fullfile(save_path,sprintf('%s_subjects_%d_%d_cross_task_%s_distances.svg',file_name,subjects(1),subjects(end),task_name_distance)));
% end
% 
% close(gcf);
% end





