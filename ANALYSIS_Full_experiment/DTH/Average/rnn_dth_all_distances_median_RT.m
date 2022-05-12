function rnn_dth_all_distances_median_RT(subjects,with_stats,with_error_bars,varargin) %distance art, distance nat, RT
%RNN_DTH_ALL_DISTANCES_MEDIAN_RT Performs the cross-modal distance-to-hyperplane analysis using
%the distances-to-hyperplane from each subject and the RTs from the RNN.
%
%Input: subjects' ID (e.g., 1:13), add stats (1 with, 0 without),
%with/without error bars (1/0), varargin: stats_type ('cluster' or 'fdr')
%
%Correlates the decision values of each subject with reaction times of each condition (60 scenes), at each timepoint,
%resulting in a plot of Spearman's correlation vs time. 
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%% Get the distances 
numTimepoints = 200;
numConditions = 60;
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;  
entropy_thresh = '0.02';
model_name = 'model_02.11_2';

%load distances
distances = NaN(max(subjects),numConditions,numTimepoints);   
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,...
        'cross_validated_dth_pseudotrials_svm_decisionValues_categorization.mat'),...
        'decisionValues_Avg');
    distances(subject,:,:) = decisionValues_Avg;   
end

%load RTs
load(sprintf('/scratch/agnek95/PDM/DATA/RNN_RTs/RNN_RTs_entropy_threshold_%s.mat',...
    entropy_thresh),'data');
RT = data';
RT_art = RT(artificial_conditions);
RT_nat = RT(natural_conditions);

%% Correlate each subject's distances with the median RT
t = 1:numTimepoints;
size_corr = [max(subjects),numTimepoints];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);

for subject = subjects
    %Find any missing scenes and exclude them
    distances_subject = squeeze(distances(subject,:,:));
    if any(ismember(artificial_conditions,find(isnan(distances_subject))))
        artificial_conditions_subject = artificial_conditions(~isnan(distances_subject(artificial_conditions,1)));
        natural_conditions_subject = natural_conditions;
    elseif any(ismember(natural_conditions,find(isnan(distances_subject))))
        artificial_conditions_subject = artificial_conditions;
        natural_conditions_subject = natural_conditions(~isnan(distances_subject(natural_conditions,1)));
    else
        artificial_conditions_subject = artificial_conditions;
        natural_conditions_subject = natural_conditions;
    end

    %Correlate RT and distances    
    all_conditions_subject = [artificial_conditions_subject natural_conditions_subject];
    correlation_art(subject,:) = arrayfun(@(x) ...
        corr(squeeze(distances_subject(artificial_conditions_subject,x)),...
        RT_art(artificial_conditions_subject),'type','Spearman'), t);
    correlation_nat(subject,:) = arrayfun(@(x) ...
        corr(squeeze(distances_subject(natural_conditions_subject,x)),...
        RT_nat(natural_conditions_subject-30),'type','Spearman'), t);
    correlation_both(subject,:) = arrayfun(@(x) ...
        corr(squeeze(distances_subject(all_conditions_subject,x)),...
        RT(all_conditions_subject),'type','Spearman'), t);
end

%% Average over participants
avg_corr_art = squeeze(nanmean(correlation_art,1));
avg_corr_nat = squeeze(nanmean(correlation_nat,1));
avg_corr_both = squeeze(nanmean(correlation_both,1));

%% Plot colours
cmap_1 = cool;
cmap_2 = summer;

color_art = cmap_1(200,:); %purple
color_nat = cmap_2(100,:); %green
color_both = 'k';

figure(abs(round(randn*10)));

%% Plot stats if needed
if with_stats   
    analysis = 'random_dth';
    %define some variables for the stats and the plot
    if isempty(varargin)
        error('Specify the type of stats');
    else
        stats_type = varargin{1};
        if strcmp(stats_type,'cluster')
            %define some variables for the stats and the plot
            permutation_stats.num_perms = 10000;
            permutation_stats.cluster_th = 0.05;
            permutation_stats.significance_th = 0.05;
            permutation_stats.tail = 'left';
        elseif strcmp(stats_type,'fdr')
            fdr_stats.num_perms = 10000;
            fdr_stats.tail = 'left';
            fdr_stats.qvalue = 0.05;
        end
    end     

    for c = 1:3
        if c == 1
            category = 'artificial';            
            plot_location = -0.22;
            color = color_art;
            if with_error_bars
                for_stats = correlation_art(subjects,:);
                data = avg_corr_art;
            end
        elseif c == 2
            category = 'natural';            
            plot_location = -0.24;
            color = color_nat;
            if with_error_bars
                for_stats = correlation_nat(subjects,:);
                data = avg_corr_nat;
            end
        elseif c == 3
            category = 'both'; 
            plot_location = -0.26;
            color = 'k';
            if with_error_bars
                for_stats = correlation_both(subjects,:);
                data = avg_corr_both;
            end
        end

        %Check if stats already exist, otherwise run the stats script
        distances_str = 'eeg';
        filename_sign = sprintf('separate_fitting_cv_%s_rnn_dth_eeg_distances_permutation_stats',model_name);       
        filename = fullfile(results_avg_dir,sprintf('%s_%d_%d_distances_%s_%s_%s.mat',filename_sign,...
            subjects(1),subjects(end),distances_str,analysis,category));
        if exist(filename,'file')
            if strcmp(stats_type,'cluster')
                load(filename,'permutation_stats');
            elseif strcmp(stats_type,'fdr')
                load(filename,'fdr_stats');
            end
        else
            if strcmp(stats_type,'cluster')
                [permutation_stats.SignificantMaxClusterWeight,permutation_stats.pValWeight,...
                    permutation_stats.SignificantMaxClusterSize,permutation_stats.pValSize] = ...
                    permutation_cluster_1sample_weight_alld(for_stats,permutation_stats.num_perms,...
                    permutation_stats.cluster_th,permutation_stats.significance_th,permutation_stats.tail); 
                save(filename,'permutation_stats');
            elseif strcmp(stats_type,'fdr')
                [fdr_stats.significant_timepoints,fdr_stats.pvalues,...
                    fdr_stats.crit_p, fdr_stats.adjusted_pvalues]...
                    = fdr_permutation_cluster_1sample_alld(for_stats,...
                    fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);
                save(filename,'fdr_stats');
            end
        end
       
        %Plot
        %1) significance
        if strcmp(stats_type,'cluster')
            st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
        elseif strcmp(stats_type,'fdr')
            st = (fdr_stats.significant_timepoints*plot_location); %depending on the stats
        end
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;
        
        %if needed:
        if with_error_bars
            %2) error bars
            stdDM = std(for_stats); 
            err = stdDM/sqrt(size(for_stats,1)); %standard deviation/sqrt of num subjects  

            %plot as a shaded area
            top_curve = data + err;
            bottom_curve = data - err;
            x2 = [1:numTimepoints, fliplr(1:numTimepoints)];
            shaded_area = [top_curve, fliplr(bottom_curve)];
            fill(x2, shaded_area, color,'FaceAlpha',0.5);
            hold on;
        end
    end
end 


%% Plot the data
plot(avg_corr_art,'LineWidth',2,'Color',color_art);
hold on;
plot(avg_corr_nat,'LineWidth',2,'Color',color_nat);
hold on;
plot(avg_corr_both,'LineWidth',2,'Color',color_both);

%% Plotting parameters
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
% legend_plot = {'Artificial scenes','Natural scenes','Both'}; 
% legend(legend_plot);
ylim([-0.3 0.3]);
yticks(-0.3:0.1:0.3);
xticks(0:20:200);
xticklabels(-200:100:800);  
xline(40,'--');
% legend_bool = 0;
% title_bool = 0;
% plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Spearman''s coefficient'); 

%% Save correlations and figures
dth_results.corr_both_categories = avg_corr_both;
dth_results.corr_artificial = avg_corr_art;
dth_results.corr_natural = avg_corr_nat;

save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
file_name = 'distance_eeg_rt_rnn';
if with_error_bars
    file_name = sprintf('error_bars_%s',file_name);
end
if ~isempty(varargin)
    if strcmp(stats_type,'cluster')
        file_name = sprintf('%s_cluster',file_name);
    elseif strcmp(stats_type,'fdr')
        file_name = sprintf('%s_fdr',file_name);
    end
end
save(fullfile(save_path,sprintf('separate_fitting_cv_rnn_dth_subjects_%d_%d_%s_%s_entropy_%s.mat',subjects(1),subjects(end),file_name,model_name,entropy_thresh)),'dth_results');
saveas(gcf,fullfile(save_path,sprintf('sf_cv_rnn_dth_subjects_%d_%d_%s_%s_entropy_%s.svg',subjects(1),subjects(end),file_name,model_name,entropy_thresh))); 
saveas(gcf,fullfile(save_path,sprintf('sf_cv_rnn_dth_subjects_%d_%d_%s_%s_entropy_%s.fig',subjects(1),subjects(end),file_name,model_name,entropy_thresh))); 

close(gcf);

end



