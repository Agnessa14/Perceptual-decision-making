function rsa_rnn_eeg(subjects,with_stats) 
%RSA_RNN_EEG Plots the RSA (already computed) between RNN and EEG with stats, if desired.
%
%Input: subjects' ID (e.g., 1:13), add stats (1 with, 0 without)
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
addpath(genpath('/home/agnek95/PyColormap4Matlab'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%% Load and plot the RSA results
layer_idxs = [1,5,7];
num_timepoints_eeg = 200;
num_timepoints_rnn = 8;
rsa_results_dir = fullfile(results_avg_dir,'RSA_matrices_eltanin');
legend_bool = 0;

for c = 1:3 %natural,artificial,all
    if c == 1
        conditions = 'artificial';
        colormap_plot = 'PuRd';
    elseif c == 2
        conditions = 'natural';
        colormap_plot = 'YlGn';
    elseif c == 3
        conditions = 'all';
        colormap_plot = 'OrRd';
    end
    filename = sprintf('Model_RDM_PCA_7_layers_8_timepoints_%s.mat',conditions);
    load(fullfile(rsa_results_dir,filename),'data');
    rsa_results = data;
    cmap = getPyPlot_cMap(colormap_plot,num_timepoints_rnn);
    for l=layer_idxs
        figure;
        legend_plot = cell(num_timepoints_rnn,1);

        for t=1:num_timepoints_rnn
            plot(squeeze(rsa_results(l,t,:)),'LineWidth',2,'Color',cmap(t,:));
            hold on;
            legend_plot{t} = sprintf('Timepoint %s',num2str(t));
        end
        if legend_bool==1
            legend(legend_plot,'Location','best');
        end
        ylim([-0.2 0.6]);
        font_size = 18;
        set(gca,'FontName','Arial','FontSize',font_size);
        xticks(0:20:200);
        xticklabels(-200:100:800);        
        onset_time = 40;
        xline(onset_time,'--');
    end
end
    


%% Load the EEG RDMs
numConditions = 60;
numTimepoints = 200;
rdm_eeg_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints);

% Loop: collect results from all subjects + avg
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'rdm_pearson_categorization.mat'),'rdm_avg');
    rdm_eeg_all_subjects(subject,:,:,:) = rdm_avg; 
%     rdm_1 = squeeze(nanmean(rdm_avg(1,:,:,:),1));
%     rdm_2 = squeeze(nanmean(rdm_avg(2,:,:,:),1));
%     [rdm_rsa(subject,:,:),rdm_1_flattened,rdm_2_flattened] = representational_SA(rdm_1,rdm_2,numTimepoints);
end 

% Remove any NaN (for non-included subjects and conditions & average over subjects
avg_rdm_eeg = squeeze(nanmean(rdm_eeg_all_subjects,1));


%% Correlate each subject's distances with the median RT
t = 1:numTimepoints;
size_corr = [max(subjects),numTimepoints];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);

if modality_distance==1
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
    
elseif modality_distance==2
    
    num_layers=size(distances,1);
    num_timepoints_rnn=size(distances,2);
    size_corr = [num_layers,num_timepoints_rnn];
    correlation_art = NaN(size_corr);
    correlation_nat = NaN(size_corr);
    correlation_both = NaN(size_corr);
    all_conditions = [artificial_conditions natural_conditions];
    
    for l = 1:num_layers
        for tp = 1:num_timepoints_rnn
            correlation_art(l,:) = arrayfun(@(x) ...
                corr(squeeze(distances(l,x,artificial_conditions)),...
                RT_art','type','Spearman'), 1:num_timepoints_rnn);
            correlation_nat(l,:) = arrayfun(@(x) ...
                corr(squeeze(distances(l,x,natural_conditions)),...
                RT_nat','type','Spearman'), 1:num_timepoints_rnn);
            correlation_both(l,:) = arrayfun(@(x) ...
                corr(squeeze(distances(l,x,all_conditions)),...
                RT','type','Spearman'), 1:num_timepoints_rnn);
        end
    end
end


%% Plot
cmap_1 = cool;
cmap_2 = summer;
cmap_3 = autumn;

% color_both = 'k';
if modality_distance==1  
    color_art = cmap_1(200,:); %purple
    color_nat = cmap_2(100,:); %green
    color_both = cmap_3(100,:);
    figure(abs(round(randn*10)));
    plot(avg_corr_art,'LineWidth',2,'Color',color_art);
    hold on;
    plot(avg_corr_nat,'LineWidth',2,'Color',color_nat);
    hold on;
    plot(avg_corr_both,'LineWidth',2,'Color',color_both);
elseif modality_distance==2
    for c = 1:3
        figure;
        legend_plot = cell(num_layers,1);

        if c==1
            corrplot = correlation_art;
            colorplot = cmap_1(1:30:num_layers*30,:);
        elseif c==2
            corrplot = correlation_nat;
            colorplot = cmap_2(1:30:num_layers*30,:);
        elseif c==3
            corrplot = correlation_both;
            colorplot = cmap_3(1:30:num_layers*30,:);
        end
        for l=1:num_layers
            plot(corrplot(l,:),'LineWidth',2,'Color',colorplot(l,:,:));
            hold on;
            legend_plot{l} = sprintf('Layer %s',num2str(l));
        end
        legend(legend_plot,'Location','best');
    end  
 
end

%% Plot stats if needed
if with_stats   
    %define some variables for the stats and the plot
    analysis = 'random_dth';
    permutation_stats.num_perms = 1000;
    permutation_stats.cluster_th = 0.05;
    permutation_stats.significance_th = 0.05;
    permutation_stats.tail = 'left';
    for c = 1:3
        if c == 1
            category = 'artificial';            
            plot_location = -0.15;
            color = color_art;
            for_stats = correlation_art(subjects,:);
        elseif c == 2
            category = 'natural';            
            plot_location = -0.16;
            color = color_nat;
            for_stats = correlation_nat(subjects,:);
        elseif c == 3
            category = 'both'; 
            plot_location = -0.17;
            color = 'k';
            for_stats = correlation_both(subjects,:);
        end

        %Check if stats already exist, otherwise run the stats script
        if modality_distance == 1
            filename_sign = 'rnn_dth_eeg_distances_permutation_stats';
        else
            filename_sign = 'rnn_dth_rnn_distances_permutation_stats_crosstask';
        end     
        filename = fullfile(results_avg_dir,sprintf('%s_%d_%d_%s_task_%s_%s.mat',filename_sign,...
            subjects(1),subjects(end),modality_distance,analysis,category));
        if exist('filename','file')
            load(filename,'permutation_stats');
        else
            [permutation_stats.SignificantMaxClusterWeight,permutation_stats.pValWeight,...
                permutation_stats.SignificantMaxClusterSize,permutation_stats.pValSize] = ...
                permutation_cluster_1sample_weight_alld(for_stats,permutation_stats.num_perms,...
                permutation_stats.cluster_th,permutation_stats.significance_th,permutation_stats.tail); 
        save(filename,'permutation_stats');
        end
        
        %Plot
%         if c < 3
        st = (permutation_stats.SignificantMaxClusterWeight*plot_location); %depending on the stats
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;
%         end
    end
end 

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
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
% legend_plot = {'Artificial scenes','Natural scenes','Both'}; 
% legend(legend_plot);
ylim([-0.2 0.2]);
% legend_bool = 0;
% title_bool = 0;
% plotting_parameters(plot_title,title_bool,legend_plot,legend_bool,font_size,'best','Spearman''s coefficient'); 

%% Save correlations and figures
dth_results.corr_both_categories = avg_corr_both;
dth_results.corr_artificial = avg_corr_art;
dth_results.corr_natural = avg_corr_nat;

save_path = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
if modality_distance==1
    file_name = 'distance_eeg_rt_rnn';
elseif modality_distance==2
    file_name = 'distance_rnn_rt_eeg';
end
model_name = 'model_22.08';

if modality_distance==1
    save(fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_entropy_%s.mat',subjects(1),subjects(end),file_name,model_name,entropy_thresh)),'dth_results');
    saveas(gcf,fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_entropy_%s.svg',subjects(1),subjects(end),file_name,model_name,entropy_thresh))); 
    saveas(gcf,fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_entropy_%s.fig',subjects(1),subjects(end),file_name,model_name,entropy_thresh))); 
else %fix this - should be a loop over layers
    for layer=1:num_layers
        save(fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_layer_%d.mat',subjects(1),subjects(end),file_name,model_name,layer)),'dth_results');
        saveas(gcf,fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_layer_%d.svg',subjects(1),subjects(end),file_name,model_name,layer))); 
        saveas(gcf,fullfile(save_path,sprintf('rnn_dth_subjects_%d_%d_%s_%s_layer_%d.fig',subjects(1),subjects(end),file_name,model_name,layer))); 
    end
close(gcf);
end





