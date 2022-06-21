function plot_topography_searchlight_dth_ak(subjects,task_distance,task_rt,with_stats)
%PLOT_TOPOGRAPHY_SEARCHLIGHT_DTH_AK Plot the distance-to-hyperplane searchlight results
%at peak. 
%
%Input: subject IDs, task_distance and task_rt (1 for categorization, 2 for distraction),
%with or without stats (1/0)
%
%Output: one topography plot
%
%Author: Agnessa Karapetian, 2022

%% default
obj = 'ctg';

%% Paths etc.
main_path = '/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/';
addpath(genpath(main_path));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%EEGLAB toolbox
eeglab_script = fullfile(main_path,'OTHER/eeglab2022.0','eeglab.m');
run(eeglab_script);

% load channel location
load(fullfile(main_path,'OTHER/adult_63channels.mat'),'channelloc'); 

%% Load Results
% numConditions = 60;
% numChannels = 63; 
task_name_distance = get_task_name(task_distance);
file_name = 'dth_searchlight_peak';
if isequal(task_distance,task_RT)
    file_name = sprintf('%s_cross_task',file_name);
end
file_name_full = sprintf('%s_subjects_%d_%d_%s.mat',file_name,subjects(1),subjects(end),task_name_distance);
load(fullfile(results_avg_dir,file_name_full),'dth_results');

%% Stats
if with_stats
    for conditions = 1:3
        if conditions == 1
            conds_name = 'artificial';    
            for_stats_data = dth_results.for_stats_corr_artificial;
            topo.(obj) = dth_results.corr_artificial;
        elseif conditions == 2
            conds_name = 'natural';    
            for_stats_data = dth_results.for_stats_corr_natural;
            topo.(obj) = dth_results.corr_natural;
        elseif conditions == 3
            conds_name = 'both';    
            for_stats_data = dth_results.for_stats_corr_both_categories;
            topo.(obj) = dth_results.corr_both_categories;
        end

        %Stat parameters
        filename = fullfile(results_avg_dir,...
            sprintf('stats_dth_searchlight_subjects_%d_%d_%s_%s.mat',subjects(1),subjects(end),conds_name,task_distance_name));
        if ~isequal(task_distance,task_name)
            filename = sprintf('%s_cross_task',filename);
        end
        if exist(filename,'file')
            load(filename,'stats_dth_sl');
        else
            stats_dth_sl.num_perms = 10000;
            stats_dth_sl.tail = 'right';
            stats_dth_sl.qvalue = 0.01;
            [stats_dth_sl.significant_channels,stats_dth_sl.pvalues,...
                stats_dth_sl.crit_p, stats_dth_sl.adjusted_pvalues]...
                = fdr_permutation_cluster_1sample_alld(for_stats_data,...
                stats_dth_sl.num_perms,stats_dth_sl.tail,stats_dth_sl.qvalue);
            save(filename,'stats_dth_sl');
        end   
    end
end

for conditions = 1:3
%% Plot searchlight results at peak 
    searchlight_patterns = topo.(obj);
    %Color map
    clim = [0,30];
    warning('off');
    if task == 1
        color_scheme = 'Blues';
    elseif task == 2
        color_scheme = 'PuRd';
    end
    color_upper = cbrewer('seq',color_scheme, 100); 
    colors =  color_upper;
    colors(colors<0) = 0; % added due to negagive values because of change in interpolation method in the newer MATLAB version

    %Set non-significant channels to off
    c=channelloc;
    if with_stats
        topomask = stats_dth_sl.significant_channels;
        for imask = 1:length(topomask)
            if topomask(imask)==0
                c(imask).labels = ' ';
            else
                c(imask).labels = 'O';
            end
        end
    end

    %Plot
    figure;
    topoplot(searchlight_patterns, c, 'colormap', colors, 'style','both','electrodes','labels','whitebk','on');
    caxis(clim);
    if task == 1 && strcmp(analysis,'object_decoding')
        peak_time = 120;
    elseif task == 1 && strcmp(analysis,'category_decoding')
        peak_time = 160;
    elseif task == 2 && strcmp(analysis,'object_decoding')
        peak_time = 110;
    elseif task == 2 && strcmp(analysis,'category_decoding')
        peak_time = 155;
    end
    title([num2str(peak_time),' ms']);


    %Color bar
    CBar_Handle = colorbar('West');
    caxis(clim);
    set(get(CBar_Handle, 'YLabel'), 'String', 'Decoding accuracy (%)-50',...
        'FontSize', 10, 'FontName', 'Arial');
    set(CBar_Handle,'Location','eastoutside');

    %Plotting parameters
    set(gca,'visible', 'off');
    set(gcf, 'color','white');
    fig_width = 600;
    fig_height = 300;
    set(gcf,'position',[100 100 fig_width fig_height]);
    set(gcf,'Renderer','painters');

    %Save
    save(fullfile(results_avg_dir,sprintf('for_stats_cat_svm_%s_subjects_%d_%d_searchlight_peak_%s.mat',analysis,subjects(1),subjects(end),task_name)),'for_stats_cat'); 
    saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_searchlight_peak_%s',analysis,subjects(1),subjects(end),task_name))); %save as matlab figure
    saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_searchlight_peak_%s.svg',analysis,subjects(1),subjects(end),task_name))); %save as svg
    close(gcf);    
end
end 
