function plot_topography_searchlight_dth_ak(subjects,task_distance,task_RT,with_stats)
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
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%EEGLAB toolbox
eeglab_script = fullfile(main_path,'OTHER/eeglab2022.0','eeglab.m');
run(eeglab_script);

% load channel location
load(fullfile(main_path,'OTHER/adult_63channels.mat'),'channelloc'); 

%% Load Results
% numConditions = 60;
% numChannels = 63; 
task_distance_name = get_task_name(task_distance);
file_name = 'dth_searchlight_peak';
if isequal(task_distance,task_RT)
    file_name = sprintf('%s_cross_task',file_name);
end
file_name_full = sprintf('%s_subjects_%d_%d_%s.mat',file_name,subjects(1),subjects(end),task_distance_name);
load(fullfile(results_avg_dir,file_name_full),'dth_results');

%% Stats
for conditions = 1:3
    if conditions == 1
        conds_name = 'artificial';    
        for_stats_data = dth_results.for_stats_corr_artificial(subjects,:);
        topo.(obj) = dth_results.corr_artificial;
        color_scheme = 'PuOr';
        func_color = 'flipud';
        if task_distance==1 && task_RT==1
            peak_time = 165;
        elseif task_distance==1 && task_RT==2
            peak_time = 110;
        elseif task_distance==2 && task_RT==2
            peak_time = 85;
        elseif task_distance==2 && task_RT==1
            peak_time = 105;
        end
    elseif conditions == 2
        conds_name = 'natural';    
        for_stats_data = dth_results.for_stats_corr_natural(subjects,:);
        topo.(obj) = dth_results.corr_natural;
        color_scheme = 'PRGn'; %green
        func_color = 'flipud';
        if task_distance==1 && task_RT==1
            peak_time = 160;
        elseif task_distance==1 && task_RT==2
            peak_time = 130;
        elseif task_distance==2 && task_RT==2
            peak_time = 110;
        elseif task_distance==2 && task_RT==1
            peak_time = 165;
        end
    elseif conditions == 3
        conds_name = 'both';    
        for_stats_data = dth_results.for_stats_corr_both_categories(subjects,:);
        topo.(obj) = dth_results.corr_both_categories;
        if task_distance==1 && task_RT==1
            color_scheme = 'PuOr'; %orange
            func_color = '';
            peak_time = 160;
        elseif task_distance==1 && task_RT==2
            color_scheme = 'YlGnBu'; %turquoise
            func_color = '';
            peak_time = 110;
        elseif task_distance==2 && task_RT==2
            color_scheme = 'RdBu';%blue
            func_color = 'flipud';
            peak_time = 110;
        elseif task_distance==2 && task_RT==1
            color_scheme = 'YlOrRd'; %yellow        
            func_color = '';
            peak_time = 165;
        end
    end
    
    if with_stats
        %Stat parameters
        filename = fullfile(results_avg_dir,...
            sprintf('stats_dth_searchlight_subjects_%d_%d_%s_%s',subjects(1),subjects(end),conds_name,task_distance_name));
        if ~isequal(task_distance,task_RT)
            filename = sprintf('%s_cross_task',filename);
        end
        filename = sprintf('%s.mat',filename);
        if exist(filename,'file')
            load(filename,'stats_dth_sl');
        else
            stats_dth_sl.num_perms = 10000;
            stats_dth_sl.tail = 'left';
            stats_dth_sl.qvalue = 0.05;
            [stats_dth_sl.significant_channels,stats_dth_sl.pvalues,...
                stats_dth_sl.crit_p, stats_dth_sl.adjusted_pvalues]...
                = fdr_permutation_cluster_1sample_alld(for_stats_data,...
                stats_dth_sl.num_perms,stats_dth_sl.tail,stats_dth_sl.qvalue);
            save(filename,'stats_dth_sl');
        end   
    end
    
    %% Plot searchlight results at peak 
    searchlight_patterns = topo.(obj);
    
    %Color map
    clim = [-0.2,0.2];
    caxis(clim);
    warning('off');
    color_upper = cbrewer2('div',color_scheme, 50); 
    colors =  color_upper;

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
    topoplot(searchlight_patterns, c, 'colormap', eval([func_color '(colors)']), 'style','both','electrodes','labels','whitebk','on');

    %Color bar
    CBar_Handle = colorbar('West');
    caxis(clim);
    set(get(CBar_Handle, 'YLabel'), 'String', 'Decoding accuracy (%)-50',...
        'FontSize', 10, 'FontName', 'Arial');
    set(CBar_Handle,'Location','eastoutside');

    %Plotting parameters
    title([num2str(peak_time),' ms']);
    set(gca,'visible', 'off');
    set(gcf, 'color','white');
    fig_width = 600;
    fig_height = 300;
    set(gcf,'position',[100 100 fig_width fig_height]);
    set(gcf,'Renderer','painters');

    file_name = 'dth_searchlight';
    if ~isequal(task_distance,task_RT)
        file_name = sprintf('%s_cross_task',file_name);
    end
    
    %Save
    keyboard;
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%s_%s_subjects_%d_%d_searchlight_peak',file_name,task_distance_name,conds_name,subjects(1),subjects(end)))); %save as matlab figure
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%s_%s_subjects_%d_%d_searchlight_peak.svg',file_name,task_distance_name,conds_name,subjects(1),subjects(end)))); %save as svg
    close(gcf);    
    
end

end 
