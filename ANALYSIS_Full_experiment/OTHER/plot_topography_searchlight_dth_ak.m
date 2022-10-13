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
close; %close EEGLAB popup

% load channel location
load(fullfile(main_path,'OTHER/adult_63channels.mat'),'channelloc'); 

%% Load Results
% numConditions = 60;
% numChannels = 63; 
task_distance_name = get_task_name(task_distance);
file_name = 'dth_searchlight_peak';
if ~isequal(task_distance,task_RT)
    file_name = sprintf('%s_cross_task',file_name);
end
file_name_full = sprintf('%s_subjects_%d_%d_%s.mat',file_name,subjects(1),subjects(end),task_distance_name);
load(fullfile(results_avg_dir,file_name_full),'dth_results');

%% Stats
for conditions = 3
    if conditions == 1
        conds_name = 'artificial';    
        for_stats_data = dth_results.for_stats_corr_artificial(subjects,:);
        topo.(obj) = dth_results.corr_artificial;

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
            peak_time = 160;
        elseif task_distance==1 && task_RT==2
            peak_time = 110;
        elseif task_distance==2 && task_RT==2
            peak_time = 110;
        elseif task_distance==2 && task_RT==1
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
        if exist('filename','file')
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
    colors = cbrewer2('div','PiYG', 50); 


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

    %Color bar - plot separately
    figure;
    CBar_Handle = colorbar('eastoutside');
    colormap(colors);
    caxis(clim);
    set(get(CBar_Handle, 'YLabel'), 'String', 'Spearman''s p',...
        'FontSize', 10, 'FontName', 'Arial');

    %Save image for the colorbar 
    keyboard; %make fullscreen
    filename_colorbar = fullfile(results_avg_dir,'colorbar_dth_searchlight_peak.svg');
    if ~exist(filename_colorbar,'file')
        saveas(gcf,filename_colorbar); %save as svg
    end
    close(gcf);
    
    %Plot - Monika/Ben Ehniger's code for good figures
    figure;
    topo_init = gca;
    [~,cdata]=topoplot(searchlight_patterns,c,'colormap',colors,'style','map','electrodes','labels','whitebk','on');
    topo_image = axes('Position',get(topo_init,'Position')); %from Monika/Ben Ehnigers - somehow works better for making figures
    uistack(topo_image,'bottom')

    %continue without colorbar because it screws everything up 
    delete(findobj(topo_init,'Type','surface'));
    h = imagesc(topo_image,cdata);
    set(topo_image,'YDir','normal');
    axis(topo_image,'square','off');
    set(h,'alphadata',~isnan(cdata));
    caxis(clim);
    
    %Plotting parameters
    title([num2str(peak_time),' ms']);
    set(gca,'visible', 'off');
    set(gcf, 'color','white');
    set(gcf,'Renderer','painters');

    %Save
    file_name = 'dth_searchlight';
    if ~isequal(task_distance,task_RT)
        file_name = sprintf('%s_cross_task',file_name);
    end
    
    keyboard;
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%s_%s_subjects_%d_%d_searchlight_peak',file_name,task_distance_name,conds_name,subjects(1),subjects(end)))); %save as matlab figure
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%s_%s_subjects_%d_%d_searchlight_peak.svg',file_name,task_distance_name,conds_name,subjects(1),subjects(end)))); %save as svg
    close(gcf);    
    
end
close;

end 
