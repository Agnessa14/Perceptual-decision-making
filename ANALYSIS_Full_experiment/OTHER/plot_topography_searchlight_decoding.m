function plot_topography_searchlight_decoding(subjects,analysis,with_stats)
%PLOT_TOPOGRAPHY_SEARCHLIGHT_DECODING Plot the decoding searchlight results
%at peak decoding accuracy. 
%
%Input: subject IDs, analysis ('object_decoding' or 'category_decoding',
%with or without stats (1/0))
%
%Output: one topography plot
%

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
close; %close EEGLAB popup
% load channel location
load(fullfile(main_path,'OTHER/adult_63channels.mat'),'channelloc'); 

%% Load Results
%Preallocate
numConditions = 60;
numChannels = 63; 

if strcmp(analysis,'object_decoding')
    decoding_accuracies_all_subjects_cat = NaN(max(subjects),numConditions,numConditions,numChannels);
    decoding_accuracies_all_subjects_dis = NaN(max(subjects),numConditions,numConditions,numChannels);
elseif strcmp(analysis,'category_decoding')
    decoding_accuracies_all_subjects_cat = NaN(max(subjects),numChannels);
    decoding_accuracies_all_subjects_dis = NaN(max(subjects),numChannels);
end

%Loop: collect results from all subjects
for subject = subjects
    subname = get_subject_name(subject);
    
    subject_results_dir = fullfile(results_dir,subname);
    if strcmp(analysis,'object_decoding')
        cat_filename = 'svm_searchlight_object_decoding_accuracy_peak_categorization.mat';
        dis_filename = 'svm_searchlight_object_decoding_accuracy_peak_fixation.mat';
        all_dimensions = repmat({':'},1,3); %conditions x conditions x channels
    elseif strcmp(analysis,'category_decoding')
        cat_filename = 'cross_validated_dth_pseudotrials_svm_decodingAccuracy_searchlight_peak_categorization.mat';
        dis_filename = 'cross_validated_dth_pseudotrials_svm_decodingAccuracy_searchlight_peak_fixation.mat';
        all_dimensions = {':'}; % channels
    end
    if strcmp(analysis,'object_decoding')
        load(fullfile(subject_results_dir,cat_filename),'decodingAccuracy_avg');
        decoding_accuracies_all_subjects_cat(subject,all_dimensions{:}) = decodingAccuracy_avg;
        load(fullfile(subject_results_dir,dis_filename),'decodingAccuracy_avg');
        decoding_accuracies_all_subjects_dis(subject,all_dimensions{:}) = decodingAccuracy_avg;
    elseif strcmp(analysis,'category_decoding')
        load(fullfile(subject_results_dir,cat_filename),'decodingAccuracy_avg');
        decoding_accuracies_all_subjects_cat(subject,all_dimensions{:}) = decodingAccuracy_avg(3,:);%at peak decoding
        load(fullfile(subject_results_dir,dis_filename),'decodingAccuracy_avg');
        decoding_accuracies_all_subjects_dis(subject,all_dimensions{:}) = decodingAccuracy_avg(4,:);
    end

end

%% For stats matrix: num subjects x num timepoints (average over conditions)
if strcmp(analysis,'object_decoding')
    for_stats_cat = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects_cat(subjects,:,:,:),2),3));
    for_stats_dis = squeeze(nanmean(nanmean(decoding_accuracies_all_subjects_dis(subjects,:,:,:),2),3));
elseif strcmp(analysis,'category_decoding')
    for_stats_cat = decoding_accuracies_all_subjects_cat(subjects,:);
    for_stats_dis = decoding_accuracies_all_subjects_dis(subjects,:);
end

%% Average over subjects 
avg_over_conditions_all_subjects_cat = squeeze(nanmean(for_stats_cat,1));
avg_over_conditions_all_subjects_dis = squeeze(nanmean(for_stats_dis,1));

%% Stats
if with_stats
    for task = 1:2
        if task == 1
            task_name = 'categorization';    
            for_stats_data = for_stats_cat-50;
            topo.(obj) = avg_over_conditions_all_subjects_cat;
        elseif task == 2
            task_name = 'fixation';   
            for_stats_data = for_stats_dis-50;
            topo.(obj) = avg_over_conditions_all_subjects_dis;
        end

        %Stat parameters
        filename = fullfile(results_avg_dir,...
            sprintf('stats_%s_searchlight_%s_subjects_%d_%d.mat',analysis,task_name,subjects(1),subjects(end)));
        if exist(filename,'file')
            load(filename,'stats_decoding_sl');
        else
            stats_decoding_sl.num_perms = 10000;
            stats_decoding_sl.tail = 'right';
            stats_decoding_sl.qvalue = 0.01;
            [stats_decoding_sl.significant_channels,stats_decoding_sl.pvalues,...
                stats_decoding_sl.crit_p, stats_decoding_sl.adjusted_pvalues]...
                = fdr_permutation_cluster_1sample_alld(for_stats_data,...
                stats_decoding_sl.num_perms,stats_decoding_sl.tail,stats_decoding_sl.qvalue);
            save(filename,'stats_decoding_sl');
        end   

        %% Plot searchlight results at peak decoding
        searchlight_patterns = topo.(obj)-50;

        %Set non-significant channels to off
        c=channelloc;
        if with_stats
            topomask = stats_decoding_sl.significant_channels;
            for imask = 1:length(topomask)
                if topomask(imask)==0
                    c(imask).labels = ' ';
                else
                    c(imask).labels = 'O';
                end
            end
        end
        
        colors = 'cool';
        clim = [0 30];
        
        %Color bar - plot separately
        figure;
        CBar_Handle = colorbar('eastoutside');
        colormap(colors);
        caxis(clim);
        set(get(CBar_Handle, 'YLabel'), 'String', 'Decoding accuracy(%)-50',...
            'FontSize', 10, 'FontName', 'Arial');

        %Save image for the colorbar 
        keyboard; %make fullscreen
        filename_colorbar = fullfile(results_avg_dir,'colorbar_decoding_searchlight_peak.svg');
        if ~exist(filename_colorbar,'file')
            saveas(gcf,filename_colorbar); %save as svg
        end
        close(gcf);
        
        %Plot using Monika Graumann/Ben Ehniger's script for a good figure
        figure;
        topo_init = gca;
        [~,cdata]=topoplot(searchlight_patterns, c, 'style','map','electrodes','labels','whitebk','on');
        colormap(colors);
        topo_image = axes('Position',get(topo_init,'Position')); %from Monika/Ben Ehnigers - somehow works better for making figures
        uistack(topo_image,'bottom')

        %continue without colorbar because it screws everything up 
        delete(findobj(topo_init,'Type','surface'));
        h = imagesc(topo_image,cdata);
        set(topo_image,'YDir','normal');
        axis(topo_image,'square','off');
        set(h,'alphadata',~isnan(cdata));        
        caxis(clim);
        
        %peak times
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

        %Plotting parameters
        set(gca,'visible', 'off');
        set(gcf, 'color','white');
        set(gcf,'Renderer','painters');

        %Save
        keyboard; %increase figure to fullscreen %sometimes issue with the saveas function - simply F10 from here for the save lines, then F5
        save(fullfile(results_avg_dir,sprintf('for_stats_cat_svm_%s_subjects_%d_%d_searchlight_peak_%s.mat',analysis,subjects(1),subjects(end),task_name)),'for_stats_cat'); 
        saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_searchlight_peak_%s',analysis,subjects(1),subjects(end),task_name))); %save as matlab figure
        saveas(gcf,fullfile(results_avg_dir,sprintf('svm_%s_subjects_%d_%d_searchlight_peak_%s.svg',analysis,subjects(1),subjects(end),task_name))); %save as svg
        close(gcf);    
    end
end 
close all;
end