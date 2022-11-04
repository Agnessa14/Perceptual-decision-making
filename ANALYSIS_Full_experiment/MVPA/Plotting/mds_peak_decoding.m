function mds_peak_decoding(subjects,task,analysis,filetype)
%MDS_PEAK_DECODING Create MDS plot from the object decoding data at
%peak decoding, and as a video across all timepoints
%
%Input: subject IDs, task (1=categorization,2=distraction), analysis
%('object_decoding' or 'category_decoding'),filetype ('image' or 'video')
%
%Output: MDS plot
%
%Author: Agnessa Karapetian, 2022

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
task_name = get_task_name(task);

%% Preallocate
numConditions = 60;
numTimepoints = 200;
conds_art = 1:numConditions/2;
conds_nat = (numConditions/2)+1:numConditions;
decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints);

%% Loop: collect results from all subjects
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)),...
            'decodingAccuracy_avg');
    decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;  
end   

%% Remove any NaN (for non-included subjects)
avg_over_conditions_all_subjects = squeeze(nanmean(decoding_accuracies_all_subjects,1));

%% Plotting parameters
dotsize = 150;
cmap_1 = cool;
cmap_2 = summer;
color_art = cmap_1(200,:); %purple
color_nat = cmap_2(100,:); %green

if strcmp(filetype,'image')

    %% 1) Plot at peak time
    if task == 1 && strcmp(analysis,'object_decoding')
        peak_time = 120;
    elseif task == 1 && strcmp(analysis,'category_decoding')
        peak_time = 160;
    elseif task == 2 && strcmp(analysis,'object_decoding')
        peak_time = 110;
    elseif task == 2 && strcmp(analysis,'category_decoding')
        peak_time = 155;
    end
    peak_time_tp = (peak_time)/5+40;
    data_peak = squeeze(avg_over_conditions_all_subjects(:,:,peak_time_tp));

    %flatten into a vector
    data_peak(isnan(data_peak)) = 0;
    data_peak_symm = data_peak+data_peak';

    %mds
    mds = cmdscale(data_peak_symm,2);
    figure(abs(round(randn*10))); %Random figure number
    points_art_peak = squeeze(mds(conds_art,:));
    points_nat_peak = squeeze(mds(conds_nat,:));
    scatter(points_art_peak(:,1),points_art_peak(:,2),dotsize,color_art,'filled');
    hold on;
    scatter(points_nat_peak(:,1),points_nat_peak(:,2),dotsize,color_nat,'filled');

    %% Save the plot
    filename = sprintf('mds_dots_peak_%s_%s',analysis,task_name);
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d',filename,subjects(1),subjects(end)))); %save as matlab figure
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d.svg',filename,subjects(1),subjects(end)))); %save as svg    
    close(gcf);

    %% Plot with images
    figure(abs(round(randn*10))); %Random figure number
    for c = 1:numConditions
        point = squeeze(mds(c,:));
        image_name = fullfile('/home/agnek95/stim_cropped',sprintf('%d.jpg',c));
        M = imread(image_name);
        imagesc('XData',[point(1)-5 point(1)+5] ,'YData', [point(2)+5 point(2)-5 ],'CData', M);
        daspect([1 1 1]); hold on
        drawnow
    end

    filename = sprintf('mds_images_peak_%s_%s',analysis,task_name);
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d',filename,subjects(1),subjects(end)))); %save as matlab figure
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d.svg',filename,subjects(1),subjects(end)))); %save as svg    
    close(gcf);
    
elseif strcmp(filetype,'video')
    
    %% 2) Plot over the whole time course and record video
    figure(abs(round(randn*10))); %Random figure number
    keyboard; %here start screen recording
    check_with_nontr = 0; %plot the procr transformed vs non transformed
    procrustes_points_art=NaN(numConditions/2,2,numTimepoints);
    procrustes_points_nat=NaN(numConditions/2,2,numTimepoints);
    points_art=NaN(numConditions/2,2,numTimepoints);
    points_nat=NaN(numConditions/2,2,numTimepoints);
    
    %Do procrustes transform
    for t = 1:numTimepoints
        data_tp = squeeze(avg_over_conditions_all_subjects(:,:,t));

        %flatten into a vector
        data_tp(isnan(data_tp)) = 0;
        data_tp_symm = data_tp+data_tp';

        %mds
        mds = cmdscale(data_tp_symm,2);    
        points_art(:,:,t) = mds(conds_art,:);
        points_nat(:,:,t) = mds(conds_nat,:);
        if t==1
            procrustes_points_art(:,:,t) = points_art(:,:,t); %no transform needed
            procrustes_points_nat(:,:,t) = points_nat(:,:,t);
        else
            [~,procrustes_points_art(:,:,t)] = procrustes(points_art(:,:,t-1),points_art(:,:,t));
            [~,procrustes_points_nat(:,:,t)] = procrustes(points_nat(:,:,t-1),points_nat(:,:,t));
        end
    end
    
    %Plot
    for t = 1:numTimepoints  
        if check_with_nontr == 1          
            scatter(points_art(:,1,t),points_art(:,2,t),dotsize,'b','filled'); %only way that worked to plot in diff colours
            hold on;
            scatter(points_nat(:,1,t),points_nat(:,2,t),dotsize,'k','filled');  
        end
        scatter(procrustes_points_art(:,1,t),procrustes_points_art(:,2,t),dotsize,color_art,'filled'); %only way that worked to plot in diff colours
        hold on;
        scatter(procrustes_points_nat(:,1,t),procrustes_points_nat(:,2,t),dotsize,color_nat,'filled');    
        title(sprintf('%d ms',(t-40)*5));
        legend({'Man-made scenes','Natural scenes'},'FontSize',12,'Location','southoutside');
        axis off;
        hold off;
        pause(0.1);
%         saveas(gcf,fullfile(results_avg_dir,'MDS',sprintf('mds_procrustes_timepoint_%d_%s_%s',t,analysis,task_name)));
    end
    close;
end
end