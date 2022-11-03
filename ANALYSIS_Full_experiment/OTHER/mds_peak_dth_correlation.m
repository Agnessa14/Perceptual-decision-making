function mds_peak_dth_correlation(subjects,task)
%MDS_PEAK_DTH_CORRELATION Create MDS plot from the object decoding data at
%peak negative correlation with distance-to-hyperplane.
%
%Input: subject IDs, task (1=categorization,2=distraction)
%
%Output: MDS plot
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
task_name = get_task_name(task);

%% Preallocate
numConditions = 60;
numTimepoints = 200;
decoding_accuracies_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints);

%% Loop: collect results from all subjects + plot each subject individually on the same plot
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,sprintf('svm_decoding_accuracy_%s.mat',task_name)),...
            'decodingAccuracy_avg');
    decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;  
end   

%% Remove any NaN (for non-included subjects)
avg_over_conditions_all_subjects = squeeze(nanmean(decoding_accuracies_all_subjects,1));

%% Load the DTH curve to check when the correlation is the most negative
load(fullfile(results_avg_dir,sprintf('cv_all_dist_med_rt_dth_subjects_%d_%d_%s.mat',subjects(1),subjects(end),task_name)),'dth_results');
min_correlation = min(dth_results.corr_both_categories);
peak_time = dth_results.corr_both_categories==min_correlation;
data_peak = squeeze(avg_over_conditions_all_subjects(:,:,peak_time));

%flatten into a vector
data_peak(isnan(data_peak)) = 0;
data_peak_symm = data_peak+data_peak';

%mds
mds = cmdscale(data_peak_symm,2);
cmap_1 = cool;
cmap_2 = summer;
color_art = cmap_1(200,:); %purple
color_nat = cmap_2(100,:); %green

figure(abs(round(randn*10))); %Random figure number
for c = 1:numConditions
    point = squeeze(mds(c,:));
    if c<31
        colorr = color_art;
    else
        colorr = color_nat;
    end
    scatter(point(1),point(2),150,colorr,'filled');
    hold on;
end

% %% Save the plot
% filename = 'mds_object_decoding_subjects';
% saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s',filename,subjects(1),subjects(end),task_name))); %save as matlab figure
% saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s.svg',filename,subjects(1),subjects(end),task_name))); %save as svg    
% close(gcf);

%% Plot with images
figure(abs(round(randn*10))); %Random figure number
for c = 1:numConditions
    point = squeeze(mds(c,:));
    image_name = fullfile('stim_cropped',sprintf('%d.jpg',c));
    M = imread(image_name);
    imagesc('XData',[point(1)-5 point(1)+5] ,'YData', [point(2)+5 point(2)-5 ],'CData', M);
    daspect([1 1 1]); hold on
    drawnow
end

filename = 'images_mds_object_decoding_subjects';
saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s',filename,subjects(1),subjects(end),task_name))); %save as matlab figure
saveas(gcf,fullfile(results_avg_dir,sprintf('%s_%d_%d_%s.svg',filename,subjects(1),subjects(end),task_name))); %save as svg    
close(gcf);

end