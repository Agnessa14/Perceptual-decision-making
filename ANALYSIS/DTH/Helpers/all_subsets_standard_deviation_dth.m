function all_subsets_standard_deviation_dth(subjects,subsets,category) %distance art, distance nat, RT
%ALL_SUBSETS_STANDARD_DEVIATION_DTH Plot the average baseline standard deviation in the
% of the DTH results, in different subsets of data.
%
%Input: subjects' ID (e.g., 1:13), subsets (e.g., [0.1 0.2])
%
%Returns a plot of subset versus baseline standard deviation. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS';
addpath(genpath(results_dir));
all_subsets = '_';
for s = 1:numel(subsets)
    all_subsets = strcat(all_subsets,num2str(subsets(s)*100),'_');
end
save_path = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';

%% Get the correlations from all subjects and subsets
%Preallocate
numTimepoints = 200;
correlations_all_subjects = NaN(numel(subjects),numTimepoints);
std_dth = NaN(numel(subsets),numTimepoints);

%Loop over subsets 
for s = 1:numel(subsets)
    for subject = subjects
        correlations_all_subjects(subject,:) = core_distance_to_hyperplane(subject,subsets(s),category,numTimepoints);    
        plot(correlations_all_subjects(subject,:),'LineWidth',2);
        hold on;
    end    
    %Plot avg curve
    plot(mean(correlations_all_subjects,1),'--','LineWidth',4);
    set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
    title(sprintf('%s percent of the data', num2str(subsets(s)*100)));
    xlabel('Timepoints')
    ylabel('Spearman''s coefficient')
    saveas(gcf,fullfile(save_path,sprintf('not_fixed_effect_dth_category_%s_subset_%s',category, num2str(subsets(s)*100))));
    saveas(gcf,fullfile(save_path,sprintf('not_fixed_effect_dth_category_%s_subset_%s.svg',category, num2str(subsets(s)*100))));
    close(gcf);
    
    %Get SD
    std_dth(s,:) = std(correlations_all_subjects);    %input is numSubjects x numTimepoints
end

%% Average over baseline timepoints
std_dth_baseline_avg = squeeze(mean(std_dth(:,1:40),2));

%Setup figure
subset_values = subsets*100;
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
plot(subset_values,std_dth_baseline_avg,'--','LineWidth',2);

% %Plot parameters
title(sprintf('Standard deviation in the baseline of the correlation between the distance-to-hyperplane and RT over different subsets of data  (N=%d)', numel(subjects)));
xticks(subset_values);
xlabel('Subset of data (%)')
ylabel('Standard deviation')

%% Save
saveas(gcf,fullfile(save_path,sprintf('standard_deviation_dth_category_%s_subsets_%s',category, all_subsets)));
saveas(gcf,fullfile(save_path,sprintf('standard_deviation_dth_category_%s_subsets_%s.svg',category, all_subsets)));
end

%     correlations_all_subjects = arrayfun(@ (x) core_distance_to_hyperplane(x,subsets(s),'artificial',numTimepoints), subjects, 'UniformOutput','false'); 