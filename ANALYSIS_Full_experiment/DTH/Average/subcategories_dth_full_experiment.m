function subcategories_dth_full_experiment(subjects,task) %distance art, distance nat, RT
%SUBCATEGORIES_DTH_FULL_EXPERIMENT Performs the distance-to-hyperplane analysis using
%the svm classifier on the subcategories (forest, beach, canyon, apartment, bedroom, highway).
%
%Input: subjects' ID (e.g., 1:13), task (1=categorization, 2=distraction)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load the overall curves
task_name = get_task_name(task);
% load(fullfile(results_avg_dir,sprintf('pseudotrials_SVM_DTH_%d_subjects_%s.mat',numel(subjects),task_name)));
% load(fullfile(results_avg_dir,sprintf('true_pseudotrials_SVM_DTH_rt_correlation_natural_%d_subjects.mat',numel(subjects))));

%% Get the distances from all subjects
%define and preallocate variables
numCategories = 2;
numSubcategories = 6;
numSubcategoriesPerCategory = numSubcategories/numCategories;
numConditions = 60;
numSubConditions = numConditions/numSubcategories;
numTimepoints = 200;
distances = NaN(numel(subjects),numConditions,numTimepoints);
RTs = NaN(max(subjects),numConditions);
RTs_art = NaN(max(subjects),numConditions/2);
RTs_nat = NaN(size(RTs_art));
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;

%get the distances and RTs for all subjects
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_name)));
    distances(subject,:,:) = decisionValues_Avg;   
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name)));

    %normalize RTs
    RTs(subject,:) = normalize(RT_per_condition);
    RTs_art(subject,:) = normalize(RT_per_condition(artificial_conditions));
    RTs_nat(subject,:) = normalize(RT_per_condition(natural_conditions));
end

%% Get the median RTs and mean distances of all subjects for each condition 
% medianRT = nanmedian(RTs,1);
medianRT_art = nanmedian(RTs_art,1);
medianRT_nat = nanmedian(RTs_nat,1);
% mean_distances = squeeze(nanmean(distances,1)); %avg over subjects

%% Analysis within each subcondition 
subcategory_names = {'apartment','bedroom','highway'; 'beach','canyon','forest'};
category_names = {'artificial','natural'};
timepoints = 1:numTimepoints;

for c = 1:numCategories
    
    figure(abs(round(randn*10)));
    set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
    legend_cell = cell(numSubcategoriesPerCategory,1);
    all_scenes = NaN(numConditions/2,1);
    
    for s = 1:numSubcategoriesPerCategory
        subcat_name = subcategory_names{c,s};
        if c == 1  
            scenes = ((s-1)*numSubConditions)+1:s*numSubConditions; 
            all_scenes(scenes) = scenes;
            RT = medianRT_art(scenes);
        elseif c == 2
            scenes = ((s+2)*numSubConditions)+1:(s+3)*numSubConditions; 
            all_scenes(scenes-30) = scenes;
            RT = medianRT_nat(scenes-30);
        end

        mean_distances = squeeze(nanmean(distances(:,scenes,:),1)); %avg over subjects
        correlation_dth_RT = arrayfun(@ (x) corr(RT',mean_distances(:,x),'type','Spearman'),timepoints);
        legend_cell{s} = subcat_name;
        plot(correlation_dth_RT,'LineWidth',2);
        hold on;
        save(fullfile(results_avg_dir,sprintf('%s_SVM_DTH_rt_correlation_subjects_%d_%d_%s.mat',...
            subcat_name,subjects(1),subjects(end),task_name)),'correlation_dth_RT');
    end
 
    %% Plotting parameters
    if c == 1
        RT_all = medianRT_art(all_scenes);
    elseif c == 2
        RT_all = medianRT_nat(all_scenes-30);
    end
    mean_distances_all = squeeze(nanmean(distances(:,all_scenes,:),1));
    correlation_dth_RT = arrayfun(@ (x) corr(RT_all',mean_distances_all(:,x),'type','Spearman'),timepoints);
    plot(correlation_dth_RT,'k --','LineWidth',2);

    xline(40,'--'); 
    title(sprintf('Correlation between distance to hyperplane and reaction time in 3 subordinate %s categories',category_names{c}));
    xlabel('Timepoints');
    ylabel('Spearman''s coefficient');
    xticks(0:10:200);
    legend_cell = [legend_cell(:)',{'all'}];
    legend(legend_cell,'FontSize',12,'Location','Best');

    %% Save figures
    saveas(gcf,fullfile(results_avg_dir,sprintf('subcategories_%s_SVM_DTH_subjects_%d_%d_%s',category_names{c},subjects(1),subjects(end),task_name))); 
    saveas(gcf,fullfile(results_avg_dir,sprintf('subcategories_%s_SVM_DTH_subjects_%d_%d_%s.svg',category_names{c},subjects(1),subjects(end),task_name)));
    close(gcf);
end

end


% apartment = 1:10;
% bedroom   = 11:20;
% highway   = 21:30;
% beach     = 31:40;
% canyon    = 41:50;
% forest    = 51:60;