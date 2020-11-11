function subcategories_svm_fixed_analysis_pilot_2(subjects) %distance art, distance nat, RT
%SUBCATEGORIES_SVM_FIXED_ANALYSIS_PILOT_2 Performs the distance-to-hyperplane analysis using
%the svm classifier on the subcategories (forest, beach, canyon, apartment, bedroom, highway).
%
%Input: subjects' ID (e.g., 1:13)
%
%Correlates the decision values with reaction times (averaged over
%participants) of each condition (60 scenes), at each timepoint, resulting in a plot of Spearman's correlation vs time. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS';
addpath(genpath(results_dir));
save_path = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';

%% Get the distances from all subjects
numTimepoints = 200;
numCategories = 2;
numSubcategories = 6;
numSubcategoriesPerCategory = numSubcategories/numCategories;
numConditions = 60;
numSubConditions = numConditions/numSubcategories;
distances = NaN(numel(subjects),numConditions,numTimepoints);
RTs = NaN(numel(subjects),numConditions);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,'true_pseudotrials_svm_decisionValues.mat'));
    distances(subject,:,:) = decisionValues_Avg;   
    load(fullfile(results_dir,subname,'RTs_correct_trials.mat'));
    RTs(subject,:) = normalize(RT_per_condition);
end

medianRT = nanmedian(RTs,1);

%% Analysis within each subcondition 
subcategory_names = {'apartment','bedroom','highway'; 'beach','canyon','forest'};
category_names = {'artificial','natural'};
timepoints = 1:numTimepoints;

for c = 1:numCategories
    
    figure(abs(round(randn*10)));
    set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen
    legend_cell = cell(numSubcategoriesPerCategory,1);

    for s = 1:numSubcategoriesPerCategory
        subcat_name = subcategory_names{c,s};
        if c == 1  
            scenes = ((s-1)*numSubConditions)+1:s*numSubConditions; 
        elseif c == 2
            scenes = ((s+2)*numSubConditions)+1:(s+3)*numSubConditions; 
        end
        disp(scenes);
        
        RT = medianRT(scenes);
        mean_distances = squeeze(nanmean(distances(:,scenes,:),1)); %avg over subjects
        correlation_dth_RT = arrayfun(@ (x) corr(RT',mean_distances(:,x),'type','Spearman'),timepoints);
        legend_cell{s} = subcat_name;
        plot(correlation_dth_RT,'LineWidth',2);
        hold on;
        save(fullfile(save_path,sprintf('%s_SVM_DTH_rt_correlation_%d_subjects.mat',...
            subcat_name,numel(subjects))),'correlation_dth_RT');
    end
 
    %% Plotting parameters
    xline(40,'--'); 
    title(sprintf('Correlation between distance to hyperplane and reaction time in 3 subordinate %s categories',category_names{c}));
    xlabel('Timepoints');
    ylabel('Spearman''s coefficient');
    xticks(0:10:200);
    legend(legend_cell,'FontSize',12);
%     set(leg, 'position', [0.7 0.2 0.1 0.01]); %put legend below graph

    %% Save figures
    saveas(gcf,fullfile(save_path,sprintf('subcategories_%s_SVM_DTH_artificial_natural_%d_subjects',category_names{c},numel(subjects)))); 
    saveas(gcf,fullfile(save_path,sprintf('subcategories_%s_SVM_DTH_artificial_natural_%d_subjects.svg',category_names{c},numel(subjects))));
    close(gcf);

end

end


% apartment = 1:10;
% bedroom   = 11:20;
% highway   = 21:30;
% beach     = 31:40;
% canyon    = 41:50;
% forest    = 51:60;