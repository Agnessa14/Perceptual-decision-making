function noise_ceiling_rsa(subjects,conditions) 
%NOISE_CEILING Calculate the noise ceiling for the representational similarity analysis
%of object decoding (categorization task). 
%
%Input: subject IDs (eg., 1:13), conditions ('all','artificial' or
%'natural')
%
%Output: 
% 1)Three 1xP noise ceiling vectors: all scenes, natural scenes and artificial scenes
% (P = # timepoints)
%
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load existing RDMs or construct them from individual ones
if strcmp(conditions,'all')
    conds = 1:60;
elseif strcmp(conditions,'artificial')
    conds = 1:30;
elseif strcmp(conditions,'natural')
    conds = 31:60;
else
    error('Wrong conditions specified');
end
numConditions = numel(conds);
numTimepoints = 200;
numSubjects = numel(subjects);
rdm_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints); 

%% Loop: collect results from all subjects 
for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'rdm_pearson_categorization.mat'),'rdm_avg');           
    rdm_all_subjects(subject,:,:,:) = rdm_avg(conds,conds,:); %average over the RDMs of each half of trials      
end   

%% Remove all subject nans 
rdm_subjects = rdm_all_subjects(subjects,:,:,:);
lower_noise_ceiling = NaN(numSubjects,numTimepoints);
upper_noise_ceiling = NaN(numSubjects,numTimepoints);

%% Average over N-1 subjects and correlate with the left-out subject 
for s = 1:numSubjects
    leftout_subject = squeeze(rdm_subjects(s,:,:,:));
    other_subjects = squeeze(nanmean(rdm_subjects(1:end~=s,:,:,:),1)); %avg over remaining subjects
    all_subjects = squeeze(nanmean(rdm_subjects,1));
    
    %take upper triangular & flatten
    leftout_subject(isnan(leftout_subject)) = 0;
    other_subjects(isnan(other_subjects)) = 0;
    all_subjects(isnan(all_subjects)) = 0;
    leftout_flattened_cell = arrayfun(@(x) squareform(leftout_subject(:,:,x)+(leftout_subject(:,:,x))'),...
        1:numTimepoints,'UniformOutput',false);
    leftout_flattened = reshape(cell2mat(leftout_flattened_cell),[],numTimepoints);
    other_flattened_cell = arrayfun(@(x) squareform(other_subjects(:,:,x)+(other_subjects(:,:,x))'),...
        1:numTimepoints,'UniformOutput',false);
    other_flattened = reshape(cell2mat(other_flattened_cell),[],numTimepoints);   
    all_flattened_cell = arrayfun(@(x) squareform(all_subjects(:,:,x)+(all_subjects(:,:,x))'),...
        1:numTimepoints,'UniformOutput',false);
    all_flattened = reshape(cell2mat(all_flattened_cell),[],numTimepoints);   

    %upper & lower bound ceilings
    for tp = 1:numTimepoints
        lower_noise_ceiling(s,tp) = corr(leftout_flattened(:,tp),other_flattened(:,tp));
        upper_noise_ceiling(s,tp) = corr(leftout_flattened(:,tp),all_flattened(:,tp));
    end
end

%% Average over iterations
noise_ceiling_lower_bound = mean(lower_noise_ceiling,1); %avg over subjects
noise_ceiling_upper_bound = mean(upper_noise_ceiling,1); %avg over subjects

%% Plot
plot(noise_ceiling_lower_bound,'DisplayName','Lower bound');
legend('-DynamicLegend');
hold all;
plot(noise_ceiling_upper_bound,'DisplayName','Upper bound');
xticks(0:20:200);
xticklabels(-200:100:800);   
ylim([-0.1 0.8]);
plot_parameters = 1;
if plot_parameters == 1
%     title(sprintf('Upper and lower bound noise ceilings for %s scenes',conditions));
    xlabel('Time (ms)');
    ylabel('Spearman''s p');
    font_size = 12;
    set(gca,'FontName','Arial','FontSize',font_size);
end
filename = sprintf('noise_ceiling_rsa_subjects_%d_%d_%s_scenes',subjects(1),subjects(end),conditions);
saveas(gcf,sprintf('%s.svg',filename));
saveas(gcf,sprintf('%s.fig',filename));
close(gcf);

%Save matrices
save(fullfile(results_avg_dir,sprintf('lower_%s',filename)),'noise_ceiling_lower_bound');
save(fullfile(results_avg_dir,sprintf('upper_%s',filename)),'noise_ceiling_upper_bound');

end
    
