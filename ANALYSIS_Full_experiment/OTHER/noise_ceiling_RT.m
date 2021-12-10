function noise_ceiling_RT(subjects) 
%NOISE_CEILING_RT Calculate the noise ceiling (lower) for the human reaction times
%of object decoding (categorization task). 
%
%Input: subject IDs (eg., 1:13), conditions ('all','artificial' or
%'natural')
%
%Output: 
% 1)Three 1xS noise ceiling vectors: all scenes, natural scenes and artificial scenes
% (S = #scenes)
%
%
%Author: Agnessa Karapetian, 2021
%

%% Paths
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load the RTs
num_subjects = numel(subjects); %included subjects are only half (the rest were used for cross validation)
filename = 'RT_all_subjects_5_35_categorization.mat';
load(fullfile(results_avg_dir,filename),'RTs');
% RTs_shuffled = RTs(num_subjects+1:end,:);

%% Average over N-1 subjects and correlate with the left-out subject 
lower_noise_ceiling = NaN(num_subjects,3);

for s = 1:num_subjects
    for c = 1:3
        if c == 1
            conds = 1:30;
        elseif c == 2
            conds = 31:60;
        elseif c == 3
            conds = 1:60;
        end
        leftout_subject = squeeze(RTs(s,conds));
        other_subjects = squeeze(nanmean(RTs(1:end~=s,conds),1)); %avg over remaining subjects
        nan_scenes = find(isnan(leftout_subject));
        if s==6
            keyboard;
        end
        if numel(nan_scenes) > 0
            leftout_subject(nan_scenes) = [];
            other_subjects(nan_scenes) = [];
        end
        lower_noise_ceiling(s,c) = corr(leftout_subject',other_subjects','type','Pearson');
    end
end

%% Average over iterations and save
noise_ceiling_lower_bound = mean(lower_noise_ceiling,1); %avg over subjects
save(fullfile(results_avg_dir,sprintf('sf_lower_bound_noise_ceiling_%s',filename)),'noise_ceiling_lower_bound');

end
    
