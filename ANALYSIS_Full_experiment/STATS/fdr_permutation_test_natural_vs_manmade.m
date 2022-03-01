function fdr_permutation_test_natural_vs_manmade(analysis_type)
%FDR_PERMUTATION_TEST_NATURAL_VS_MANMADE Perform fdr correction stats to calculate the
%significance of the difference between natural and man-made scenes (neural or behavioural). 
%Random effects: the subject-level data are randomly multiplied by 1 or -1 
%(sign permutation test) to create the permutation samples.
%
%Input: 
% - analysis type ('rsa' or 'rt')
%
%Output: For RSA, 1xP vector of significance (1) or not (0) per layer and RNN timestep, where P is the number
%of timepoints, along with the pvalues and the critical p-value; for RT,
%one significance value.
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));

                 %%%%% CALCULATING THE GROUND TRUTH AND PERMUTATION SAMPLES P-VALUES %%%%%
%% Load data 
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
if strcmp(analysis_type,'rsa')
    filename = 'all_subjects_all_tps_rsa_rnn';
    load(fullfile(results_avg_dir,sprintf('%s_artificial',filename)),'rsa_all_subs'); 
    data_artificial_scenes = rsa_all_subs;
    load(fullfile(results_avg_dir,sprintf('%s_natural',filename)),'rsa_all_subs');
    data_natural_scenes = rsa_all_subs;
elseif strcmp(analysis_type,'rt')
    filename = 'correlation_RT_human_RNN_cross-validated.mat';
    load(fullfile(results_avg_dir,'02.11_2_rnn/Model_RDM_redone',filename),'data');
    data_artificial_scenes = data(:,2);
    data_natural_scenes = data(:,3);
end

%% 1) Sign permutation test: randomly multiply by 1/-1 the subject-level artificial scenes' data and calculate the artificial-natural difference
size_data = size(data_artificial_scenes,2:numel(size(data_artificial_scenes)));
cln = repmat({':'},size_data); %to select the N-1 dimensions
num_subjects = size(data_artificial_scenes,1);
numPermutations = 10000;
q_value = 0.05;

%calculate ground difference
samples_plus_ground_difference = NaN([numPermutations,size_data]);
samples_plus_ground_difference(1,cln{:}) = abs(mean(data_artificial_scenes,1)-mean(data_natural_scenes,1));

for perm = 2:numPermutations    
    if ~mod(perm,100)
        fprintf('Permutation sample %d \n',perm);
    end   
    
    %create samples by randomly multiplying each subject's data by 1 or -1
    random_vector = single(sign(rand(num_subjects,1)-0.5));
    sample_art = repmat(random_vector,[1,size_data]).*data_artificial_scenes;
    sample_nat = repmat(random_vector,[1,size_data]).*data_natural_scenes;

    %get the difference of each sample
%     samples_plus_ground_difference(perm,cln{:}) = abs(mean(sample,1)-nanmean(data_natural_scenes,1));
    samples_plus_ground_difference(perm,cln{:}) = abs(mean(sample_art,1)-mean(sample_nat,1));
    
end

%% 2) Calculate the p-value of the ground truth and of the permuted samples
%right-tailed
p_ground_and_samples = (numPermutations+1 - tiedrank(samples_plus_ground_difference)) / numPermutations;

%% 3) Perform FDR correction
pvalues = squeeze(p_ground_and_samples(1,cln{:}));
[stats_natural_vs_manmade.SignificantVariables,stats_natural_vs_manmade.crit_p,~,stats_natural_vs_manmade.adjusted_pvalues] = fdr_bh(pvalues,q_value,'pdep');
save(fullfile(results_avg_dir,sprintf('stats_natural_vs_manmade_%s',analysis_type)),'stats_natural_vs_manmade');
end  
                   