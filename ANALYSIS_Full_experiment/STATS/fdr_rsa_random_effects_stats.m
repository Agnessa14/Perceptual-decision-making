function [SignificantVariables, crit_p, adjusted_pvalues] = fdr_rsa_random_effects_stats(true_rsa_results,numPermutations,tail,q_value)
%FDR_RSA_RANDOM_EFFECTS_STATS Perform fdr correction stats to calculate the
%significance of the timepoints in the RSA. Random effects: the
%subject-level data are randomly multiplied by 1 or -1 (sign permutation test) to create the permutation samples.
%
%Input: 
% - true_rsa_results: the true correlation, for each subject, between rdm1&rdm2 (num subjects x
% num timepoints)
% - number of permutations for the stats (ex:10000)
% - tail: 'both', 'right' or 'left' (default right)
% - q-value: for fdr (default: 0.05)
%
%Output: 1xP vector of significance (1) or not (0), where P is the number
%of timepoints, along with the pvalues and the critical p-value.
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));

                 %%%%% CALCULATING THE GROUND TRUTH AND PERMUTATION SAMPLES P-VALUES %%%%%

%% 1) Sign permutation test: randomly multiply by 1/-1 the subject-level RSA results and calculate the t-statistic
numTimepoints = size(true_rsa_results,2);
samples_plus_ground_tstatistic = NaN(numPermutations,numTimepoints);

%first perm is ground test statistic
samples_plus_ground_tstatistic(1,:) = mean(true_rsa_results,1) ./ std(true_rsa_results); 

for perm = 2:numPermutations    
    if ~mod(perm,100)
        fprintf('Permutation sample %d \n',perm);
    end   
    
    %create samples by randomly multiplying each subject's data by 1 or -1
    random_vector = single(sign(rand(30,1)-0.5));
    sample = repmat(random_vector,1,numTimepoints).*true_rsa_results;
    
    %get the  test statistic (mean ./ std) of each sample
    samples_plus_ground_tstatistic(perm,:) = mean(sample,1) ./ std(sample);
    
end

%% 2) Calculate the p-value of the ground truth and of the permuted samples
if strcmp(tail,'right')
    p_ground_and_samples = (numPermutations+1 - tiedrank(samples_plus_ground_tstatistic)) / numPermutations;
else
    error('Wrong tail');
end 

%% 3) Perform FDR correction
pvalues = squeeze(p_ground_and_samples(1,:));
[SignificantVariables,crit_p,~,adjusted_pvalues] = fdr_bh(pvalues,q_value,'pdep');

end  
                   