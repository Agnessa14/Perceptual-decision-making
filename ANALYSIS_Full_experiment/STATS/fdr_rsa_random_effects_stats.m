function [SignificantVariables, pvalues, crit_p, adjusted_pvalues] = fdr_rsa_random_effects_stats(rdm_1,rdm_2,true_rsa_rdm,numPermutations,tail,q_value)
%FDR_RSA_RANDOM_EFFECTS_STATS Perform fdr correction stats to calculate the
%significance of the timepoints in the RSA. Random effects: the
%subject-level data are randomly multiplied by 1 or -1 (sign permutation test) to create the permutation samples.
%
%Input: 
% - rdm 1 & rdm 2: the data that goes into RSA (num subjects x num pairwise
% combinations x num timepoints)
% - true_rsa_rdm: the true correlation between rdm1&rdm2
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

%% 1) Randomly multiply by 1/-1 the subject-level RDMs and calculate the Spearman's correlation with the other RDM at each timepoint
numTimepoints = size(true_rsa_rdm,1);
all_rsa_rdm = NaN(max(subjects),numPermutations,numTimepoints);
for perm = 2:numPermutations
    if ~mod(perm,100)
        fprintf('Calculating the correlation %d \n',perm);
    end
    
    random_vector = single(sign(rand(30,1)-0.5));
    random_data = random_vector*rdm_1;
    for subject = subjects
        all_rsa_rdm(subject,1,:,:) = true_rsa_rdm;
        all_rsa_rdm(subject,perm,:,:) = representational_SA_rnn(random_data,rdm_2);   
    end
    all_rsa_rdm = all_rsa_rdm(subjects,:,:,:);
    avg_rsa = squeeze(mean(all_rsa_rdm,1));
end

%% 2) Calculate the p-value of the ground truth and of the permuted samples
if strcmp(tail,'right')
    p_ground_and_samples = (numPermutations+1 - tiedrank(avg_rsa)) / numPermutations;
else
    error('Wrong tail');
end 

%% 3) Perform FDR correction
pvalues = squeeze(p_ground_and_samples(1,:,:));
[SignificantVariables,crit_p,~,adjusted_pvalues] = fdr_bh(pvalues,q_value,'pdep');

end  
                   