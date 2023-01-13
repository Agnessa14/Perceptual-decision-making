function [SignificantVariables, pvalues, crit_p, adjusted_pvalues] = fdr_rsa_rnn_eeg(rdm_1,rdm_2,true_rsa_rdm,numPermutations,tail,q_value)
%FDR_RSA_RNN_EEG Perform fdr correction stats to calculate the
%significance of the timepoints in the RSA. The EEG RDM used is the average
%over subjects.
%
%Input: 
% - rdm 1 & rdm 2: the data that goes into RSA (num conditions x num conditions x num timepoints)
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

%% 1)Permute the subject-level RDMs N times and calculate the Spearman's correlation with the other RDM at each timepoint
numTimepoints = size(rdm_1,3);
all_rsa_rdm = NaN(numPermutations,numTimepoints);
all_rsa_rdm(1,:) = true_rsa_rdm;
for perm = 2:numPermutations
    if ~mod(perm,100)
        fprintf('Calculating the correlation %d \n',perm);
    end
    %flatten and permute EEG RDM
%     rdm_1(isnan(rdm_1)) = 0;
%     rdm_flattened_cell_1 = arrayfun(@(x) squareform(rdm_1(:,:,x)+(rdm_1(:,:,x))'),...
%         1:numTimepoints,'UniformOutput',false);
%     rdm_flattened_1 = reshape(cell2mat(rdm_flattened_cell_1),[],numTimepoints);
%     permuted_rdm_1 = rdm_flattened_1(randperm(size(rdm_flattened_1,1)),:);
    random_order = randperm(size_rdm_1,1);
    permuted_rdm_1 = rdm_1(random_order,random_order,:);
    
    %RSA
    all_rsa_rdm(perm,:) = representational_SA_rnn(permuted_rdm_1,rdm_2); %modify the RSA function    
end

%% 2) Calculate the p-value of the ground truth and of the permuted samples
if strcmp(tail,'right')
    p_ground_and_samples = (numPermutations+1 - tiedrank(all_rsa_rdm)) / numPermutations;
else
    error('Wrong tail');
end 

%% 3) Perform FDR correction
pvalues = squeeze(p_ground_and_samples(1,:,:));
[SignificantVariables,crit_p,~,adjusted_pvalues] = fdr_bh(pvalues,q_value,'pdep');

end  
                   