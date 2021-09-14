function [rdm_rsa,rdm_flattened_eeg,rdm_flattened_rnn] = representational_SA_rnn(rdm_eeg,rdm_rnn)
%REPRESENTATIONAL_SA_RNN Perform the representational similarity analysis on
%two representational dissimilarity matrices of choice, over time. FOr the
%analysis with RNN.
%
%Input: rdm_eeg, rdm_rnn (EEG: either original NxNxP matrix, or its
%flattened upper diagonal, containing 1-Pearson's coefficient values, where N is the number of conditions and P is the number of
%EEG timepoints; RNN: NxN matrix at a given layer)
%
%Output: rdm_rsa, a PxP matrix of 1-Spearman's coefficient values
%
%
%Author: Agnessa Karapetian, 2021
%

%%  Reshape the matrices: take only the upper diagonal, in vector form
%EEG 
if find(isnan(rdm_eeg)) >0 %full matrix version
    numTimepoints_eeg = size(rdm_eeg,3);
    rdm_eeg(isnan(rdm_eeg)) = 0;
    rdm_flattened_cell_eeg = arrayfun(@(x) squareform(rdm_eeg(:,:,x)+(rdm_eeg(:,:,x))'),...
                1:numTimepoints,'UniformOutput',false);
    rdm_flattened_eeg = reshape(cell2mat(rdm_flattened_cell_eeg),[],numTimepoints_eeg);
else
    numTimepoints_eeg = size(rdm_eeg,2);
    rdm_flattened_eeg = rdm_eeg;
end

%RNN
rdm_flattened_rnn = squareform(rdm_rnn);

%% Perfom RSA at each EEG timepoint
rdm_rsa = NaN(1,numTimepoints_eeg);
for tp = 1:numTimepoints_eeg
    rdm_rsa(tp) = corr(rdm_flattened_eeg(:,tp),rdm_flattened_rnn','type','Spearman');
end

end