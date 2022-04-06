function [dataMatrixPT,numPseudotrials] = create_pseudotrials(numTrialsBin,dataMatrix) 
%CREATE_PSEUDOTRIALS Transform the data matrix by collapsing the data
%across trials into pseudotrials (bins containing a number of trials).
%
%Input: 
%   -numTrialsBin: number of trials per bin/pseudotrial (integer)
%   -dataMatrix: matrix of data containing voltage information for each condition, trial, electrode and timepoint,
%   as obtained by create_data_matrix.
%
%Output: 
%   -dataMatrixPT: NxPTxExT matrix of data reorganized into pseudotrials.
%
%
numConditions = size(dataMatrix,1);
minNumTrials  = size(dataMatrix,2);
numElectrodes = size(dataMatrix,3);
numTimepoints = size(dataMatrix,4);
numPseudotrials=round(minNumTrials/numTrialsBin); %must be an integer; number of pseudo trials
dataMatrixPT=NaN(numConditions,numPseudotrials,numElectrodes,numTimepoints); %pre-allocate memory
for step=1:numPseudotrials-1 %average by steps
    trial_selector=(1+(step-1)*numTrialsBin):(numTrialsBin+(step-1)*numTrialsBin); %select trials to be averaged
    dataMatrixPT(:,step,:,:)= nanmean(dataMatrix(:,trial_selector,:,:),2); %assign pseudo trial to dataMatrixPT
end
dataMatrixPT(:,numPseudotrials,:,:) = nanmean(dataMatrix(:,(1+(numPseudotrials-1)*numTrialsBin):end,:,:),2);