function [dataMatrix,trialsMatAll] = create_data_matrix_subset(numConditions,triggers,minNumTrials,data)
%CREATE_DATA_MATRIX_SUBSET Construct the matrix of data containing voltage
%information for each condition, trial, electrode and timepoint, only using a specified subset of data. 
%
%Input:  
%   -numConditions: number of conditions (integer)
%   -triggers: list of triggers as they appear in the experiment (Tx1,
%   where T is the number of overall trials)
%   -minNumTrials: the overall minimum number of trials as obtained by the
%   min_number_trials function
%   -data: TxExTP matrix containing the EEG data
%   -subset: value indicating the proportion of data to keep (0 to 1)
%
%Output: NxMxExTP matrix containing the data in the desired format.
%
numElectrodes = size(data,2);
numTimepoints = size(data,3);
dataMatrix = NaN(numConditions,minNumTrials,numElectrodes,numTimepoints);
trialsMatAll  = NaN(numConditions,minNumTrials);
for c = 1:numConditions
    trialsMat = find(triggers==c);
    trialsMat = trialsMat(1:minNumTrials);
    trialsMat = trialsMat(randperm(numel(trialsMat))); %randomize
    trialsMatAll(c,:) = trialsMat;
    dataMatrix(c,:,:,:) = data(trialsMat,:,:); 
end