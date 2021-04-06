function [dataMatrix] = create_data_matrix(numConditions,triggers,minNumTrials,data)
%CREATE_DATA_MATRIX Construct the matrix of data containing voltage
%information for each condition, trial, electrode and timepoint. 
%
%Input:  
%   -numConditions: number of conditions (integer)
%   -triggers: list of triggers as they appear in the experiment (Tx1,
%   where T is the number of overall trials)
%   -minNumTrials: the overall minimum number of trials as obtained by the
%   min_number_trials function
%   -data: TxExTP matrix containing the EEG data
%
%Output: NxMxExTP matrix containing the data in the desired format.
%

rng('shuffle');
numElectrodes = size(data,2);
numTimepoints = size(data,3);
dataMatrix = NaN(numConditions,minNumTrials,numElectrodes,numTimepoints);
for c = 1:numConditions
    trialsMat = find(triggers==c);
    trialsMat_randomized = trialsMat(randperm(numel(trialsMat))); %randomize
    if numel(trialsMat_randomized)<minNumTrials
        continue;
    else
        trialsMat_minimum = trialsMat_randomized(1:minNumTrials);
        dataMatrix(c,:,:,:) = data(trialsMat_minimum,:,:); 
    end
end