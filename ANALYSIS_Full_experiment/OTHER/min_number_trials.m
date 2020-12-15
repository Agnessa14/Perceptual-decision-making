function [N, numTrials] = min_number_trials(conditionsList,numConditions)
%MIN_NUMBER_TRIALS Find the minimum number of trials for decoding or DTH.
%
%Input: 
%   -conditionsList (Tx1 list of conditions as they appear in the experiment, one condition per trial (T is the number of trials))
%   -numConditions (integer number of conditions)
%
%Output: 
%   -N: overall minimum number of trials (integer)
%   -numTrials: Nx1 vector including the minimum num of trials for each
%   condition (N is the number of conditions)
%

numTrials = [];
for cond = 1:numConditions
    numTrials(cond) = numel(find(conditionsList==cond));
end
N = min(numTrials(numTrials>0));