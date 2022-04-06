function x = multivariate_noise_normalization(dataMatrix) 
%MULTIVARIATE_NOISE_NORMALIZATION Perform MVNN on the data.
%
%Input: 
%   -finalMatrix: NxMxExTP matrix containing EEG data, where N is the
%   number of conditioins, M is the number of trials, E is the number of
%   electrodes and TP is the number of timepoints.
%
%Output: 
%   -x: matrix of the same size as the dataMatrix, but with normalized
%   values.

%data: condition x sensor x time  
%1. calculate sigma: covariance matrix for each timepoint and each condition, then avg across time & conditions
numConditions = size(dataMatrix,1);
numTrials     = size(dataMatrix,2);
numElectrodes = size(dataMatrix,3);
numTimepoints = size(dataMatrix,4);
sigma_covariance = zeros(numTimepoints,numConditions,numElectrodes,numElectrodes); 
for t = 1:numTimepoints
    for c = 1:numConditions
        X = (squeeze(dataMatrix(c,:,:,t))); 
        [sigma,~] = covCor(X);
        sigma_covariance(t,c,:,:) = sigma;
    end
end
sigma_covariance = squeeze(mean(mean(sigma_covariance,1),2)); 

%2. sigma^-1/2
inverted_sigma_covariance = sigma_covariance^(-1/2);

%3. mmn_x = inverted sigma^-1/2 * x  (whiten data)
x = NaN(size(dataMatrix));
for t = 1:numTimepoints %and for each condition
    for c = 1:numConditions
        for tr = 1:numTrials
            X = squeeze(dataMatrix(c,tr,:,t))'; 
            x(c,tr,:,t) = X*inverted_sigma_covariance;
        end
    end
end
