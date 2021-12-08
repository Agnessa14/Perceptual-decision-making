 function [SignificantVariables, pvalues, crit_p, adjusted_pvalues] = fdr_permutation_cluster_1sample_alld(data,nperm,tail,q_value,StatMapPermPV)

% Performs one-sided (>0) or two-sided pemutation test + FDR correction on the 'data'.
% The data array can have any number of dimensions (1D, 2D, 3D, 4D etc) in
% the order [observations x variable1 x variable2 x ...]
% Data from each observation are randomly multiplied by +-1 to create
% permutation samples. These samples are converted to pvalues and then the
% fdr_bh function is applied to identify statistically significant points.
%

% INPUT:
%   data: observations x variable1 x variable2 x variable3 x ... (supports all dimensions)
%   nperm: number of permutations
%   q_value: alpha level for FDR correction
%   tail: string 'right' or 'both'. Default is 'both'.
%   StatMapPermPV: (optional) permutation pvalue map (see output)
%
% OUTPUT:
%   SignificantVariables:  times with significant values
%   pvalues: p-values of the data (ground truth)
%   crit_p: critical (statictically decisive) p-value in FDR
%   adjusted_pvalues: FDR corrected p-values of the ground truth
%
% Author: Dimitrios Pantazis, December 2015
% Modified: Agnessa Karapetian, June 2021


%decide one-sided (right) or two-sided (both) test
if ~exist('tail') || strcmp(tail,'right') %if two tail t-test
    func = '';
else %if two-sided t-test
    func = 'abs';
end

%initialize
N = ndims(data); %data is observations x variable1 x variable2 x ... (supports all dimensions)
nobservations = size(data,1);

for n = 2:N
    nvariable(n-1) = size(data,n);
end

cln = repmat({':'},1,N-1); %to select the N-1 dimensions

%create permutation samples and convert them to pvalues (perms x variable1 x variable2 x ...)
if ~exist('StatMapPermPV') %if pvalues have not been precomputed     
    StatMapPerm = single(zeros([nperm nvariable]));
    
    %first permutation sample is original data
    StatMapPerm(1,cln{:}) = mean(data,1) ./ std(data);    

    %perform permutations
    for i = 2:nperm %par
        if ~rem(i,100)
            disp(['Create permutation samples: ' num2str(i) ' out of ' num2str(nperm)]);
        end
        perm = single(sign(rand(nobservations,1)-0.5));
        data_perm = repmat(perm,1,nvariable(1)) .* data; %        
        StatMapPerm(i,cln{:}) = mean(data_perm,1) ./ std(data_perm);         
    end    

    %convert to pvalues
    eval([ 'StatMapPermPV = (nperm+1 - tiedrank(' func '(StatMapPerm)))/nperm;' ]);    
end

clear StatMapPerm;

%Perform fdr
pvalues = squeeze(StatMapPermPV(1,cln{:}));
[SignificantVariables,crit_p,~,adjusted_pvalues] = fdr_bh(pvalues,q_value,'pdep');

 end