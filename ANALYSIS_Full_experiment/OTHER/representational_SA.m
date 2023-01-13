function [rdm_rsa,rdm_1_flattened,rdm_2_flattened] = representational_SA(rdm_1,rdm_2,numTimepoints)
%REPRESENTATIONAL_SA Perform the representational similarity analysis on
%two representational dissimilarity matrices of choice, over time.
%
%Input: rdm_1, rdm_2 (either as NxNxP matrices containing 1-Pearson's coefficient
%values, where N is the number of conditions and P is the number of
%timepoints, or as flattened versions of said matrices), numTimepoints is P
%
%Output: rdm_rsa, a PxP matrix of 1-Spearman's coefficient values
%
%
%Author: Agnessa Karapetian, 2021
%

%%  Reshape the matrices: take only the upper diagonal, in vector form
for r = 1:2
    if r == 1
        rdm = rdm_1;
    elseif r == 2
        rdm = rdm_2;
    end
    if find(isnan(rdm)) >0 
        rdm(isnan(rdm)) = 0;
        rdm_flattened_cell = arrayfun(@(x) squareform(rdm(:,:,x)+(rdm(:,:,x))'),...
            1:numTimepoints,'UniformOutput',false);
        rdm_flattened = reshape(cell2mat(rdm_flattened_cell),[],numTimepoints);
    elseif numel(size(rdm)) == 2 && (size(rdm,2) == numTimepoints) %if it's a flattened rdm 
        rdm_flattened = rdm;
    end
    if r == 1 
        rdm_1_flattened = rdm_flattened;
    elseif r == 2
        rdm_2_flattened = rdm_flattened;
    end
end


%% Perfom RSA at each combination of timepoints
rdm_rsa = NaN(numTimepoints,numTimepoints);
for tp1 = 1:numTimepoints
    for tp2 = 1:numTimepoints
        rdm_rsa(tp1,tp2) = corr(rdm_1_flattened(:,tp1),rdm_2_flattened(:,tp2),'type','Spearman');
    end
end

end