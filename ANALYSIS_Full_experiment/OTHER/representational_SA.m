function rdm_rsa = representational_SA(rdm_1,rdm_2,numTimepoints)
%REPRESENTATIONAL_SA Perform the representational similarity analysis on
%two representational dissimilarity matrices of choice, over time.
%
%Input: rdm_1, rdm_2 (NxNxP matrices containing 1-Pearson's coefficient
%values, where N is the number of conditions and P is the number of
%timepoints), numTimepoints is P 
%
%Output: rdm_rsa, a PxP matrix of 1-Spearman's coefficient values
%
%
%Author: Agnessa Karapetian, 2021
%

%%  Reshape the matrices: take only the upper diagonal, in vector form
rdm_1(isnan(rdm_1)) = 0;
rdm_2(isnan(rdm_2)) = 0;

rdm_1_flattened_cell = arrayfun(@(x) squareform(rdm_1(:,:,x)+(rdm_1(:,:,x))'),...
    1:numTimepoints,'UniformOutput',false);
rdm_2_flattened_cell = arrayfun(@(x) squareform(rdm_2(:,:,x)+(rdm_2(:,:,x))'),...
    1:numTimepoints,'UniformOutput',false);

rdm_1_flattened = reshape(cell2mat(rdm_1_flattened_cell),[],numTimepoints);
rdm_2_flattened = reshape(cell2mat(rdm_2_flattened_cell),[],numTimepoints);

%% Perfom RSA at each combination of timepoints
rdm_rsa = NaN(numTimepoints,numTimepoints);
for tp1 = 1:numTimepoints
    for tp2 = 1:numTimepoints
        rdm_rsa(tp1,tp2) = 1-corr(rdm_1_flattened(:,tp1),rdm_2_flattened(:,tp2),'type','Spearman');
    end
end

end