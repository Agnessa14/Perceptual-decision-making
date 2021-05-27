function [significantVarWei,significantVarMax,pValWei,pValMax] = rsa_weighted_cluster_perm_stats(subjects,rdm_cat_flattened,rdm_dis_flattened,true_correlation,tail,if_save,numPermutations)
%RSA_WEIGHTED_CLUSTER_PERM_STATS Perform weighted cluster permutation stats to calculate the
%significance of the timepoints for the object decoding RSA.
%
%Input:
% - subject IDs 
% - rdm_cat_flattened & rdm_dis_flattened: upper diagonal in vector form of the RDM NxNxP matrices of 1-Pearson's coefficient values,
%   where N is the number of conditions and P is the number of timepoints,
%   for each of the two tasks
% - true correlation: PxP vector of correlations between rdm_cat_flattened
%   and rdm_dis_flattened; ground truth
% - tail: 'left', only consider correlations smaller than the ground truth
%   as significant; 'both', consider both smaller & larger correlations
% - save (save the structure (1) or not (0))
% - number of permutations for the stats (ex:10000)
%
%Output: 
% - significanceVarWei: 1xP vector of significance (1) or not (0) based on
%   cluster weight (not sure if this part of the script is correct)
% - significanceVarMax: 1xP vector of significance (1) or not based on
%   cluster size
% - pValWei: p-value of the max cluster in the data based on weight
% - pValMax: p-value of the max cluster in the data based on size
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

                 %%%%% CALCULATING THE GROUND TRUTH AND PERMUTATION SAMPLES P-VALUES %%%%%

%% 1)Permute the categorization RDM X times and calculate the correlation with the distraction RDM at each timepoint
numTimepoints = 200;
sample_correlations = NaN(numPermutations,numTimepoints,numTimepoints);

for perm = 800:numPermutations
    if ~mod(perm,100)
        fprintf('Calculating the correlation %d \n',perm);
    end
    permuted_rdm_cat = rdm_cat_flattened(randperm(size(rdm_cat_flattened,1)),:);
    for tp1 = 1:numTimepoints
        for tp2 = 1:numTimepoints
            sample_correlations(perm,tp1,tp2) = 1-corr(permuted_rdm_cat(:,tp1),rdm_dis_flattened(:,tp2),'type','Spearman');
        end
    end
end

%% 2) Calculate the p-value of the ground truth and of the permuted samples
all_correlations = NaN(numPermutations+1,numTimepoints,numTimepoints);
all_correlations(1,:,:) = true_correlation;
all_correlations(2:end,:,:) = sample_correlations;

if strcmp(tail,'left')
    p_ground_and_samples = (numPermutations+1 - tiedrank(all_correlations*-1)) / numPermutations;
elseif strcmp(tail,'both')
    p_ground_and_samples = (numPermutations+1 - tiedrank(abs(all_correlations))) / numPermutations;
end 
                %%%%% CLUSTER-BASED PT: ATTRIBUTE SIGNIFICANCE TO TIMEPOINTS %%%%%

%% 1) Find maximum weighted cluster (highest sum of correlations) in the ground truth and the permutation samples
cluster_th = 0.05;
clusterMaxSize = NaN(numPermutations,1);
clusterMaxWei = NaN(numPermutations,1);

%ground truth
p_ground = squeeze(p_ground_and_samples(1,:,:));
if strcmp(tail,'left')
    [clusterMaxSize_ground, clusterMaxWei_ground, clusters, clustersize, clusterweight]  = find_clusters_weight_alld(true_correlation*-1, p_ground<=cluster_th);
elseif strcmp(tail,'both')
    [clusterMaxSize_ground, clusterMaxWei_ground, clusters, clustersize, clusterweight]  = find_clusters_weight_alld(abs(true_correlation), p_ground<=cluster_th);
end

%permutation samples
p_samples = p_ground_and_samples(2:end,:,:);
for perm = 1:numPermutations
    if ~mod(perm,100)
        fprintf('Permutation %d \n',perm);
    end
    if strcmp(tail,'left')
        [clusterMaxSize(perm), clusterMaxWei(perm)] = find_clusters_weight_alld(squeeze(sample_correlations(perm,:,:)*-1), squeeze(p_samples(perm,:)<=cluster_th));
    elseif strcmp(tail,'both')
        [clusterMaxSize(perm), clusterMaxWei(perm)] = find_clusters_weight_alld(squeeze(abs(sample_correlations(perm,:,:))), squeeze(p_samples(perm,:)<=cluster_th));
    end
end           

%combine ground and permutation samples
clusterMaxSize_all = [clusterMaxSize_ground;clusterMaxSize];
clusterMaxWei_all = [clusterMaxWei_ground;clusterMaxWei];

%% 2) Compare ground truth cluster to the permutation samples and identify significant timepoints, if any
significance_th = 0.05;

%find cluster threshold
clusterMaxSize_sorted = sort(clusterMaxSize_all, 'descend');
clusterMaxWei_sorted = sort(clusterMaxWei_all, 'descend');
th_max = clusterMaxSize_sorted( numPermutations*significance_th );
th_wei = clusterMaxWei_sorted( numPermutations*significance_th );

%preallocate significant variables
significantVarMax = zeros(numTimepoints,numTimepoints);
significantVarWei = zeros(numTimepoints,numTimepoints);

%apply threshold on found clusters
significantVarMax([clusters{clustersize>th_max}]) = 1;
significantVarWei([clusters{clusterweight>th_wei}]) = 1;
if ~isempty (clustersize)
    pValMax = ( find(clusterMaxSize_sorted == clusterMaxSize_ground, 1, 'first') ) / length(clusterMaxSize_sorted);
    pValWei = ( find(clusterMaxWei_sorted == clusterMaxWei_ground, 1, 'first') ) / length(clusterMaxWei_sorted);
    disp(['pValue(max) = ', num2str(pValMax), '; pValue(weighted) = ', num2str(pValWei), '.']);
else
    pValMax = NaN;
    pValWei = NaN;
    disp('No cluster found.');
end

%% Save        
if if_save == 1 && ~isempty(clustersize)
    permutation_stats.SignificantMaxClusterSize = significantVarMax;
    permutation_stats.SignificantMaxClusterWeight = significantVarWei;
    permutation_stats.pValueClusterSize = pValMax;
    permutation_stats.pValueWeight = pValWei;
    permutation_stats.info.num_permutations = numPermutations;
    permutation_stats.info.cluster_th = cluster_th;
    permutation_stats.info.significance_th = significance_th;
    permutation_stats.info.tail = tail;
    
    save(fullfile(results_avg,sprintf('permutation_stats_rsa_cross_task_subjects_%d_%d',...
        subjects(1),subjects(end)) ),'permutation_stats');
end

end  
                   