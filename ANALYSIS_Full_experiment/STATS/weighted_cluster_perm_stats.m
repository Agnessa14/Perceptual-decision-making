function [significantVarWei,significantVarMax,pValWei,pValMax] = weighted_cluster_perm_stats(medianRT,mean_distances,true_correlation,task_distance,task_RT,category,tail,if_save,numPermutations,method)
%WEIGHTED_CLUSTER_PERM_STATS Perform weighted cluster permutation stats to calculate the
%significance of the timepoints in the distance-to-hyperplane analysis.
%
%Input: 
% - medianRT: matrix of reaction times, either averaged over subjects (fixed effects), or per subject (random effects)
% - task for the EEG data (1 = categorization, 2 = distraction),
% - task for the reaction times (1 = categorization, 2 = distraction)
% - category ('artificial', 'natural', 'average' or 'both')
% - save (save the structure (1) or not (0))
% - number of permutations for the stats (ex:10000)
% - method of analysis ('fixed', averaging over all RTs and
% distances & correlating, or 'random', averaging over all RTs &
% correlating with each subject's distances)
%
%Output: 1xP vector of significance (1) or not (0), where P is the number
%of timepoints.
%
%Author: Agnessa Karapetian, 2021

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

                 %%%%% CALCULATING THE GROUND TRUTH AND PERMUTATION SAMPLES P-VALUES %%%%%

%% 1)Fixed effects: Permute the objects' RTs 10 000 times and calculate the correlation at each timepoint
task_distance_name = get_task_name(task_distance);
numTimepoints = 200;
sample_correlations = NaN(numPermutations,numTimepoints);
if strcmp(analysis,'random')
   subject_corr = NaN(size(medianRT,1),numTimepoints);
end
for perm = 1:numPermutations
    if ~mod(perm,100)
        fprintf('Calculating the correlation %d \n',perm);
    end
    if strcmp(analysis,'fixed')
        permuted_RTs = medianRT(randperm(numel(medianRT)));
        sample_correlations(perm,:) = arrayfun(@(x) corr(mean_distances(:,x),permuted_RTs','type','Spearman'),t);
    elseif strcmp(analysis,'random')
        permuted_RTs = medianRT(:,randperm(size(medianRT,2))); %permute for each subject
        for subject = subjects
            subject_corr(subject,:) = arrayfun(@(x) corr(squeeze(mean_distances(subject,:,x))',permuted_RTs','type','Spearman'), t);
        end
        sample_correlations(perm,:) = squeeze(mean(subject_corr,1));
    end     
end

%% 2) Calculate the p-value of the ground truth and of the permuted samples
all_correlations = [true_correlation;sample_correlations];
if strcmp(tail,'left')
    p_ground_and_samples = (numPermutations+1 - tiedrank(all_correlations*-1)) / numPermutations;
elseif strcmp(tail,'both')
    p_ground_and_samples = (numPermutations+1 - tiedrank(abs(all_correlations))) / numPermutations;
end 
                %%%%% CLUSTER-BASED PT: ATTRIBUTE SIGNIFICANCE TO TIMEPOINTS %%%%%

%% 1) Find maximum weighted cluster (highest sum of correlations) in the ground truth and the permutation samples
%find maximum cluster size and maximum weighted cluster for all permutation samples
cluster_th = 0.05;
clusterMaxSize = NaN(numPermutations,1);
clusterMaxWei = NaN(numPermutations,1);

%ground truth
p_ground = squeeze(p_ground_and_samples(1,:));
if strcmp(tail,'left')
    [clusterMaxSize_ground, clusterMaxWei_ground, clusters, clustersize, clusterweight]  = find_clusters_weight_alld(true_correlation*-1, p_ground<=cluster_th);
elseif strcmp(tail,'both')
    [clusterMaxSize_ground, clusterMaxWei_ground, clusters, clustersize, clusterweight]  = find_clusters_weight_alld(abs(true_correlation), p_ground<=cluster_th);
end

%permutation samples
p_samples = p_ground_and_samples(2:end,:);
for perm = 1:numPermutations
    if ~mod(perm,100)
        fprintf('Permutation %d \n',perm);
    end
    if strcmp(tail,'left')
        [clusterMaxSize(perm), clusterMaxWei(perm)] = find_clusters_weight_alld(squeeze(sample_correlations(perm,:)*-1), squeeze(p_samples(perm,:)<=cluster_th));
    elseif strcmp(tail,'both')
        [clusterMaxSize(perm), clusterMaxWei(perm)] = find_clusters_weight_alld(squeeze(abs(sample_correlations(perm,:))), squeeze(p_samples(perm,:)<=cluster_th));
    end
end           

%combine ground and permutation samples
clusterMaxSize_all = [clusterMaxSize_ground;clusterMaxSize];
clusterMaxWei_all = [clusterMaxWei_ground;clusterMaxWei];

%% 3) Compare ground truth cluster to the permutation samples and identify significant timepoints, if any
significance_th = 0.05;

%find cluster threshold
clusterMaxSize_sorted = sort(clusterMaxSize_all, 'descend');
clusterMaxWei_sorted = sort(clusterMaxWei_all, 'descend');
th_max = clusterMaxSize_sorted( numPermutations*significance_th );
th_wei = clusterMaxWei_sorted( numPermutations*significance_th );

%preallocate significant variables
significantVarMax = zeros(numTimepoints,1);
significantVarWei = zeros(numTimepoints,1);

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
    
    if isequal(task_distance,task_RT)
        filename = 'dth_permutation_stats';
    else
        filename = 'dth_permutation_stats_crosstask';
    end
    
    if strcmp(method,'random')
        filename = sprintf('random_%s',filename);
    end
    
    save(fullfile(results_avg,sprintf('%s_%s_%s_distance_subjects_%d_%d',...
        filename,category,task_distance_name,subjects(1),subjects(end)) ),'permutation_stats');
end

end  
                   