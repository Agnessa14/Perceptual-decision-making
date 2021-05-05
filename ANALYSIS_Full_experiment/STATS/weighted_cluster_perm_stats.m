function [significantVarWei,significantVarMax,pValWei,pValMax] = weighted_cluster_perm_stats(subjects,task_distance,task_RT,category,if_save,numPermutations,method)
%WEIGHTED_CLUSTER_PERM_STATS Perform weighted cluster permutation stats to calculate the
%significance of the timepoints in the distance-to-hyperplane analysis.
%
%Input: subject IDs, task (1 = categorization, 2 = distraction),
%category ('artificial', 'natural', 'average' or 'both'), save (save the
%structure (1) or not (0)), number of permutations for the stats (ex:
%10000), method of calculating the DTH ('fixed', averaging over all RTs and
%distances & correlating, or 'non-fixed', averaging over all RTs &
%correlating with each subject's distances)
%
%Output: 1xP vector of significance (1) or not (0), where P is the number
%of timepoints.
%
%
%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Get the distances from all subjects
task_distance_name = get_task_name(task_distance);
task_RT_name = get_task_name(task_RT);
numTimepoints = 200;
numConditions = 60;
sorted_subjects = sort(subjects); %order by ID
distances = NaN(sorted_subjects(end),numConditions,numTimepoints);
RTs = NaN(numel(subjects),numConditions);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_distance_name)));
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_RT_name)));
    distances(subject,:,:) = decisionValues_Avg;   
    RTs(subject,:) = normalize(RT_per_condition);
end

%% Get the median RTs and mean distances of all subjects for each condition 
medianRT = nanmedian(RTs,1);
mean_distances = squeeze(nanmean(distances,1)); %avg over subjects

%% Only keep the RTs and distances of the category for the analysis
if strcmp(category,'artificial')
    conditions = 1:30;
elseif strcmp(category,'natural')
    conditions = 31:60;
elseif strcmp(category,'average') || strcmp(category,'both')
    conditions = 1:60;
end

medianRT = medianRT(conditions);
mean_distances = mean_distances(conditions,:);

%% Calculate the true correlation value
t = 1:numTimepoints;
true_correlation = arrayfun(@(x) corr(mean_distances(:,x),medianRT','type','Spearman'),t);

                 %%%%% CALCULATING THE GROUND TRUTH AND PERMUTATION SAMPLES P-VALUES %%%%%

%% 1) Permute the objects' RTs 10 000 times and calculate the correlation at each timepoint
sample_correlations = NaN(numPermutations,numTimepoints);
for perm = 1:numPermutations
    if ~mod(perm,100)
        fprintf('Calculating the correlation %d \n',perm);
    end
    permuted_RTs = medianRT(randperm(numel(medianRT)));
    sample_correlations(perm,:) = arrayfun(@(x) corr(mean_distances(:,x),permuted_RTs','type','Spearman'),t);
end

%% 2) Calculate the p-value of the ground truth and of the permuted samples
tail = 'left'; %because we are interested in the strength of negative correlations
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
    
    if strcmp(method,'non-fixed')
        filename = sprintf('non_fixed_%s',filename);
    end
    
    save(fullfile(results_avg,sprintf('%s_%s_%s_distance_subjects_%d_%d',...
        filename,category,task_distance_name,subjects(1),subjects(end)) ),'permutation_stats');
end

end  
                   