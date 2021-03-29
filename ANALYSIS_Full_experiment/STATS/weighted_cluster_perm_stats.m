function significant_timepoints = weighted_cluster_perm_stats(subjects,task,category)
%WEIGHTED_CLUSTER_PERM_STATS Perform weighted cluster permutation stats to calculate the
%significance of the timepoints in the distance-to-hyperplane analysis.
%
%Input: subject IDs, with_stats (1 plot with stats, 0 plot without),
%category ('artificial', 'natural', 'average' or 'both')
%
%Output: 1xP vector of significance (1) or not (0), where P is the number
%of timepoints.
%
%
%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';

%% Get the distances from all subjects
task_name = get_task_name(task);
numTimepoints = 200;
numConditions = 60;
sorted_subjects = sort(subjects); %order by ID
distances = NaN(sorted_subjects(end),numConditions,numTimepoints);
RTs = NaN(numel(subjects),numConditions);

for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,sprintf('dth_pseudotrials_svm_decisionValues_%s.mat',task_name)));
    load(fullfile(results_dir,subname,sprintf('RTs_correct_trials_%s.mat',task_name)));
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

                 %%%%% PERMUTATION TEST: CALCULATING THE GROUND TRUTH P-VALUES %%%%%
%% 1) Permute the objects' RTs 10 000 times and calculate the correlation at each timepoint
numPermutations = 10000;
sample_correlations = NaN(numPermutations,numTimepoints);
for perm = 1:numPermutations
    if ~mod(perm,100)
        fprintf('Calculating the correlation %d \n',perm);
    end
    permuted_RTs = medianRT(randperm(numel(medianRT)));
    sample_correlations(perm,:) = arrayfun(@(x) corr(mean_distances(:,x),permuted_RTs','type','Spearman'),t);
end

%% 2) Calculate the p-value of the true value WRT the permuted samples
p_ground_truth = NaN(numTimepoints,1);
 
for tp = 1:numTimepoints
    %calculate the b value: num of permutations larger than the ground truth 
    b = numel(find( abs(sample_correlations(:,tp)) > abs(true_correlation(tp)) )); %https://benediktehinger.de/blog/science/permutation-test-for-matlab/
    p_ground_truth(tp) = (b+1) / (numPermutations+1); 
end
%p_ground_truth = (numPermutations+1 - tiedrank(abs([true_correlation;all_correlations]))) / numPermutations;

           
                %%%%% CLUSTER-BASED PT: ATTRIBUTE SIGNIFICANCE TO TIMEPOINTS %%%%%
     
%% 1) Calculate the p-values of each permutation sample (at each timepoint)
p_samples = NaN(numPermutations,numTimepoints);
% all_correlations = [true_correlation;sample_correlations];
for perm = 1:numPermutations
%     p_samples(perm,tp) = (numPermutations+1 - tiedrank(abs(sample_correlations(perm,:)))) / numPermutations;
    if ~mod(perm,100)
        fprintf('Calculating the p-value %d \n',perm);
    end
    for tp = 1:numTimepoints
        b = numel(find( abs(sample_correlations(:,tp)) > abs(sample_correlations(perm,tp)) )); 
        p_samples(perm,tp) = (b+1) / (numPermutations+1); 
    end
end
% p_samples(perm,tp) = arrayfun(@(x,y) ...
%     (numPermutations+1 - tiedrank(abs([true_correlation(x);all_correlations(y,x)]))) / numPermutations, [1:numTimepoints,1:numPermutations]);

%% 2) Find maximum weighted cluster (highest sum of correlations) in the ground truth and the permutation samples
%find maximum cluster size and maximum weighted cluster for all permutation samples
cluster_th = 0.05;
clusterMaxSize = NaN(numPermutations,1);
clusterMaxWei = NaN(numPermutations,1);

%ground truth
[clusterMaxSize_ground, clusterMaxWei_ground, clusters, clustersize, clusterweight]  = find_clusters_weight_alld(p_ground_truth, p_ground_truth<=cluster_th);

%permutation samples
for perm = 1:numPermutations
    if ~mod(perm,100)
        fprintf('Permutation %d \n',perm);
    end
    [clusterMaxSize(perm), clusterMaxWei(perm)] = find_clusters_weight_alld(squeeze(p_samples(perm,:)), squeeze(p_samples(perm,:)<=cluster_th));
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
if ~isempty(clustersize)
    
end

end  
                   