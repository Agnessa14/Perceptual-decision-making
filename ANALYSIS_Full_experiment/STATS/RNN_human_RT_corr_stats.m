function RNN_human_RT_corr_stats(numPermutations,alpha,model_type)
%RNN_HUMAN_RT_CORR_STATS Run the stats for the correlation between human
%and RNN reaction times (right-tailed). 
%
%Returns a correlation value and whether it's significant or not. 
%
%Input: numPermutations (e.g., 1000), alpha value for the FDR correction
%(e.g., 0.05), model_type ('b', 'bl' or b_d')
%
%Author: Agnessa Karapetian, 2021
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load subject-level RNN-human correlations
if strcmp(model_type,'bl')
    load(fullfile(results_avg_dir,'02.11_2_rnn/Model_RDM_redone','correlation_RT_human_RNN_cross-validated.mat'),'data');
elseif strcmp(model_type,'b_d')
    load(fullfile('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/DNN/','correlation_RT_human_b_d_net_7_layers_cross-validated.mat'),'data');
elseif strcmp(model_type,'b')
    load(fullfile('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/DNN/','correlation_RT_human_b_net_cross-validated.mat'),'data');
end

rnn_human_corr_all = data;
rng('shuffle');
pvalues = NaN(3,1);

for c = 1:3
    correlations_conds = rnn_human_corr_all(:,c); 
    samples_plus_ground_tstatistic = NaN(numPermutations,1);
    samples_plus_ground_tstatistic(1) = mean(correlations_conds)./ std(correlations_conds); 

    for perm = 2:numPermutations    
        if ~mod(perm,100)
            fprintf('Permutation sample %d \n',perm);
        end   
        random_vector = single(sign(rand(size(correlations_conds,1),1)-0.5));  %create samples by randomly multiplying each subject's data by 1 or -1
        sample = random_vector.*correlations_conds;
        samples_plus_ground_tstatistic(perm) = mean(sample) ./ std(sample);
    end

    %% 2) Calculate the p-value of the ground truth and of the permuted samples
    p_ground_and_samples = (numPermutations+1 - tiedrank(samples_plus_ground_tstatistic)) / numPermutations;
    pvalues(c) = p_ground_and_samples(1); 
 
end

%% 3) Add FDR correction and check significance
qvalue = 0.05;
[significance,crit_p,~,~] = fdr_bh(pvalues,qvalue,'pdep');

%% Save
stats_RT_corr.numPerm = numPermutations;
stats_RT_corr.alpha = alpha;
stats_RT_corr.pvalues = pvalues;
stats_RT_corr.significance = significance;
stats_RT_corr.correlation = mean(rnn_human_corr_all,1);
stats_RT_corr.crit_p = crit_p;
save(fullfile(results_avg_dir,sprintf('stats_%s_human_RT_correlation',model_type)),'stats_RT_corr');

end