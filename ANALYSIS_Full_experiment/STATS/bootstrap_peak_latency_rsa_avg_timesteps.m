function bootstrap_peak_latency_rsa_avg_timesteps(subjects,conditions_int)
%BOOTSTRAP_PEAK_LATENCY_RSA_AVG_TIMESTEPS Apply bootstrapping to calculate the 95% confidence
%interval of peak latency difference (for RSA between RNN and
%EEG-categorization task). Results are for all RNN timepoints & for the median. 
%
%Input: subject IDs (e.g., 1:13), conditions (1(all),2(artificial) or
%3(natural)
%
%Output: saved boostrapped peak latencies (in ms), 95% confidence interval of the difference
%
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
rnn_dir = '02.11_2_rnn/Input_RDM';

                                    %%%%% SETUP %%%%%
%% Load the subject-level RDMs
if conditions_int==1
    conds = 1:60;
    conditions = 'all';
elseif conditions_int==2
    conds = 1:30;
    conditions = 'artificial';
elseif conditions_int==3
    conds = 31:60;
    conditions='natural';
end

numConditions = numel(conds);
numTimepoints = 200;
numTimepointsRNN = 8;
layers_idx = [1,4,7]; %RNN layers of interest
rdm_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints); 

for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'rdm_pearson_categorization.mat'),'rdm_avg');           
    rdm_all_subjects(subject,:,:,:) = rdm_avg(conds,conds,:);      
end   

rdm_subjects = rdm_all_subjects(subjects,:,:,:); %remove NaN subjects

%% Load the RNN RDMs & RSA results
%results for 3 layers and all timepoints
filename = sprintf('02_11_2_redone_rsa_plus_noise_ceilings_subject_level_-100ms_%s_subjects_5_35.mat',conditions);
load(fullfile(results_avg_dir,filename),'rsa_results');

%Preallocate
rdm_rnn = NaN(numel(layers_idx),numTimepointsRNN,numConditions,numConditions);
peak_latency_ground_truth = NaN(numel(layers_idx),numTimepointsRNN);

for layer = layers_idx
    index_layer = find(layers_idx==layer);
    %RDMs
    for t = 1:numTimepointsRNN
        load(fullfile(results_avg_dir,rnn_dir,sprintf('ReLU_Layer_%d_Time_%d_Input_RDM.mat',layer-1,t-1)),'data');
        rdm_rnn(index_layer,t,:,:) = data(conds,conds);
        for c1 = 1:numConditions
            for c2 = 1:numConditions
                if c1==c2
                    rdm_rnn(index_layer,t,c1,c2) = 0; %Make sure the diagonal is all 0
                end
            end
        end
        %RSA results: get ground truth peak latency (for each timepoint)
        rsa_results_temp = squeeze(rsa_results(index_layer,t,:))';
        peak_latency_ground_truth(index_layer,t) = (find(squeeze(rsa_results_temp)==max(rsa_results_temp))-40)*5;    
    end
end

                                  %%%%% BOOTSTRAPPING %%%%%
                                  
%% 0) Prepare required variables
num_bootstrap_samples = 1000;
num_datasets = size(rdm_subjects,1);
peak_latency = NaN(numel(layers_idx),numTimepointsRNN,num_bootstrap_samples);
peak_latency_median_tps = NaN(numel(layers_idx),num_bootstrap_samples);
max_dataset = num_datasets; %for the random dataset index generator
min_dataset = 1;

%% 1) Create the bootstrap samples & calculate the peak difference
rng(0,'twister');

%for each bootstrap sample, create a new dataset (N=30) and
%calculate the peak latency
for bs = 1:num_bootstrap_samples
    rsa = NaN(num_datasets,numel(layers_idx),numTimepointsRNN,numTimepoints);   
    
    %create sample and perform RSA on each dataset
    for d = 1:num_datasets
        idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
        dataset = squeeze(rdm_subjects(idx,:,:,:));        
        for layer = layers_idx
            index_layer = find(layers_idx==layer);            
            for t = 1:numTimepointsRNN                       
                rsa(d,index_layer,t,:) = representational_SA_rnn(dataset,squeeze(rdm_rnn(index_layer,t,:,:)));
            end
        end
    end
    avg_rsa_datasets = squeeze(mean(rsa,1)); %avg over datasets
    avg_rsa_tps = squeeze(median(avg_rsa_datasets,2)); %avg over timepoints
    
    %calculate sample peak latency
    for layer = layers_idx
        index_layer = find(layers_idx==layer);            
        peak_latency_median_tps(index_layer,bs) = (find(avg_rsa_tps(index_layer,:)==max(avg_rsa_tps(index_layer,:)),1)-40)*5;
        for t = 1:numTimepointsRNN
            peak_latency(index_layer,t,bs) = (find(avg_rsa_datasets(index_layer,t,:)==max(avg_rsa_datasets(index_layer,t,:)),1)-40)*5;
        end
    end
end

%% 2) Calculate difference between intervals for bootstrap samples 
%for each RNN timepoint 
layer1_layer4_all_tps = abs(squeeze(peak_latency(1,:,:))-squeeze(peak_latency(2,:,:)));
layer4_layer7_all_tps = abs(squeeze(peak_latency(2,:,:))-squeeze(peak_latency(3,:,:)));
layer1_layer7_all_tps = abs(squeeze(peak_latency(1,:,:))-squeeze(peak_latency(3,:,:)));

%for the median
layer1_layer4_median_tps = abs(squeeze(peak_latency_median_tps(1,:))-squeeze(peak_latency_median_tps(2,:)));
layer4_layer7_median_tps = abs(squeeze(peak_latency_median_tps(2,:))-squeeze(peak_latency_median_tps(3,:)));
layer1_layer7_median_tps = abs(squeeze(peak_latency_median_tps(1,:))-squeeze(peak_latency_median_tps(3,:)));

%% 3) Get 95% confidence interval for peaks and for peak difference
%all timepoints
%peaks 
CI_peaks_all_tps = NaN(numel(layers_idx),numTimepointsRNN,2);

for layer = layers_idx
    index_layer = find(layers_idx==layer);
    for t = 1:numTimepointsRNN
        CI_peaks_all_tps(index_layer,t,1) = prctile(squeeze(peak_latency(index_layer,t,:)),2.5);
        CI_peaks_all_tps(index_layer,t,2) = prctile(squeeze(peak_latency(index_layer,t,:)),97.5);
    end
end

%peak differences
CI_diff_l1_l4_all_tps = NaN(numTimepointsRNN,2);
CI_diff_l4_l7_all_tps = NaN(numTimepointsRNN,2);
CI_diff_l1_l7_all_tps = NaN(numTimepointsRNN,2);

for t = 1:numTimepointsRNN
    CI_diff_l1_l4_all_tps(t,1) = prctile(squeeze(layer1_layer4_all_tps(t,:)),2.5);
    CI_diff_l1_l4_all_tps(t,2) = prctile(squeeze(layer1_layer4_all_tps(t,:)),97.5);
    CI_diff_l4_l7_all_tps(t,1) = prctile(squeeze(layer4_layer7_all_tps(t,:)),2.5);
    CI_diff_l4_l7_all_tps(t,2) = prctile(squeeze(layer4_layer7_all_tps(t,:)),97.5);
    CI_diff_l1_l7_all_tps(t,1) = prctile(squeeze(layer1_layer7_all_tps(t,:)),2.5);
    CI_diff_l1_l7_all_tps(t,2) = prctile(squeeze(layer1_layer7_all_tps(t,:)),97.5);        
end

%median of timepoints
%peaks
CI_peaks_median_tps = NaN(numel(layers_idx),2);

for layer = layers_idx
    index_layer = find(layers_idx==layer);
    CI_peaks_median_tps(index_layer,1) = prctile(squeeze(peak_latency_median_tps(index_layer,:)),2.5);
    CI_peaks_median_tps(index_layer,2) = prctile(squeeze(peak_latency_median_tps(index_layer,:)),97.5);
end

%peak differences
CI_diff_l1_l4_median_tps = NaN(1,2);
CI_diff_l4_l7_median_tps = NaN(1,2);
CI_diff_l1_l7_median_tps = NaN(1,2);

CI_diff_l1_l4_median_tps(1) = prctile(layer1_layer4_median_tps,2.5);
CI_diff_l1_l4_median_tps(2) = prctile(layer1_layer4_median_tps,97.5);
CI_diff_l4_l7_median_tps(1) = prctile(layer4_layer7_median_tps,2.5);
CI_diff_l4_l7_median_tps(2) = prctile(layer4_layer7_median_tps,97.5);
CI_diff_l1_l7_median_tps(1) = prctile(layer1_layer7_median_tps,2.5);
CI_diff_l1_l7_median_tps(2) = prctile(layer1_layer7_median_tps,97.5);        


%% Save as structure
bootstrap_peak_latencies_rsa.peak_latency_bs_all_tps = squeeze(mean(peak_latency,3));
bootstrap_peak_latencies_rsa.peak_latency_bs_median_tps = squeeze(mean(peak_latency_median_tps,2));
bootstrap_peak_latencies_rsa.peak_latency_true_all_tps = peak_latency_ground_truth;
bootstrap_peak_latencies_rsa.peak_latency_true_median_tps = median(peak_latency_ground_truth,2);

bootstrap_peak_latencies_rsa.CI_peaks_all_tps = CI_peaks_all_tps;
bootstrap_peak_latencies_rsa.CI_diff_l1_l4_all_tps = CI_diff_l1_l4_all_tps;
bootstrap_peak_latencies_rsa.CI_diff_l4_l7_all_tps = CI_diff_l4_l7_all_tps;
bootstrap_peak_latencies_rsa.CI_diff_l1_l7_all_ps = CI_diff_l1_l7_all_tps;

bootstrap_peak_latencies_rsa.CI_peaks_median_tps = CI_peaks_median_tps;
bootstrap_peak_latencies_rsa.CI_diff_l1_l4_median_tps = CI_diff_l1_l4_median_tps;
bootstrap_peak_latencies_rsa.CI_diff_l4_l7_median_tps = CI_diff_l4_l7_median_tps;
bootstrap_peak_latencies_rsa.CI_diff_l1_l7_median_ps = CI_diff_l1_l7_median_tps;

save(fullfile(results_avg_dir,sprintf('bootstrap_peak_latencies_rsa_all_and_median_tps_%s_subjects_%d_%d',conditions,subjects(1),subjects(end))),'bootstrap_peak_latencies_rsa');

end
