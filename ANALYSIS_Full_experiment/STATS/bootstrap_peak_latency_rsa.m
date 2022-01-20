function bootstrap_peak_latency_rsa(subjects,conditions_int)
%BOOTSTRAP_PEAK_LATENCY_RSA Apply bootstrapping to calculate the 95% confidence
%interval of peak latency difference (for RSA between RNN and EEG-categorization task).
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
                                  %%%%% BOOTSTRAPPING %%%%%

                                  
%% 0) Prepare required variables
num_bootstrap_samples = 1000;
num_datasets = size(rdm_subjects,1);
peak_latency = NaN(numel(layers_idx),numTimepointsRNN,num_bootstrap_samples);
max_dataset = num_datasets; %for the random dataset index generator
min_dataset = 1;

%% 1) Create the bootstrap samples & calculate the peak difference
rng(0,'twister');
peak_latency_ground_truth = NaN(numel(layers_idx),numTimepointsRNN);
%load the RSA results
filename = sprintf('02_11_2_redone_rsa_plus_noise_ceilings_subject_level_-100ms_%s_subjects_5_35.mat',conditions);

for layer = layers_idx
    index_layer = find(layers_idx==layer);
    load(fullfile(results_avg_dir,filename),'rsa_results');

    %bootstrap samples collected for each timepoint (RNN) separately
    for t = 1:numTimepointsRNN         
        rsa_results_tl = squeeze(rsa_results(index_layer,t,:))';
        
        %calculate ground truth peak latency 
        peak_latency_ground_truth(index_layer,t) = (find(squeeze(rsa_results_tl)==max(rsa_results_tl))-40)*5;     
        
        %load RNN RDM
        load(fullfile(results_avg_dir,rnn_dir,sprintf('ReLU_Layer_%d_Time_%d_Input_RDM.mat',layer-1,t-1)),'data');
        rdm_rnn = data(conds,conds);
        
        %Make sure the diagonal is all 0
        for c1 = 1:numConditions
            for c2 = 1:numConditions
                if c1==c2
                    rdm_rnn(c1,c2) = 0;
                end
            end
        end
        
        %for each bootstrap sample, create a new dataset (N=30) and
        %calculate the peak latency
        for bs = 1:num_bootstrap_samples
            rsa = NaN(num_datasets,numTimepoints);
            for d = 1:num_datasets
                idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
                dataset = squeeze(rdm_subjects(idx,:,:,:));
                rsa(d,:) = representational_SA_rnn(dataset,rdm_rnn);
            end
            avg_rsa = mean(rsa,1);
            %Find peak latency
            peak_latency(index_layer,t,bs) = (find(avg_rsa==max(avg_rsa),1)-40)*5;    
        end
        
        %average over bootstrap samples - not for analysis, just to know
        avg_peak_latency_bs = squeeze(mean(peak_latency,3));     
    end
end

%% 2) Calculate difference between intervals for bootstrap samples averaged over RNN timepoints
peak_latency_avg = mean(peak_latency,2);
layer1_layer4 = abs(squeeze(peak_latency_avg(1,:))-squeeze(peak_latency_avg(2,:)));
layer4_layer7 = abs(squeeze(peak_latency_avg(2,:))-squeeze(peak_latency_avg(3,:)));
layer1_layer7 = abs(squeeze(peak_latency_avg(1,:))-squeeze(peak_latency_avg(3,:)));

%% 3) Get 95% confidence interval for peaks and for peak difference
%peaks
CI_peaks = NaN(numel(layers_idx),2);

for layer = layers_idx
    index_layer = find(layers_idx==layer);
    CI_peaks(index_layer,1) = prctile(squeeze(peak_latency_avg(index_layer,:)),2.5);
    CI_peaks(index_layer,2) = prctile(squeeze(peak_latency_avg(index_layer,:)),97.5);
end

%peak differences
CI_diff_l1_l4 = NaN(1,2);
CI_diff_l4_l7 = NaN(1,2);
CI_diff_l1_l7 = NaN(1,2);

CI_diff_l1_l4(t,1) = prctile(layer1_layer4,2.5);
CI_diff_l1_l4(t,2) = prctile(layer1_layer4,97.5);
CI_diff_l4_l7(t,1) = prctile(layer4_layer7,2.5);
CI_diff_l4_l7(t,2) = prctile(layer4_layer7,97.5);
CI_diff_l1_l7(t,1) = prctile(layer1_layer7,2.5);
CI_diff_l1_l7(t,2) = prctile(layer1_layer7,97.5);        

%% Save as structure
bootstrap_peak_latencies_rsa.peak_latency_bs = avg_peak_latency_bs;
bootstrap_peak_latencies_rsa.peak_latency_true = peak_latency_ground_truth;
bootstrap_peak_latencies_rsa.CI_peaks = CI_peaks;
bootstrap_peak_latencies_rsa.CI_diff_l1_l4 = CI_diff_l1_l4;
bootstrap_peak_latencies_rsa.CI_diff_l4_l7 = CI_diff_l4_l7;
bootstrap_peak_latencies_rsa.CI_diff_l1_l7 = CI_diff_l1_l7;

save(fullfile(results_avg_dir,sprintf('bootstrap_peak_latencies_rsa_%s_subjects_%d_%d',conditions,subjects(1),subjects(end))),'bootstrap_peak_latencies_rsa');

end
