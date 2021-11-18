function bootstrap_peak_latency_rsa(subjects)
%BOOTSTRAP_PEAK_LATENCY_RSA Apply bootstrapping to calculate the 95% confidence
%interval of peak latency difference (for RSA between RNN and EEG-categorization task).
%
%Input: subject IDs (e.g., 1:13)
%
%Output: saved boostrapped peak latencies (in ms), 95% confidence interval of the difference
%
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS/';

                                    %%%%% SETUP %%%%%
%% Pre-allocate + setup for loading
% numTimepoints = 200;
% decoding_accuracies_all = NaN(sorted_subjects(end),numTimepoints);
% if strcmp(analysis,'object_decoding')
%     filename = 'svm_decoding_accuracy';
% elseif strcmp(analysis,'category_decoding')
%     filename = 'svm_artificial_vs_natural_decoding_accuracy';
% end
% 
%% Load the subject-level RDMs
numConditions = 60;
numTimepoints = 200;
layers_idx = [1,4,7]; %RNN layers of interest
numTimepointsRNN = 8;
rdm_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepoints); 

for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'rdm_pearson_categorization.mat'),'rdm_avg');           
    rdm_all_subjects(subject,:,:,:) = rdm_avg;      
end   

rdm_subjects = rdm_all_subjects(subjects,:,:,:); %remove NaN subjects
                                  %%%%% BOOTSTRAPPING %%%%%

                                  
%% 0) Prepare required variables
num_bootstrap_samples = 1000;
num_datasets = size(rdm_subjects,1);

%% 1) Create the bootstrap samples & calculate the peak difference
rng(0,'twister');
max_dataset = num_datasets; %for the random dataset index generator
min_dataset = 1;
rnn_dir = '02.11_2_rnn/Input_RDM';

for layer = layers_idx
    for t = 1:numTimepointsRNN
        
        %load RNN RDM
        load(fullfile(rnn_dir,sprintf('ReLU_Layer_%d_Time_%d.mat',layer,t)),'data');
        rdm_rnn = data;
        %Make sure the diagonal is all 0
        for c1 = 1:num_conditions
            for c2 = 1:num_conditions
                if c1==c2
                    rdm_rnn(:,:,c1,c1) = 0;
                end
            end
        end
%         rdm_rnn(isnan(rdm_1f)) = 0;
%         rdm_dnn_flattened_cell = squareform(avg_datasets(:,:,x)+(avg_datasets(:,:,x))'),...
%             1:numTimepoints,'UniformOutput',false);
%         rdm_eeg_flattened_1 = reshape(cell2mat(rdm_flattened_cell_1),[],numTimepoints);

        %create the samples from the EEG RDMs
        peak_latency_1 = NaN(num_bootrstrap_samples,1);
        peak_latency_2 = NaN(num_bootrstrap_samples,1);
        peak_latency_diff = NaN(num_bootrstrap_samples,1);

        for bs = 1:num_bootstrap_samples
            datasets = NaN(size(rdm_subjects));
            for d = 1:num_datasets
                idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
                datasets(d,:,:,:) = rdm_subjects(idx,:,:,:);
            end
            avg_datasets_rdm_eeg = squeeze(mean(datasets,1));
            
            %Redo the RSA on the average of the datasets of the sample,
            %then calculate the peak latencies
%             avg_datasets(isnan(rdm_1f)) = 0;
%             rdm_eeg_flattened_cell = arrayfun(@(x) squareform(avg_datasets(:,:,x)+(avg_datasets(:,:,x))'),...
%                 1:numTimepoints,'UniformOutput',false);
%             rdm_eeg_flattened = reshape(cell2mat(rdm_flattened_cell_1),[],numTimepoints);
%                 peak_latency_all_samples(bs) = find(avg_datasets==max(avg_datasets),1);
            rsa = representational_SA_rnn(avg_datasets_rdm_eeg,rdm_rnn);        
            %Find both peak latencies
            corr_sorted = sort(rsa,'descend');
            
            %peak latency should not be at the end of the trial
%             %peak latency 1
%             if find(rsa==max(rsa)) < 500
%                 peak_latency_1 = find(rsa==max(rsa));
%             else
%                 if find(rsa==corr_sorted(2)) < 500
%                     peak_latency_1 = find(rsa==corr_sorted(2));
%                 else
%                     error('Something''s wrong');
%                 end
%             end
%             
%             %peak latency 2
%             if find(rsa==corr_sorted(2)) < 500 && find(rsa==corr_sorted(2)) ~= peak_latency_1
%                 peak_latency_2 = find(rsa==corr_sorted(2));
%             else
%                 if find(rsa==corr_sorted(3)) < 500
%                     peak_latency_2 = corr_sorted(3);
%                 else
%                     error('Something''s wrong');
%                 end
%             end
            
 
            
            i = 1;
            while (find(rsa==corr_sorted(i)) - 40)*5 >= 500
                i = i+1;
            end
            peak_latency_1(bs) = find(rsa==corr_sorted(i));
            
            j = 2;
            while (find(rsa==corr_sorted(j)) - 40)*5 >= 500
                j = j+1;
            end
            peak_latency_2(bs) = (find(rsa==corr_sorted(j)) - 40)*5;
            
            %Calculate peak difference
            peak_latency_diff(bs) = abs(peak_latency_1-peak_latency_2);
            
            
        end
    end
end


%% 2) Get 95% confidence interval for the peak difference
confidence_interval = NaN(2,1);
confidence_interval(1) = prctile(peak_latency_diff,2.5);
confidence_interval(2) = prctile(peak_latency_diff,97.5);a

%% Bonus: Get mean bootstrapped peak latencies
peak_latency_1_avg = round(mean(peak_latency_1));
peak_latency_2_avg = round(mean(peak_layency_2));
peak_latency_diff_avg = round(mean(peak_layency_diff));

%% Save as structure
bootstrap_peak_latencies_rsa.peak_latency_1 = peak_latency_1_avg;
bootstrap_peak_latencies_rsa.peak_latency_2 = peak_latency_2_avg;
bootstrap_peak_latencies_rsa.peak_latency_diff = peak_latency_diff_avg;
bootstrap_peak_latencies_rsa.CI = confidence_interval;
save(fullfile(results_avg_dir,sprintf('bootstrap_peak_latencies_rsa_subjects_%d_%d',subjects(1),subjects(end))),'bootstrap_peak_latencies_rsa');

end