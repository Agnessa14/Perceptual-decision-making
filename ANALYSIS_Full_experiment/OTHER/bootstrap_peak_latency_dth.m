function bootstrap_peak_latency_dth(task_distance,task_RT)
%BOOTSTRAP_PEAK_LATENCY Apply bootstrapping to calculate the 95% confidence
%interval of peak distance-to-hyperplane latency.
%
%Input: task of EEG distances (1=categorization,2=distraction), task of RT
%(same), scene conditions ('artificial','natural' or 'both')
%
%

%% Paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';
                                    %%%%% SETUP %%%%%
%% Load all subjects' data
task_name_distance = get_task_name(task_distance);

for c = 1:3 %artificial, natural, both
    switch c
        case 1
            conditions = 'artificial';
            conds_short = conditions(1:3);
        case 2 
            conditions = 'natural';
            conds_short = conditions(1:3);
        case 3
            conditions = 'both';
            conds_short = conditions;
    end
    filename = sprintf('all_subjects_dth_%s_%s',task_name_distance,conditions);
    if task_distance~=task_RT
        filename = sprintf('all_subjects_dth_%s_distances_%s',task_name_distance,conditions);
    end
    load(fullfile(results_avg,filename));
    all_sub_dth = eval(sprintf('all_sub_dth_%s',conds_short));
     
    %%%%% BOOTSTRAPPING %%%%%

    %% 0) Prepare required variables
    num_bootstrap_samples = 1000;
    num_datasets = size(all_sub_dth,1); %num subjects
    peak_latency_all_samples = NaN(num_bootstrap_samples,1);

    %% 1) Create the bootstrap samples
    rng(0,'twister');
    max_dataset = num_datasets; %for the random dataset index generator
    min_dataset = 1;
    for bs = 1:num_bootstrap_samples
        datasets = NaN(size(all_sub_dth));
        for d = 1:num_datasets
            idx = round((max_dataset-min_dataset).*rand(1,1) + min_dataset); %pick one random number between one and num_datasets
            datasets(d,:) = all_sub_dth(idx,:);
        end
        avg_datasets = squeeze(mean(datasets,1));
        peak_latency_all_samples(bs) = find(avg_datasets==min(avg_datasets),1);
    end

    %% 2) Get ground mean bootstrapped peak latency & max value *actually min bcs we are interested in negative correlations*
    mean_all_sub_dth = mean(all_sub_dth,1);
    bootstrap_dth.peak_value_ground = min(mean_all_sub_dth);
    bootstrap_dth.peak_latency_ground = (find(mean_all_sub_dth==min(mean_all_sub_dth))-40)*5;
    
    %% 3) Get 95% confidence interval
    bootstrap_dth.confidence_interval = NaN(2,1);
    bootstrap_dth.confidence_interval(1) = (prctile(peak_latency_all_samples,2.5)-40)*5;
    bootstrap_dth.confidence_interval(2) = (prctile(peak_latency_all_samples,97.5)-40)*5;
    fprintf('Peak value: %d \nPeak latency: %d \nConfidence interval: [%d %d]\n',...
        bootstrap_dth.peak_value_ground,bootstrap_dth.peak_latency_ground,bootstrap_dth.confidence_interval(1),bootstrap_dth.confidence_interval(2)); 
    %% Save
    save(fullfile(results_avg,sprintf('bootstrap_%s',filename')),'bootstrap_dth');
end


end



