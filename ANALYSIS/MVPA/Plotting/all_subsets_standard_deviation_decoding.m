function all_subsets_standard_deviation_decoding(subjects,subsets) %distance art, distance nat, RT
%ALL_SUBSETS_STANDARD_DEVIATION_DECODING Plot the average baseline standard deviation in the
% of the object decoding results, in different subsets of data.
%
%Input: subjects' ID (e.g., 1:13), subsets (e.g., [0.1 0.2])
%
%Returns a plot of subset versus baseline standard deviation. 
%
%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS';
addpath(genpath(results_dir));

%% Get the standard deviation for all subsets
%Preallocate
numTimepoints = 200;
numConditions = 60;
std_dth = NaN(numel(subsets),numTimepoints);

%Loop over subsets & subjects
for s = 1:numel(subsets)
    decoding_accuracies = NaN(numel(subjects),numConditions,numConditions,numTimepoints);

    for subject = subjects
        subname = get_subject_name(subject);
        load(fullfile(results_dir,subname,sprintf('subset_%s_svm_decoding_accuracy.mat',num2str(subsets(s)))));
        decoding_accuracies(subject,:,:,:) = decodingAccuracy_avg;   
    end    
    
    % Get the standard deviation for all subjects and timepoints
    std_dth(s,:) = std(squeeze(nanmean(nanmean(decoding_accuracies,2),3)));    
end

%Average over baseline timepoints
std_dth_baseline_avg = squeeze(mean(std_dth(:,1:40),2));

%% Plot
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen 
plot(std_dth_baseline_avg,'LineWidth',2);

% %Plot parameters
title(sprintf('Standard deviation in the baseline over different subsets of data  (N=%d)', numel(subjects)));
xlabel('Subset of data (%)')
ylabel('Standard deviation')

%% Save
all_subsets = '_';
for s = 1:numel(subsets)
    all_subsets = strcat(all_subsets,num2str(subsets(s)*100),'_');
end
save_path = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';
save(fullfile(save_path,sprintf('standard_deviation_subsets_%s',all_subsets)),'std_dth');
end