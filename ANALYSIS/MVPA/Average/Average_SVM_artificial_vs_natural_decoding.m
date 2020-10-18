function Average_SVM_artificial_vs_natural_decoding(subjects) 
%AVERAGE_SVM_ARTIFICIAL_VS_NATURAL_DECODING Average over subject-level category decoding
%(artificial scenes vs natural scenes) performed with SVM.  
%
%Input: subject IDs (e.g., 1:4)
%
%Returns a plot of the average category decoding over timepoints. 
%

%Add paths
addpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS/OTHER');
% addpath(genpath('/home/agnek95/SMST/PDM/ANALYSIS/'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS/';
addpath(genpath(results_dir));

%Define variables
numTimepoints = 200;
allSubjects = NaN(numel(subjects),numTimepoints);
legend_cell = cell(numel(subjects) + 1,1);

%Prepare figure
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen

% Loop: create a matrix containing the decoding accuracy matrices of each subject
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,'pseudotrials_svm_artificial_vs_natural_decoding_accuracy.mat'));
    plot(decodingAccuracy_avg,'LineWidth',2);
    hold on;
    legend_cell{subject} = sprintf('Subject %d',subject);
    allSubjects(subject,:) = decodingAccuracy_avg; 
end

avg_decoding_accuracy = squeeze(mean(allSubjects,1)); % avg over subjects

%Plot
plot(avg_decoding_accuracy,'--','LineWidth',2);
title(sprintf('Decoding of artificial vs natural scenes per timepoint (N=%d)',numel(subjects)));
legend_cell{end} = sprintf('Average of %d subjects',numel(subjects));
legend(legend_cell);
xlabel('Timepoint');
ylabel('Decoding accuracy (%)');
leg = legend(legend_cell,'FontSize',12);
set(leg, 'position', [0.7 0.3 0.1 0.01]); %put legend below graph

%Save
results_avg = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';
saveas(gcf,fullfile(results_avg,sprintf('svm_artificial_vs_natural_decoding_accuracy_%d_subjects',numel(subjects))));
saveas(gcf,fullfile(results_avg,sprintf('svm_artificial_vs_natural_decoding_accuracy_%d_subjects.svg',numel(subjects))));

end