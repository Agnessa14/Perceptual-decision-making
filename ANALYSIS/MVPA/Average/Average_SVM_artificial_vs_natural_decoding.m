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
addpath(genpath('/home/agnek95/SMST/PDM/ANALYSIS/'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS/';
addpath(genpath(results_dir));

%Define variables
numTimepoints = 200;
allSubjects = NaN(numel(subjects),numTimepoints);

% Loop: create a matrix containing the decoding accuracy matrices of each subject
for subject = subjects
    subname = get_subject_name(subject);
    load(fullfile(results_dir,subname,'svm_artificial_vs_natural_decoding_accuracy.mat'));
    allSubjects(subject,:) = decodingAccuracy_Avg; 
end


end