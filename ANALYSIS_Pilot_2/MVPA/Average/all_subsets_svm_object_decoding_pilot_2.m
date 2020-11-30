function all_subsets_svm_object_decoding_pilot_2(subjects,subsets) %distance art, distance nat, RT
%ALL_SUBSETS_SVM_OBJECT_DECODING_PILOT_2 Plots the averaged (over participantd) object decoding
%performed on different subsets of data. 
%
%Input: subjects' ID (e.g., 1:13), subsets (e.g., [0.1 0.3])
%
%Takes the object decoding for each participant and averages it, for each subset, 
%resulting in N curves, N being the number of subsets. 
%
%% Setup 
% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS'));
results_dir = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS';
addpath(genpath(results_dir));

%Setup the subset names
% subset_str = cell(numel(subsets),1);
% 
% for s = 1:numel(subsets)
%     subset_str{s} = num2str(subsets(s));
% end

%Setup some variables
numTimepoints = 200;
numConditions = 60;
legend_cell = cell(numel(subsets),1);

%Setup figure
figure(abs(round(randn*10)));
set(gcf, 'Position', get(0, 'Screensize')); %make fullscreen

%% Plot the decoding accuracies from all subjects & all subsets
%Collect the DAs
for s = 1:numel(subsets)
    decoding_accuracies_all_subjects = NaN(numel(subjects),numConditions,numConditions,numTimepoints);

    for subject = subjects
        subname = get_subject_name(subject);
        load(fullfile(results_dir,subname,sprintf('subset_%s_svm_decoding_accuracy.mat',num2str(subsets(s)))));
        decoding_accuracies_all_subjects(subject,:,:,:) = decodingAccuracy_avg;   
    end  
    
    avg_over_conditions_all_subjects = squeeze(nanmean(decoding_accuracies_all_subjects,1:3));
    
    %Plot
    plot(avg_over_conditions_all_subjects,'LineWidth',2);
    hold on;
    legend_cell{s} = sprintf('%s percent ',num2str(subsets(s)*100));
end
    
%Plot parameters
title(sprintf('Scene decoding per timepoint for various subsets of data (N=%d)', numel(subjects)));
xline(40,'--');
legend_cell{numel(subsets)+1} = 'Stimulus onset';
l = legend(legend_cell,'FontSize',12);
% set(l, 'position', [0.7 0.2 0.1 0.01]); %put legend below graph
% axis([0,200,-0.8,0.8]);
xticks(0:10:200)
xlabel('Timepoint')
ylabel('Decoding accuracy (%)')

%% Save figure
all_subsets = '_';
for s = 1:numel(subsets)
    all_subsets = strcat(all_subsets,num2str(subsets(s)*100),'_');
end
save_path = '/home/agnek95/SMST/PDM_PILOT_2/RESULTS_AVG/';
saveas(gcf,fullfile(save_path,sprintf('subset%sSVM_object_decoding_%d_subjects',all_subsets,numel(subjects)))); %turn subset into percentage: otherwise thinks it's a file extension
saveas(gcf,fullfile(save_path,sprintf('subset%sSVM_object_decoding_%d_subjects.svg',all_subsets,numel(subjects))));

end