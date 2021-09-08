function rnn_dth_one_layer
%RNN_DTH_ONE_LAYER Performs distance-to-hyperplane  using
%the distances-to-hyperplane using RNN dats. Plots the last layer.
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));

%% DTH
numConditions = 60;
artificial_conditions = 1:numConditions/2;
natural_conditions = (numConditions/2)+1:numConditions;  
entropy_thresh = '0.09';

%Load RTs
load(sprintf('/scratch/agnek95/PDM/DATA/RNN_RTs/reaction_time_entropy_th_%s_model_22.08.mat',...
    entropy_thresh),'data');
RT = data';
RT_art = RT(artificial_conditions);
RT_nat = RT(natural_conditions);
    
%load distances
load('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS/rnn_distances_all_scenes_model_22.08.mat','data');
distances = data;
     
%% Correlate each subject's distances with the median RT
num_layers=size(distances,1);
num_timepoints_rnn=size(distances,2);
size_corr = [num_layers,num_timepoints_rnn];
correlation_art = NaN(size_corr);
correlation_nat = NaN(size_corr);
correlation_both = NaN(size_corr);
all_conditions = [artificial_conditions natural_conditions];

for l = 1:num_layers
    for tp = 1:num_timepoints_rnn
        correlation_art(l,:) = arrayfun(@(x) ...
            corr(squeeze(distances(l,x,artificial_conditions)),...
            RT_art,'type','Spearman'), 1:num_timepoints_rnn);
        correlation_nat(l,:) = arrayfun(@(x) ...
            corr(squeeze(distances(l,x,natural_conditions)),...
            RT_nat,'type','Spearman'), 1:num_timepoints_rnn);
        correlation_both(l,:) = arrayfun(@(x) ...
            corr(squeeze(distances(l,x,all_conditions)),...
            RT,'type','Spearman'), 1:num_timepoints_rnn);
    end
end


%% Surface plot
for p = [1,3] %1:3
    if p==1
        conditions = 'all';
        corr_plot = correlation_both;
    elseif p==2
        conditions = 'artificial';
        corr_plot = correlation_art;
    elseif p==3
        conditions = 'natural';
        corr_plot = correlation_nat;
    end
    surf(corr_plot);
    colorbar;
    title(sprintf('Correlation between RT and distances-to-the-hyperplane in an RNN (%s scenes)',conditions));
    xlabel('Timepoint');
    ylabel('Layer');
    zlabel('Spearman''s r');
    saveas(gcf,fullfile('SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/',sprintf('rnn_dth_surface_plot_%s.png',conditions)));
    close(gcf);
end







