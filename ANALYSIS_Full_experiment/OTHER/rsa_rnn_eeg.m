function rsa_rnn_eeg(subjects,with_stats) 
%RSA_RNN_EEG Plots the RSA (already computed) between RNN and EEG with stats, if desired.
%
%Input: subjects' ID (e.g., 1:13), add stats (1 with, 0 without)
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS'));
addpath(genpath('/scratch/agnek95/PDM/DATA/RNN_RTs'));
addpath(genpath('/home/agnek95/PyColormap4Matlab'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%% Load the RDMs
%Define some variables
num_conditions = 60;
num_timepoints_rnn = 8;
layer_idxs = [1,5,7];
legend_bool = 0;

%Stats parameters
num_permutations = 1000;
tail = 'right';
q_value = 0.05;

%Load the EEG RDM
load(fullfile(results_avg_dir,sprintf('average_rdm_categorization_subjects_%d_%d',...
    subjects(1),subjects(end))),'rdm_cat');
rdm_eeg = rdm_cat;

%Load the RNN RDMs
rdm_rnn = NaN(numel(layer_idxs),num_timepoints_rnn,num_conditions,num_conditions);
rdm_dir = fullfile(results_avg_dir,'Input_RDM_RNN');
for l=layer_idxs
    for t=1:num_timepoints_rnn
        load(fullfile(rdm_dir,sprintf('ReLU_Layer_%d_Time_%d_Input_RDM_PCA.mat',...
            l-1,t-1)),'data');
        rdm_rnn(l,t,:,:) = data;
        
        %Make sure the diagonal is all 0
        for c1 = 1:num_conditions
            for c2 = 1:num_conditions
                if c1==c2
                    rdm_rnn(:,:,c1,c1) = 0;
                end
            end
        end
    end
end

%% Load and plot the RSA results
rsa_results_dir = fullfile(results_avg_dir,'RSA_matrices_eltanin');
plot_location = -0.1:-0.01:-0.17;

for c = 1:3 %natural,artificial,all
    if c == 1
        conditions = 'artificial';
        conds = 1:30;
        colormap_plot = 'PuRd';
    elseif c == 2
        conditions = 'natural';
        conds = 31:60;
        colormap_plot = 'YlGn';
    elseif c == 3
        conditions = 'all';
        conds = 1:60;
        colormap_plot = 'OrRd';
    end
    filename = sprintf('Model_RDM_PCA_7_layers_8_timepoints_%s.mat',conditions);
    load(fullfile(rsa_results_dir,filename),'data');
    rsa_results = data;
    
    cmap = getPyPlot_cMap(colormap_plot,num_timepoints_rnn);
    for l=layer_idxs
        figure;
        legend_plot = cell(num_timepoints_rnn,1);

        for t=1:num_timepoints_rnn
            %Plot the results
            plot(squeeze(rsa_results(l,t,:)),'LineWidth',2,'Color',cmap(t,:));
            hold on;
            legend_plot{t} = sprintf('Timepoint %s',num2str(t));
            
            %Stats
            if with_stats
                %Check if exist
                filename_sign = 'rsa_rnn_eeg_stats';   
                filename = fullfile(results_avg_dir,sprintf('%s_%s_subjects_%d_%d_layer_%d_tp_%d.mat',...
                    filename_sign,conditions,subjects(1),subjects(end),l,t));
                if exist(filename,'file')
                    load(filename,'rsa_rnn_eeg_stats');
                else
                    [rsa_rnn_eeg_stats.SignificantVariables, rsa_rnn_eeg_stats.pvalues, rsa_rnn_eeg_stats.crit_p,...
                        rsa_rnn_eeg_stats.adjusted_pvalues] = fdr_rsa_rnn_eeg(rdm_eeg(conds,conds,:),...
                        squeeze(rdm_rnn(l,t,conds,conds)),rsa_results(l,t,:),num_permutations,tail,q_value);
                    save(filename,'rsa_rnn_eeg_stats');
                end
                
                %Plot the stats
                st = (rsa_rnn_eeg_stats.SignificantVariables*plot_location(t)); %depending on the stats
                st(st==0) = NaN;
                plot(st,'*','Color',cmap(t,:)); 
                hold on;
            end
        end
        if legend_bool==1
            legend(legend_plot,'Location','best');
        end
        ylim([-0.2 0.6]);
        font_size = 18;
        set(gca,'FontName','Arial','FontSize',font_size);
        xticks(0:20:200);
        xticklabels(-200:100:800);        
        onset_time = 40;
        xline(onset_time,'--');
        
        %Save plot
        filename_plot = fullfile(results_avg_dir,sprintf('%s_%s_subjects_%d_%d_layer_%d',...
            filename_sign,conditions,subjects(1),subjects(end),l));
        saveas(gcf,sprintf('%s.svg',filename_plot)); 
        saveas(gcf,sprintf('%s.fig',filename_plot)); 
    end
end

end





