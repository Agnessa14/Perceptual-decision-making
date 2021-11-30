function rsa_rnn_eeg_subject_by_subject(subjects,with_stats) 
%RSA_RNN_EEG_SUBJECT_BY_SUBJECT Compute and plots the RSA between RNN and EEG with stats, if desired.
%RSA is computed by correlating each of the subject-specific EEG RDMs with
%the RNN RDM and averaging over.
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
numTimepointsRNN = 8;
numTimepointsEEG = 200;
layers_idx = [1,4,7];
legend_bool = 0;

%Stats parameters
% num_permutations = 1000;
% tail = 'right';
% q_value = 0.05;
% permutation_stats.num_perms = 1000;
% permutation_stats.cluster_th = 0.05;
% permutation_stats.significance_th = 0.05;

%Load subject-level RDMs 
numConditions = 60;
rdm_eeg_all_subjects = NaN(max(subjects),numConditions,numConditions,numTimepointsEEG); 

for subject = subjects
    subname = get_subject_name(subject);
    subject_results_dir = fullfile(results_dir,subname);
    load(fullfile(subject_results_dir,'rdm_pearson_categorization.mat'),'rdm_avg');           
    rdm_eeg_all_subjects(subject,:,:,:) = rdm_avg;      
end   

%Prepare for loading the RNN Input RDMs
rdm_dir = fullfile(results_avg_dir,'02.11_2_rnn/Input_RDM');


%% Calculate and plot the RSA results
% rsa_results_dir = fullfile(results_avg_dir,'02.11_2_rnn/Model_RDM_redone');
rsa_results = NaN(max(subjects),3,numel(layers_idx),numTimepointsRNN,numTimepointsEEG); %artificial, natural, all
plot_location = -0.1:-0.02:-0.24;
model_name = '02_11_2';

for c = 1:3 %artificial,natural,all
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
    
    %get colormap from Python
    cmap = getPyPlot_cMap(colormap_plot,numTimepointsRNN);
         
    %load noise ceilings (lower and upper bound)
    filename_ceiling_l = sprintf('lower_noise_ceiling_rsa_subjects_%d_%d_%s_scenes',subjects(1),subjects(end),conditions);
    load(fullfile(results_avg_dir,filename_ceiling_l),'noise_ceiling_lower_bound');    
    filename_ceiling_u = sprintf('upper_noise_ceiling_rsa_subjects_%d_%d_%s_scenes',subjects(1),subjects(end),conditions);
    load(fullfile(results_avg_dir,filename_ceiling_u),'noise_ceiling_upper_bound');
   
    
    for l=layers_idx
        %plot the noise ceilings as a shaded area
        figure;
        dark_grey = [0.5 0.5 0.5];
        light_grey = [0.8 0.8 0.8];
        x_var = [1:numTimepointsEEG, fliplr(1:numTimepointsEEG)];
        shaded_area = [noise_ceiling_upper_bound, fliplr(noise_ceiling_lower_bound)];
        fill(x_var, shaded_area, light_grey,'FaceAlpha',0.3);
        hold on;
        plot(noise_ceiling_lower_bound,'LineWidth',2,'Color',dark_grey);
        plot(noise_ceiling_upper_bound,'LineWidth',2,'Color',dark_grey);
        hold on;
        legend_plot = cell(numTimepointsRNN,1);

        for t=1:numTimepointsRNN
            
            %load the RNN RDM
            load(fullfile(rdm_dir,sprintf('ReLU_Layer_%d_Time_%d_Input_RDM.mat',...
                l-1,t-1)),'data');
            rdm_rnn = data(conds,conds);
            
            %Make sure the diagonal is all 0
            for c1 = 1:numel(conds)
                for c2 = 1:numel(conds)
                    if c1==c2
                        rdm_rnn(c1,c2) = 0;
                    end
                end
            end
            
            %RSA: calculate correlation between RNN and EEG
            for subject = subjects
                rdm_eeg = squeeze(rdm_eeg_all_subjects(subject,conds,conds,:));
                rsa_results(subject,c,l,t,:) = representational_SA_rnn(rdm_eeg,rdm_rnn); %modify the RSA function 
            end  
            
            %plot the RSA averaged over subjects
            rsa_results = rsa_results(subjects,:,:,:,:);
            plot(squeeze(mean(rsa_results(:,c,l,t,:),1)),'LineWidth',2,'Color',cmap(t,:));
            hold on;
            
            %Stats: FDR-corrected
            if with_stats
                
                %Check if exist
                filename_sign = sprintf('%s_rsa_rnn_eeg_stats',model_name);   
                filename = fullfile(results_avg_dir,sprintf('%s_%s_subjects_%d_%d_layer_%d_tp_%d.mat',...
                    filename_sign,conditions,subjects(1),subjects(end),l,t));
                if exist(filename,'file')
                    load(filename,'fdr_stats');
                else %if not, run them
%                     [rsa_rnn_eeg_stats.SignificantVariables, rsa_rnn_eeg_stats.pvalues, rsa_rnn_eeg_stats.crit_p,...
%                         rsa_rnn_eeg_stats.adjusted_pvalues] = fdr_rsa_rnn_eeg(rdm_eeg(conds,conds,:),...
%                         squeeze(rdm_rnn(l,t,conds,conds)),rsa_results(l,t,:),num_permutations,tail,q_value);
%                     save(filename,'rsa_rnn_eeg_stats');
                    fdr_stats.num_perms = 1000;
                    fdr_stats.tail = 'right';
                    fdr_stats.qvalue = 0.05;
%                     [fdr_stats.significant_timepoints,fdr_stats.pvalues,...
%                         fdr_stats.crit_p, fdr_stats.adjusted_pvalues]...
%                         = fdr_permutation_cluster_1sample_alld(for_stats_data,...
%                         fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);
                     [fdr_stats.significant_timepoints,fdr_stats.pvalues,...
                        fdr_stats.crit_p, fdr_stats.adjusted_pvalues]...
                        = fdr_rsa_random_effects_stats(rdm_eeg_all_subjects,rdm_rnn,...
                        fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);
                    
                    save(filename,'fdr_stats');
                end
                
                %Plot the stats
                st = (fdr_stats.SignificantVariables*plot_location(t)); %depending on the stats
                st(st==0) = NaN;
                plot(st,'*','Color',cmap(t,:)); 
                hold on;
            end
        end
        
        %Plot parameters
        if legend_bool==1
            legend(legend_plot,'Location','best');
        end
        ylim([-0.3 0.8]);
        yticks(-0.2:0.2:0.8);
        font_size = 18;
        set(gca,'FontName','Arial','FontSize',font_size);
        xticks(0:20:200);
        xticklabels(-200:100:800);        
        onset_time = 40;
        xline(onset_time,'--');
        
        %Save plot
%         filename_part = 'rsa_plus_noise_ceilings_subject_level';
%         filename_plot = fullfile(results_avg_dir,sprintf('%s_%s_%s_subjects_%d_%d_layer_%d',...
%             model_name,filename_part,conditions,subjects(1),subjects(end),l));
%         saveas(gcf,sprintf('%s.svg',filename_plot)); 
%         saveas(gcf,sprintf('%s.fig',filename_plot)); 
    end
end

end





