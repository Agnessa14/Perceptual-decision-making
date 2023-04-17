function rsa_feedforward_eeg_subject_by_subject(subjects,model_type,with_stats,with_error_bars,rsa_in_matlab) 
%RSA_FEEDFORWARD_EEG_SUBJECT_BY_SUBJECT Compute and plot the RSA between DNN and EEG with stats, if desired.
%RSA is computed by correlating each of the subject-specific EEG RDMs with
%the RNN RDM and averaging over. Average over 3 layers
%
%Input: subjects' ID (e.g., 1:13), model_type ('b_d' or 'b'), add stats (1 with, 0 without), plot
%with/without error bars (1/0),rsa_in_matlab (1 for performed in Matlab, 0
%for Python)
%
%Author: Agnessa Karapetian, 2021

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
addpath(genpath('/home/agnek95/PyColormap4Matlab'));
results_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS';
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
addpath(genpath(results_dir));

%% Load the RDMs
%Define some variables
numTimepointsEEG = 200;
if strcmp(model_type,'b_d')
    layers_idx = [1,7,14];
elseif strcmp(model_type,'b')
    layers_idx = [1,4,7];
end

legend_bool = 0;

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
if strcmp(model_type,'b_d')
    rdm_dir = fullfile(results_avg_dir,'B_D_net/Input_RDM_pearson');
elseif strcmp(model_type,'b')
    rdm_dir = fullfile(results_avg_dir,'B_net/Input_RDM_pearson');
end
    
%% Calculate and plot the RSA results
plot_location = -0.1:-0.04:-0.2;
% model_name = '02_11_2';
plot_whole_epoch = 0;
rsa_results = NaN(numel(layers_idx),numTimepointsEEG); %artificial, natural, all

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
    cmap = getPyPlot_cMap(colormap_plot,10); %randomly generate 10 colours
         
    %load noise ceilings (lower and upper bound)
    filename_ceiling_l = sprintf('lower_noise_ceiling_rsa_subjects_%d_%d_%s_scenes',subjects(1),subjects(end),conditions);
    load(fullfile(results_avg_dir,filename_ceiling_l),'noise_ceiling_lower_bound');    
    filename_ceiling_u = sprintf('upper_noise_ceiling_rsa_subjects_%d_%d_%s_scenes',subjects(1),subjects(end),conditions);
    load(fullfile(results_avg_dir,filename_ceiling_u),'noise_ceiling_upper_bound');  
    figure;
    
    %preallocate
    rsa_all_layers = NaN(max(subjects),numel(layers_idx),numTimepointsEEG);

    for l=layers_idx
        index_layer = layers_idx==l;
        color = cmap(find(index_layer==1)*3,:);
        %plot the noise ceilings as a shaded area
        dark_grey = [0.5 0.5 0.5];
        light_grey = [0.8 0.8 0.8];
        x_var = [1:numTimepointsEEG, fliplr(1:numTimepointsEEG)];
        shaded_area = [noise_ceiling_upper_bound, fliplr(noise_ceiling_lower_bound)];
        fill(x_var, shaded_area, light_grey,'FaceAlpha',0.3);
        hold on;
        plot(noise_ceiling_lower_bound,'LineWidth',2,'Color',dark_grey);
        plot(noise_ceiling_upper_bound,'LineWidth',2,'Color',dark_grey);
        hold on;
        %legend_plot = cell(numTimepointsDNN,1);

        %load the DNN RDM
        if rsa_in_matlab
            load(fullfile(rdm_dir,sprintf('ReLU_Layer_%d_Input_RDM_pearson.mat',...
                l-1)),'input_rdm');
            rdm_dnn = input_rdm(conds,conds);
        else
            load(fullfile(rdm_dir,sprintf('ReLU_Layer_%d_Input_RDM_pearson.mat',...
                l-1)),'data');
            rdm_dnn = data(conds,conds);
        end
        %Make sure the diagonal is all 0
        for c1 = 1:numel(conds)
            for c2 = 1:numel(conds)
                if c1==c2
                    rdm_dnn(c1,c2) = 0;
                end
            end
        end

        %RSA: calculate correlation between RNN and EEG
        rsa_results_layer = NaN(max(subjects),numTimepointsEEG);
        for subject = subjects
            rdm_eeg = squeeze(rdm_eeg_all_subjects(subject,conds,conds,:));
            results_rdm = representational_SA_rnn(rdm_eeg,rdm_dnn); %modify the RSA function 
            rsa_results_layer(subject,:) = results_rdm;
            rsa_all_layers(subject,index_layer,:) = results_rdm;
        end  
  
        %Average over subjects & plot
        rsa_results_layer=rsa_results_layer(subjects,:);
        rsa_results_avg = squeeze(mean(rsa_results_layer,1));
%         rsa_results(index_layer,:) = rsa_results_avg;
        if ~with_stats
            plot(rsa_results_avg,'LineWidth',2,'Color',color);
            hold on;
        
        %Stats: FDR-corrected
        else 
            %Check if exist
            filename = fullfile(results_avg_dir,sprintf('%s_net_rsa_stats_%s_subjects_%d_%d_layer_%d.mat',...
                model_type,conditions,subjects(1),subjects(end),l));
            if exist(filename,'file')
                load(filename,'fdr_stats');
            else %if not, run them
                fdr_stats.num_perms = 10000;
                fdr_stats.tail = 'right';
                fdr_stats.qvalue = 0.05;
                [fdr_stats.significant_timepoints,fdr_stats.crit_p,fdr_stats.adjusted_pvalues] = ...
                    fdr_rsa_random_effects_stats(rsa_results_layer,fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);                   
                save(filename,'fdr_stats');
            end

            %Plot the stats
            st = (fdr_stats.significant_timepoints*plot_location(index_layer)); 
            st(st==0) = NaN;
            plot(st,'*','Color',cmap(find(index_layer==1)*3,:)); 
            hold on;
            
            %if needed: error bars
            if with_error_bars
                %calculate error bars
                stdDM = std(rsa_results_layer); 
                err = stdDM/sqrt(size(rsa_results_layer,1)); %standard deviation/sqrt of num subjects  

                %plot as a shaded area
                top_curve = rsa_results_avg + err;
                bottom_curve = rsa_results_avg - err;
                x2 = [1:numTimepointsEEG, fliplr(1:numTimepointsEEG)];
                shaded_area = [top_curve, fliplr(bottom_curve)];
                fill(x2, shaded_area, color,'FaceAlpha',0.2);
                hold on;
            end
            
            %Plot the data
            plot(rsa_results_avg,'LineWidth',2,'Color',color);
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

    if ~plot_whole_epoch 
        xlim([20 200]);
    end

    %Save plot
    filename_part = sprintf('rsa_%s_net',model_type);
    if ~plot_whole_epoch 
        filename_part = sprintf('%s_-100ms',filename_part);
    end
    if with_error_bars
        filename_part = sprintf('error_bars_%s',filename_part);
    end
    filename_plot = fullfile(results_avg_dir,sprintf('%s_%s_subjects_%d_%d',...
        filename_part,conditions,subjects(1),subjects(end)));
    saveas(gcf,sprintf('%s.svg',filename_plot)); 
    saveas(gcf,sprintf('%s.fig',filename_plot)); 
    
    save(sprintf('%s',filename_plot),'rsa_results');
    rsa_all_subs = rsa_all_layers(subjects,:,:);
    save(fullfile(results_avg_dir,sprintf('all_subjects_all_layers_rsa_%s_net_%s',model_type,conditions)),'rsa_all_subs');
end

end


