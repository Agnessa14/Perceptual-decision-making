function plot_RT_rnn_vs_human(model_type)
%PLOT_RT_RNN_VS_HUMAN Create a scatterplot of reaction times.
%
%Returns a scatterplot of human reaction times in the categorization task vs
%RNN reaction times. 
%
%Input: model_type ('b', b_d' or 'bl')
%
%Author: Agnessa Karapetian, 2021
%

%% Add paths
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/';

%% Load the correlations and the RTs for the noise ceiling calculation
load(fullfile(results_avg_dir,'RT_all_subjects_5_35_categorization.mat'),'RTs');

%BLnet
load(fullfile(results_avg_dir,'02.11_2_rnn/Model_RDM_redone','correlation_RT_human_RNN_cross-validated.mat'),'data');
correlations_RT_BLnet=data;

%Feedforward model
if strcmp(model_type,'b_d')
    load(fullfile('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/DNN/','correlation_RT_human_b_d_net_7_layers_cross-validated_readouts.mat'),'data');
elseif strcmp(model_type,'b')
    load(fullfile('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/DNN/','correlation_RT_human_b_net_cross-validated_readouts.mat'),'data');
end

correlations_RT = data;
num_subjects = size(correlations_RT,1);
RT_noise_ceiling_lower = NaN(3,1);
xticklabels = cell(3,1);
h = NaN(2,3);
location_x = NaN(3,1);
location_x_diff = NaN(3,1);
label_x = NaN(3,1);
label_x_diff = NaN(3,1);


for pl=1:2 %1=corr results, 2=diff. plot 
    figure;

    for c = 1:3
        if c == 1
            xticklabels{c} = 'All scenes';
            conds = 1:60;
            location_x(c) = 0.5;
            label_x(c) = 0.5;
            location_x_diff(c) = 2.0;
            label_x_diff(c) = 2.0;
        elseif c==2
            xticklabels{c} = 'Man-made scenes';
            conds = 1:30;
            location_x(c) = 1.5;
            label_x(3) = 1.5;
            location_x_diff(c) = 3.0;
            label_x_diff(c) = 3.0;
        elseif c==3
            xticklabels{c} = 'Natural scenes';
            conds = 31:60;
            location_x(c) = 1;
            label_x(2) = 1;
            location_x_diff(c) = 2.5;
            label_x_diff(c) = 2.5;
        end

        corr_conditions = correlations_RT(:,c);
        average_corr = mean(corr_conditions,1);

        %Calculate noise ceiling
        if pl==1
            noise_ceil_temp = NaN(1,num_subjects);
            for subject = 1:num_subjects
                conds_bool = ~isnan(RTs(subject,conds));
                noise_ceil_temp(subject) = corr(nanmean(RTs(1:end~=subject,conds(conds_bool)),1)',RTs(subject,conds(conds_bool))','type','Pearson');
            end
            RT_noise_ceiling_lower(c) = mean(noise_ceil_temp);
        end

        %Plot all subjects, average and noise ceiling
        if pl==1
            color_subjects = [0.75 0.75 0.75];   
            color_avg = [0 0.5 0.5];
            color_noise_ceiling = 'k';
            plot_data_sub=corr_conditions;
            plot_avg=average_corr;
        elseif pl==2 
            color_subjects = [0.35 0.35 0.35];
            color_avg = 'y';
            diff_models=correlations_RT_BLnet-correlations_RT;
            plot_data_sub=diff_models(:,c);
            plot_avg=mean(plot_data_sub,1);
        end
        
        x = repmat(location_x(c),1,num_subjects);
        h(c,1) = scatter(x,plot_data_sub,[],color_subjects,'filled','SizeData',100,'DisplayName','Single subject'); %subject-level data
        hold on;
        h(c,2) = scatter(location_x(c),plot_avg,[],color_avg,'filled','SizeData',200,'DisplayName','Average'); %avg data
        if pl==1
            line([location_x(c)-0.05,location_x(c)+0.05],[RT_noise_ceiling_lower(c),RT_noise_ceiling_lower(c)],'Color',color_noise_ceiling); %noise ceiling
        end

    end

    %% Plot parameters
    % xlabel('Corr human-RNN RTs (Pearson''s r)');
    % ylabel('Conditions');
    
    xlim([0.2 3.2]);
    ylim([-0.4 0.8])
    % legend(h(1,:),'Location','southeast');
    % set(gca,'xtick',location_x,'xticklabel',xticklabels);
    set(gca,'xtick',label_x, 'xticklabel','');
    font_size = 18;
    set(gca,'FontName','Arial','FontSize',font_size);

    %% Save
    %RT plots
    if pl==1
        plot_name=sprintf('corr_RT_human_%s_line_ceil_readouts',model_type);
    elseif pl==2
        plot_name=sprintf('corr_RT_human_%s_diff_readouts',model_type);
    end
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s.fig',plot_name)));
    saveas(gcf,fullfile(results_avg_dir,sprintf('%s.svg',plot_name)));
    close(gcf);

    %Diff waves
    if pl==2 && ~strcmp(model_type,'bl')
        save(fullfile(results_avg_dir,sprintf('diff_wave_corr_RT_human_%s_allsubs_readouts',model_type)),'diff_models');
    end
    %Noise ceilings
    % save(fullfile(results_avg_dir,'corr_RT_RNN_human_noise_ceiling_lower'),'RT_noise_ceiling_lower');
end

end