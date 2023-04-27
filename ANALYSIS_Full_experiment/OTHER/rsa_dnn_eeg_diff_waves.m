function rsa_dnn_eeg_diff_waves(model_type,with_stats,with_error_bars)
%RSA_DNN_EEG_DIFF_WAVES Plot the difference waves (BLNet-model) for the RSA
%analysis.
%
%Input: model_type (model compared to BLNet: 'b' or 'b_d), with_stats (1/0)
%
%Output: 1 plot per scene condition (all/natural/manmade) containing 3
%difference waves (1 per layer), with/without stats.
%
%% Load the RSA results for both models
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
plot_location = -0.13:-0.015:-0.2; %where stats will be
for c=1:3
    if c == 1
        conditions = 'artificial';
    elseif c == 2
        conditions = 'natural';
    elseif c == 3
        conditions = 'all';
    end
    %feedforward model
    load(fullfile(results_avg_dir,sprintf('all_subjects_all_layers_rsa_%s_net_%s',model_type,conditions)),'rsa_all_subs');
    rsa_ff_model_all_subs=rsa_all_subs; 
   
    %bl-net
    load(fullfile(results_avg_dir,sprintf('all_subjects_all_tps_rsa_rnn_%s_10k_perms',conditions)),'rsa_all_subs')
    rsa_bl_net_all_subs=squeeze(median(rsa_all_subs,3)); %median over timesteps
    
    %Define some variables
    numSubjects=size(rsa_ff_model_all_subs,1);
    numLayers=size(rsa_ff_model_all_subs,2);
    numTimepointsEEG=size(rsa_ff_model_all_subs,3);
    
    %Plot diff wave
    diff_waves=rsa_bl_net_all_subs-rsa_ff_model_all_subs;
    diff_waves_avg=squeeze(mean(diff_waves,1));
    
    figure;
     
    for layer=1:numLayers
        if layer==1
            color=[0.8 0.8 0.8];
        elseif layer==2
            color=[0.5 0.5 0.5];
        elseif layer==3
            color=[0 0 0];
        end
        diff_l=squeeze(diff_waves(:,layer,:));
        diff_l_avg=squeeze(diff_waves_avg(layer,:));
        peak_value=max(diff_l_avg);
        peak_latency=find(diff_l_avg==max(diff_l_avg));
        fprintf('The peak value is %f at %d ms \n',peak_value,(peak_latency-40)*5);
        %Perform and plot stats if needed
        if ~with_stats
            plot(diff_l_avg,'LineWidth',2,'Color',color);
            hold on;    
        
        %Stats: FDR-corrected
        else
            %Check if exist
            filename = fullfile(results_avg_dir,sprintf('%s_net_diff_wave_stats_%s_layer_%d.mat',...
                model_type,conditions,layer));
            if exist(filename,'file')
                load(filename,'fdr_stats');
            else %if not, run them
                fdr_stats.num_perms = 10000;
                fdr_stats.tail = 'both';
                fdr_stats.qvalue = 0.05;
                [fdr_stats.significant_timepoints,fdr_stats.crit_p,fdr_stats.adjusted_pvalues] = ...
                    fdr_rsa_random_effects_stats(diff_l,fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);
                save(filename,'fdr_stats');
            end

            %Plot the stats
            st = (fdr_stats.significant_timepoints*plot_location(layer));
            st(st==0) = NaN;
            plot(st,'*','Color',color);
            hold on;

            %if needed: error bars
            if with_error_bars
                
                %calculate error bars
                stdDM = std(diff_l);
                err = stdDM/sqrt(numSubjects); %standard deviation/sqrt of num subjects

                %plot as a shaded area
                top_curve = diff_l_avg + err;
                bottom_curve = diff_l_avg - err;
                x2 = [1:numTimepointsEEG, fliplr(1:numTimepointsEEG)];
                shaded_area = [top_curve, fliplr(bottom_curve)];
                fill(x2, shaded_area, color,'FaceAlpha',0.2);
                hold on;
            end

            %Plot the data
            plot(diff_l_avg,'LineWidth',2,'Color',color);
            hold on;
        end
    end
    
    %Plot parameters
    ylim([-0.2 0.2]);
%     yticks(-0.2:0.2:0.8);
    font_size = 18;
    set(gca,'FontName','Arial','FontSize',font_size);
    xticks(0:20:200);
    xticklabels(-200:100:800);        
    onset_time = 40;
    xline(onset_time,'--');
    
    plot_whole_epoch=0;
    if ~plot_whole_epoch 
        xlim([20 200]);
    end

    filename=fullfile(results_avg_dir,sprintf('diff_wave_%s_net_vs_bl_net_%s',model_type,conditions));
    save(filename,'diff_waves');
    saveas(gcf,sprintf('%s.fig',filename));
    saveas(gcf,sprintf('%s.svg',filename));
end
close all;
end


