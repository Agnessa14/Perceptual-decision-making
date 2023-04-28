function dnn_dth_diff_wave(with_stats,model_type) 
%DNN_DTH_DIFF_WAVE Calculates the difference wave for the DTH analysis
%between BLnet and a feedforward model. 
%
%Input: add stats (1 with, 0 without), model_type ('b' or 'b_d')
%
%Results in one difference wave plot (3 curves - all, natural and man-made
%scenes)
%
%Author: Agnessa Karapetian, 2023

%% Define some variables
results_avg_dir = '/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG';
subjects=[5:26,28:35];
numTimepoints=200;
diff_wave_all_subs=NaN(numel(subjects),3,numTimepoints);

if with_stats
    stats_location=-0.15:-0.02:-0.19;
end

figure;

%Loop over conditions & calculate the diff wave
for c=1:3
    if c==1 %artificial
        color=[0.8 0.8 0.8];
        conditions='artificial';
        var_name='correlation_art';
        
    elseif c==2 %natural
        color=[0.5 0.5 0.5];
        conditions='natural';
        var_name='correlation_nat';
    elseif c==3 %both
        color=[0 0 0];
        conditions='both';
        var_name='correlation_both';
    end
   
    %load the results
    load(fullfile(results_avg_dir,sprintf('dth_results_%s_net_all_subjects_%s',model_type,conditions)),var_name); %feedforward
    results_ff=eval(var_name);
    load(fullfile(results_avg_dir,sprintf('dth_results_bl_net_all_subjects_%s',conditions)),var_name); %blnet
    results_blnet=eval(var_name);
    
    %calculate and plot the difference wave
    diff_wave_all_subs_cond=(results_blnet(subjects,:)*-1)-(results_ff(subjects,:)*-1);
    diff_wave_all_subs(:,c,:)=diff_wave_all_subs_cond;
    diff_wave_avg=mean(diff_wave_all_subs_cond,1);
    plot(diff_wave_avg,'LineWidth',2,'Color',color);
    hold on;
    
    %Stats if needed
    if with_stats
        fdr_stats.num_perms = 10000;
        fdr_stats.tail = 'right';
        fdr_stats.qvalue = 0.05;
  
        %Check if stats already exist, otherwise run the stats script
        filename = sprintf('dth_diff_wave_bl_net_minus_%s_net_%s',model_type,conditions);       
        filename_full = fullfile(results_avg_dir,sprintf('%s_stats',filename));
        if exist(filename_full,'file')
            load(filename_full,'fdr_stats');
        else
            [fdr_stats.significant_timepoints,fdr_stats.pvalues,...
                fdr_stats.crit_p, fdr_stats.adjusted_pvalues]...
                = fdr_permutation_stats(diff_wave_all_subs_cond,...
                fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);
            save(filename_full,'fdr_stats');
        end
       
        % Plot significance
        st = (fdr_stats.significant_timepoints*stats_location(c)); 
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;
    end
end 


%% Plotting parameters
font_size = 18;
set(gca,'FontName','Arial','FontSize',font_size);
ylim([-0.3 0.3]);
yticks(-0.3:0.1:0.3);
xticks(0:20:200);
xticklabels(-200:100:800);  
xline(40,'--');
xlabel('Time (ms)');
ylabel('Difference (Spearman''s p)');

%% Save difference wave and figures
save(filename_full,'diff_wave_all_subs');
saveas(gcf,sprintf('%s.fig',filename_full))
saveas(gcf,sprintf('%s.svg',filename_full))
close(gcf);

end



