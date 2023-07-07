function dnn_dth_diff_wave(with_stats,with_error_bars,model_type) 
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
    stats_location=-0.22:-0.02:-0.26;
end

figure;

%Loop over conditions & calculate the diff wave
for c=1:3
    if c==1 %artificial
        color=[0.8 0.8 0.8];
        color_manm=color;
        conditions='artificial';
        var_name='correlation_art';
        
    elseif c==2 %natural
        color=[0.5 0.5 0.5];
        color_nat=color;
        conditions='natural';
        var_name='correlation_nat';
    elseif c==3 %both
        color=[0 0 0];
        color_both=color;
        conditions='both';
        var_name='correlation_both';
    end
   
    %load the results
    load(fullfile(results_avg_dir,sprintf('dth_results_%s_net_all_subjects_%s_readouts',model_type,conditions)),var_name); %feedforward
    results_ff=eval(var_name);
    load(fullfile(results_avg_dir,sprintf('dth_results_bl_net_all_subjects_%s',conditions)),var_name); %blnet
    results_blnet=eval(var_name);
    
    %calculate and plot the difference wave
%     diff_wave_all_subs_cond=(results_blnet(subjects,:)*-1)-(results_ff(subjects,:)*-1);
    diff_wave_all_subs_cond=(results_blnet(subjects,:))-(results_ff(subjects,:));
    diff_wave_all_subs(:,c,:)=diff_wave_all_subs_cond;
    diff_wave_avg=mean(diff_wave_all_subs_cond,1);
    
    if c==1
        diff_wave_manm=diff_wave_avg;
    elseif c==2
        diff_wave_nat=diff_wave_avg;
    elseif c==3
        diff_wave_both=diff_wave_avg;
    end
    
    %Stats if needed
    if with_stats  
        %Check if stats already exist, otherwise run the stats script
        filename = sprintf('dth_diff_wave_bl_net_minus_%s_net_readouts',model_type);       
        filename_stats = fullfile(results_avg_dir,sprintf('%s_%s_stats.mat',filename,conditions));
        if exist(filename_stats,'file')
            load(filename_stats,'fdr_stats');
        else
            fdr_stats.num_perms = 10000;
            fdr_stats.tail = 'both';
            fdr_stats.qvalue = 0.05;
            [fdr_stats.significant_timepoints,fdr_stats.pvalues,...
                fdr_stats.crit_p, fdr_stats.adjusted_pvalues]...
                = fdr_permutation_stats(diff_wave_all_subs_cond,...
                fdr_stats.num_perms,fdr_stats.tail,fdr_stats.qvalue);
            save(filename_stats,'fdr_stats');
        end
       
        % Plot significance
        st = (fdr_stats.significant_timepoints*stats_location(c)); 
        st(st==0) = NaN;
        plot(st,'*','Color',color); 
        hold on;
    end
    
    %error bars
    if with_error_bars
        %2) error bars
        stdDM = std(diff_wave_all_subs_cond); 
        err = stdDM/sqrt(numel(subjects)); %standard deviation/sqrt of num subjects  

        %plot as a shaded area
        top_curve = diff_wave_avg + err;
        bottom_curve = diff_wave_avg - err;
        x2 = [1:numTimepoints, fliplr(1:numTimepoints)];
        shaded_area = [top_curve, fliplr(bottom_curve)];
        fill(x2, shaded_area, color,'FaceAlpha',0.5);
        hold on;
    end
end 
plot(diff_wave_manm,'LineWidth',2,'Color',color_manm);
hold on;
plot(diff_wave_nat,'LineWidth',2,'Color',color_nat);
hold on;
plot(diff_wave_both,'LineWidth',2,'Color',color_both);
hold on;

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
save(fullfile(results_avg_dir,filename),'diff_wave_all_subs');
saveas(gcf,fullfile(results_avg_dir,sprintf('%s_readouts.fig',filename)));
saveas(gcf,fullfile(results_avg_dir,sprintf('%s_readouts.svg',filename)));
close(gcf);

end



