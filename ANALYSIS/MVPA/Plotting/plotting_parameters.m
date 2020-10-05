function plotting_parameters(plot_title,plot_legend,onset_time)
%PLOT_DATA Plot input data.
%
%Input: title of the plot, legend of the plot, stimulus onset (in seconds)
%
xlim([0 201]);
xlabel('Timepoints');
ylabel('Decoding accuracy');   
title(plot_title);
legend(plot_legend);
xline(onset_time,'--');