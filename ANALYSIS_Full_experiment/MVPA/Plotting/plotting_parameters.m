function plotting_parameters(plot_title,plot_legend,onset_time,font_size,legend_location)
%PLOT_DATA Plot input data.
%
%Input: title of the plot, legend of the plot, stimulus onset (in seconds)
%
xlim([0 201]);
xlabel('Timepoints');
ylabel('Decoding accuracy');   
title(plot_title);
xline(onset_time,'--');
legend(plot_legend,'Location',legend_location,'FontSize',font_size);