function plotting_parameters(plot_title,title_bool,plot_legend,legend_bool,font_size,legend_location,label_y)
%PLOTTING_PARAMETERS Specify plot parameters.
%
%Input: title of the plot, title bool (with or without title), legend of the plot, stimulus onset, font size,
%legend location, label of the y axis
%
set(gca,'FontName','Arial')
xticks(0:20:200);
xticklabels(-200:100:800);
xlim([0 201]);
xlabel('Time (ms)');
ylabel(label_y);   
if title_bool
    title(plot_title);
end
onset_time = 40;
xline(onset_time,'--');
if legend_bool
    legend(plot_legend,'Location',legend_location,'FontSize',font_size);
end