function plot_topography_searchlight (varargin)
%Author: Siying Xie (2022)

% default
obj = 'ctg';
% popu = 'ift';
chan = 'whol23c';
% meas = 'lda';

% custom setting
% for i = 1:2:length(varargin)
%     param = varargin{i};
%     value = varargin{i+1};
%     if ~ischar(param)
%         error('Flag arguments must be strings')
%     end
%     param = lower(param);
%     switch param
%         case 'popu'
%             popu = value;
%             if ~strcmp(popu,'ift') && ~strcmp(popu,'adt')
%                 error('Value of ''popu'' must be ''ift'' or ''adt''.');
%             end
%         case 'obj'
%             obj = value;
%             if ~strcmp(obj,'exp') && ~strcmp(obj,'ctg')
%                 error('Value of ''obj'' must be ''exp'' or ''ctg''.');
%             end
%         case 'meas'
%             meas = value;
%         case 'chan'
%             chan = value;
%     end
% end

% define path 
%EEGLAB toolbox
toolboxDir = '/home/siyingxie/gitLib';
%Data
scratchDir = '/scratch/siyingxie/projects/IFTADT';

addpath(genpath(fullfile(toolboxDir, 'HelperFunctions')));
addpath(genpath(fullfile(toolboxDir, 'eeglab')));

% define measurement
% if strcmp(meas, 'corr')
%     measureName = '@cosmo_correlation_measure';
% end
% if strcmp(meas, 'lda')
%     measureName = '@cosmo_crossvalidation_measure';
% end

%% Load Results
% load channel location
load( fullfile(scratchDir, 'data_info/adult_64channels.mat'),'c64'); channelloc = c64;

%load results
sl_result_dir = fullfile(scratchDir, 'Data/Result_searchlight', [popu, '_', obj, measureName]);
files = dir(fullfile(sl_result_dir, '*.mat'));
for subn = 1:length(files)
    load(fullfile(sl_result_dir,files(subn).name),'sl_tl_ft');
    topo.([popu, '_', obj,'_', meas, '_measure'])(subn,:,:) = sl_tl_ft.avg;
end

% Figure
%% instead of this, select the data from 11 timepoints: -195 ms, -95 ms, 5 ms,...,905 ms
% times = -200:2:998;
% time_bins = [ -100, 0, 200, 400, 600, 800];
% for tbi = 1:(length(time_bins)-1) % n-1 bins
%     bin = [time_bins(tbi) time_bins(tbi+1)]';
%     [tp2use] = dsearchn(times', bin);
%     res2plot = topo.([popu, '_', obj,'_', meas, '_measure'])(:,:,tp2use(1):tp2use(2));
%     searchlight_patterns(tbi,:)= squeeze(nanmean(nanmean(res2plot,1),3));
% end
% bin = [time_bins(length(time_bins)) 998]'; % the last bin
% [tp2use] = dsearchn(times', bin);
% res2plot = topo.([popu, '_', obj,'_', meas, '_measure'])(:,:,tp2use(1):tp2use(2));
% searchlight_patterns(length(time_bins),:) = squeeze(nanmean(nanmean(res2plot,1),3));

% if strcmp('lda', meas)
%     searchlight_patterns = (searchlight_patterns-0.25)*100;
% end

%% prepare color map - could be just 40-100 ms
% set anchor
if strcmp('corr', meas)
    MID=round(abs(max(max(searchlight_patterns))),2);
    mid=round(abs(min(min(searchlight_patterns))),2);
    MID_portion = round(MID/(MID+mid),2);
%     mid_portion = round(mid/(MID+mid),2);
    clim = [-mid MID];
end
if strcmp('lda', meas)
    mid = 0;
    MID = round(abs(max(max(searchlight_patterns))),2);
    MID_portion = 1;
%     mid_portion = 0.5;
    clim=[mid, MID];
end
warning('off');
color_upper = cbrewer('seq', 'YlOrRd', 100*MID_portion); %maybe rather blue-green or magenta, to fit the decoding plot
% color_lower = flipud(cbrewer('seq', 'Blues', 100*mid_portion));
% colors = cat(1, color_lower, color_upper);
colors =  color_upper;

%% plt the 11 plots
figure;
for tbi = 1:length(time_bins) %1:11 
    subplot(2,4,tbi);
    topoplot(squeeze(searchlight_patterns(tbi,:)), channelloc, 'colormap', colors, 'style' ,'both');
    caxis(clim);
    
    titleBins = {'-100to0 ms', '0to200 ms','200to400 ms','400to600 ms','600to800 ms', '800to1000ms'};
    title(titleBins{1,tbi});
end

subplot(2,4,7);
CBar_Handle = colorbar('West');
caxis(clim);

set(get(CBar_Handle, 'YLabel'), ...
    'String', ['Decoding accuracy' newline ...
    'minus chance (%)'], 'FontSize', 10, 'FontName', 'Arial');


set(gca,'visible', 'off');
set(gcf, 'color','white');

fig_width = 600;
fig_height = 300;
set(gcf,'position',[100 100 fig_width fig_height]);
set(gcf,'Renderer','painters');

%% save
% save figure
figureSavePath =  fullfile('~/Documents/Projects/Infant/Figure_pro_v2.0','topog', [num2str(fig_width),'x', num2str(fig_height),'_figSize']);
if ~exist(figureSavePath,'dir');mkdir(figureSavePath);end
figureName =  ['sl_', popu,'_',obj, '_', chan, '_', meas];
print (gcf, fullfile(figureSavePath, figureName), '-dpng', '-r300');
% print (gcf, fullfile(figureSavePath, figureName), '-dpsc', '-r300');
end