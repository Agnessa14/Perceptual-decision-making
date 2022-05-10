function [neighbourhood_map, all_missing_channel_ids, data_mod] = correct_neighbourhood_map(channels_labels_in_data, data)
% This function is used to correct the neighbourhood matrix for grouping
% the electrodes with nearest k-neighbours
%
%
% Input: channels_labels_in_data : The channels in the data (Array of
%                                  channel names with size Num channels x 1)
%        data : The EEG data (num_conditions x num_conditions x num_channels x num_timepoints)
%
% Output: neighbourhood_map : A matrix (double) with a dimension of number
%                             of channels x k (number of neighbours to be
%                             grouped)
%         all_missing_channel_ids : An array with the indexes of the
%                                   missing channels in the neighbourhoods
%                                   file containing all channels returns
%                                   empty if there are no missing channels
%         data_mod : The modified data matrix with the missing channel data
%                added as NaN at the missing channel ids.
% Author: Muthukumar Pandaram

% Load channel locations
load('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/adult_63channels.mat')
all_channel_labels = extractfield(channelloc,'labels')';

% Check for missing channels
channels_missing = setdiff(all_channel_labels, channels_labels_in_data);
% disp('The missing channels are: ', channels_missing);

% load Monika's neighbourhood (neighbourhoods, double)
load('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/EEG_neighbourhoods.mat')
neighbourhood_map = neighbourhoods;

%Find the number id for the missing channel labels and remove them from the
%neighbourhoods map

all_missing_channel_ids = double.empty;
data_mod = data;

for i = 1:length(channels_missing)
    all_missing_channel_ids(i,1) = find(strcmp(all_channel_labels,channels_missing(i)));
end
 
all_missing_channel_ids = sort(all_missing_channel_ids);

for i = 1:length(channels_missing)
    for j = 1:size(neighbourhood_map,1)
        ind=neighbourhood_map(j,:)== all_missing_channel_ids(i,1);
        neighbourhood_map(j, ind) = NaN;        
    end
    substitute_missing_data = NaN(size(data,1),size(data,2),1,size(data,4));
    data_mod = cat(3,data_mod(:,:,1:all_missing_channel_ids(i,1)-1,:),substitute_missing_data,data_mod(:,:,all_missing_channel_ids(i,1):end,:));
end

if ~isempty(channels_missing)
    neighbourhood_map(all_missing_channel_ids,:) = NaN;    
end

end