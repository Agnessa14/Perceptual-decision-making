function reconstructing_data(subn,~)
% RECONSTRUCTING_DATA For subjects 2,3,4 that are missing data, reconstructing the vmrk file with the correct triggers.
%
%Input: subject ID (e.g., 1), task (1=categorization, 2=distraction)
%
%
%

%% Add paths and rewrite subject & task name
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
subname = get_subject_name(subn);
% task_name = get_task_name_capitalized(task);
data_path = '/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT';
addpath(genpath(fullfile(data_path,subname)));
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

%%  Load  behavioural data
mainDirectoryBeh = dir(fullfile(data_path,subname,'/Behavioural/','*.mat')); 

%% Fix order of blocks since they do not load in the correct order
filenamesAll = extractfield(mainDirectoryBeh,'name');
for m = 1:numel(mainDirectoryBeh)  
    parsedName = split(filenamesAll{m},'_');
    parsedBlocknum = split(parsedName{4},'.');
    mainDirectoryBeh(m).blockNum = str2double(parsedBlocknum{1});     
end

%reorder by block number
T = struct2table(mainDirectoryBeh); 
sortedT = sortrows(T, 'blockNum'); 
mainDirectoryBeh = table2struct(sortedT); 
filenamesSorted = extractfield(mainDirectoryBeh,'name');

%double check that the order is correct
blocks = 1:numel(mainDirectoryBeh);
blockorderarray = arrayfun(@(x) mainDirectoryBeh(x).blockNum,blocks);
if ~isequal(blockorderarray,blocks)
    error('Incorrect block order: Check and fix manually');
end

%% Take the task-related trials only
behav.triggers = [];
behav.RT = [];
behav.points = [];
behav.reconstructed_data = [];

for block = 1:numel(filenamesSorted)
    block_name = filenamesSorted{block};
%     if contains(block_name,task_name)
    t = 1;
    load(block_name);
    behav.triggers = [behav.triggers;data.triggers];
    behav.RT = [behav.RT; data.rt]; 
    behav.points = [behav.points; data.points]; 

    %get the timing+triggers from  behavioural data
    timing_minus_onset = [data.tStimStart-data.tOnset data.timePressed-data.tOnset];
    timing_minus_onset(1,1) = 0;
    reshaped_rowmaj = reshape(timing_minus_onset.',[],1);
    with_triggers = [reshaped_rowmaj NaN(size(reshaped_rowmaj))];
    for s=1:size(reshaped_rowmaj,1)
        if mod(s,2) %odd
            with_triggers(s,2) = data.triggers(t);
            t=t+1;
        else
            with_triggers(s,2) = 222;
        end
    end

    nan_trials = isnan(with_triggers(:,1));
    with_triggers_nonans = with_triggers;
    with_triggers_nonans(nan_trials,:) = [];
    reconstructed_data_block = [0 200; with_triggers_nonans; with_triggers_nonans(end,1)+1 244];
    behav.reconstructed_data = [behav.reconstructed_data; reconstructed_data_block];
%     end
end

%% Create the cell array that will be turned into a vmrk file
num_rows = size(behav.reconstructed_data,1)+1;
all_rows = 1:size(behav.reconstructed_data,1);
rows=2:num_rows;
vmrk_cell = cell(num_rows,5);

%Fill the first row manually
vmrk_cell{1,1} = 'Mk=New Segment,';
vmrk_cell{1,2} = 1;
vmrk_cell{1,3} = 1;
vmrk_cell{1,4} = 0;
vmrk_cell{1,5} = 20201126171218962590;

%Column 1
vmrk_cell(rows,1) = arrayfun(@(x) sprintf( 'M%d',x),all_rows, 'UniformOutput',false);

%Column 2
if subn==2
    onset_eeg = 202565; %taken manually from the vmrk file
end

vmrk_cell(rows,2) = arrayfun(@(x) behav.reconstructed_data(x,1)+onset_eeg,all_rows','UniformOutput',false);

%Column 3: Triggers
vmrk_cell(rows,3) = arrayfun(@(x) sprintf( 'S%d',behav.reconstructed_data(x,2) ),all_rows','UniformOutput',false);

%Columns 4 and 5
for row=rows
    vmrk_cell{row,4} = 1;
    vmrk_cell{row,5} = 0;
end

keyboard;


