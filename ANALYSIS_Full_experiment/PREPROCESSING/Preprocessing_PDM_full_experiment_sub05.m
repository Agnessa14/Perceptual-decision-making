function Preprocessing_PDM_full_experiment_sub05(subn,task) 
% PREPROCESSING_PDM_full_experiment Preprocess the EEG data for the PDM experiment (post-pilot data) for one
% participant and one task.
%
%Input: subject ID (e.g., 1), task (1=categorization, 2=distraction)
%
%Steps include:
%-removing extra EEG data
%-removing the paperclip trials and response triggers
%-filtering 
%-downsampling
%-getting rid of jump artifacts
%-transform into timelock.m
%
%Returns a matlab structure timelock.m containing the list of triggers and
%preprocessed data
%

%% Add paths and rewrite subject & task name
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment'));
subname = get_subject_name(subn);
task_name = get_task_name_capitalized(task);
data_path = '/scratch/agnek95/PDM/DATA/DATA_FULL_EXPERIMENT';
addpath(genpath(fullfile(data_path,subname)));
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

%%  Load EEG and behavioural data
fileNameEEG = fullfile(data_path, subname,'/EEG/', 'pdm_full_experiment_sub05.eeg'); 
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

%% Split into two groups
filenamesSorted_1_17 = filenamesSorted(1:17);
filenamesSorted_18_20 = filenamesSorted(18:20);

%% Take the task-related trials only: blocks 1-17
%Behavioural data
behav_1_17.triggers = [];
behav_1_17.RT = [];
behav_1_17.points = [];

for block = 1:numel(filenamesSorted_1_17)
    block_name = filenamesSorted_1_17{block};
    if contains(block_name,task_name)
        load(block_name);
        behav_1_17.triggers = [behav_1_17.triggers;data.triggers];
        behav_1_17.RT = [behav_1_17.RT; data.rt]; 
        behav_1_17.points = [behav_1_17.points; data.points]; 
    end
end

%EEG data
%Define Events
cfg_1_17=[];
cfg_1_17.dataset=fileNameEEG;
cfg_1_17.trialdef.eventtype='Stimulus';
cfg_1_17.trialdef.prestim=0.2;
cfg_1_17.trialdef.poststim=0.8;
cfg_1_17=ft_definetrial(cfg_1_17);

%Define the columns of cfg
eegtriggers_1_17 = cfg_1_17.trl(:,4);
beginningEpoch_1_17 = cfg_1_17.trl(:,1); %beginning of each trial relative to the beginning of the raw data
endEpoch_1_17 = cfg_1_17.trl(:,2); %end of each trial; column 2-column 1 = duration of epoch
offsetTrigger_1_17 = cfg_1_17.trl(:,3); %offset of the trigger with respect to the trial - defined by cfg.trialdef.prestim

numConditions = 60;
if task == 1
    task_triggers_1_17 = (1:numConditions); 
elseif task == 2
    task_triggers_1_17 = (1:numConditions)+100;
end

triggers_eeg_1_17 = [];
trials_remaining_1_17 = [];
for t = 1:numel(eegtriggers_1_17)
    if ismember(eegtriggers_1_17(t),task_triggers_1_17)
        trials_remaining_1_17 = [trials_remaining_1_17;t];
        triggers_eeg_1_17 = [triggers_eeg_1_17;eegtriggers_1_17(t)];
    end
end

if task == 2
    triggers_eeg_1_17 = triggers_eeg_1_17-100;
end

%Remove the same trials from the other columns of cfg
beginningEpoch_1_17 = beginningEpoch_1_17(trials_remaining_1_17);
endEpoch_1_17 = endEpoch_1_17(trials_remaining_1_17);
offsetTrigger_1_17 = offsetTrigger_1_17(trials_remaining_1_17);

%% Remove extra triggers
%remove  paperclips from behavioural data
paperclips = behav_1_17.triggers==999;
behav_1_17.triggers = behav_1_17.triggers(~paperclips);
behav_1_17.RT = behav_1_17.RT(~paperclips);
behav_1_17.points = behav_1_17.points(~paperclips);

for i = 1:numel(behav_1_17.triggers)
    while round(triggers_eeg_1_17(i)) ~= round(behav_1_17.triggers(i))
        triggers_eeg_1_17(i) = [];
        beginningEpoch_1_17(i)=[];
        endEpoch_1_17(i)=[];
        offsetTrigger_1_17(i)=[];
    end
end

if isequal(round(triggers_eeg_1_17),round(behav_1_17.triggers)) 
    disp('all good!');
else
    warning('problem with the triggers: check manually');
    keyboard;
end

%put back into the configuration file
cfg_1_17.trl = [beginningEpoch_1_17 endEpoch_1_17 offsetTrigger_1_17 triggers_eeg_1_17];

%% Take the task-related trials only: blocks 1-17
%Behavioural data
behav_18_20.triggers = [];
behav_18_20.RT = [];
behav_18_20.points = [];

for block = 1:numel(filenamesSorted_18_20)
    block_name = filenamesSorted_18_20{block};
    if contains(block_name,task_name)
        load(block_name);
        behav_18_20.triggers = [behav_18_20.triggers;data.triggers];
        behav_18_20.RT = [behav_18_20.RT; data.rt]; 
        behav_18_20.points = [behav_18_20.points; data.points]; 
    end
end

%EEG data
fileNameEEG_18_20 = fullfile(data_path, subname,'/EEG/', 'pdm_full_experiment_sub05_runs_18_20.eeg');

%Define Events
cfg_18_20=[];
cfg_18_20.dataset=fileNameEEG_18_20;
cfg_18_20.trialdef.eventtype='Stimulus';
cfg_18_20.trialdef.prestim=0.2;
cfg_18_20.trialdef.poststim=0.8;
cfg_18_20=ft_definetrial(cfg_18_20);

%Define the columns of cfg
eegtriggers_18_20 = cfg_18_20.trl(:,4);
beginningEpoch_18_20 = cfg_18_20.trl(:,1); %beginning of each trial relative to the beginning of the raw data
endEpoch_18_20 = cfg_18_20.trl(:,2); %end of each trial; column 2-column 1 = duration of epoch
offsetTrigger_18_20 = cfg_18_20.trl(:,3); %offset of the trigger with respect to the trial - defined by cfg.trialdef.prestim

numConditions = 60;
if task == 1
    task_triggers_18_20 = (1:numConditions); 
elseif task == 2
    task_triggers_18_20 = (1:numConditions)+100;
end

triggers_eeg_18_20 = [];
trials_remaining_18_20 = [];
for t = 1:numel(eegtriggers_18_20)
    if ismember(eegtriggers_18_20(t),task_triggers_18_20)
        trials_remaining_18_20 = [trials_remaining_18_20;t];
        triggers_eeg_18_20 = [triggers_eeg_18_20;eegtriggers_18_20(t)];
    end
end

if task == 2
    triggers_eeg_18_20 = triggers_eeg_18_20-100;
end

%Remove the same trials from the other columns of cfg
beginningEpoch_18_20 = beginningEpoch_18_20(trials_remaining_18_20);
endEpoch_18_20 = endEpoch_18_20(trials_remaining_18_20);
offsetTrigger_18_20 = offsetTrigger_18_20(trials_remaining_18_20);

%% Remove extra triggers
%remove  paperclips from behavioural data
paperclips = behav_18_20.triggers==999;
behav_18_20.triggers = behav_18_20.triggers(~paperclips);
behav_18_20.RT = behav_18_20.RT(~paperclips);
behav_18_20.points = behav_18_20.points(~paperclips);

for i = 1:numel(behav_18_20.triggers)
    while round(triggers_eeg_18_20(i)) ~= round(behav_18_20.triggers(i))
        triggers_eeg_18_20(i) = [];
        beginningEpoch_18_20(i)=[];
        endEpoch_18_20(i)=[];
        offsetTrigger_18_20(i)=[];
    end
end

if isequal(round(triggers_eeg_18_20),round(behav_18_20.triggers)) 
    disp('all good!');
else
    warning('problem with the triggers: check manually');
    keyboard;
end

%put back into the configuration file
cfg_18_20.trl = [beginningEpoch_18_20 endEpoch_18_20 offsetTrigger_18_20 triggers_eeg_18_20];

%% Rest of preprocessing
%Filter and read
cfg_1_17.dataset=fileNameEEG;
cfg_1_17.lpfilter = 'yes';
cfg_1_17.lpfreq = 50;
cfg_1_17.demean='yes'; %apply baseline correction
cfg_1_17.baselinewindow=[-0.2,0];
data_1_17=ft_preprocessing(cfg_1_17);

%Filter and read
cfg_18_20.dataset=fileNameEEG;
cfg_18_20.lpfilter = 'yes';
cfg_18_20.lpfreq = 50;
cfg_18_20.demean='yes'; %apply baseline correction
cfg_18_20.baselinewindow=[-0.2,0];
data_18_20=ft_preprocessing(cfg_18_20);

%Resample data to 200hz
cfg_1_17=[];
cfg_1_17.resamplefs=200;
data_1_17=ft_resampledata(cfg_1_17,data_1_17);

%Resample data to 200hz
cfg_18_20=[];
cfg_18_20.resamplefs=200;
data_18_20=ft_resampledata(cfg_18_20,data_18_20);

% Get rid of jumps
%The data are preprocessed (again) with the following configuration parameters,
%which are optimal for identifying jump artifacts.
cfg_1_17.artfctdef.jump.medianfilter  = 'yes';
cfg_1_17.artfctdef.jump.medianfiltord = 9;
cfg_1_17.artfctdef.jump.absdiff       = 'yes';

%  Artifacts are identified by means of thresholding the z-transformed value
%  of the preprocessed data.
cfg_1_17.artfctdef.jump.channel       = {'eeg'}; %Nx1 cell-array with selection of channels, see FT_CHANNELSEL
cfg_1_17.artfctdef.jump.cutoff        = 20; %z-value at which to threshold (default = 20)
[cfg_1_17, ~] = ft_artifact_jump(cfg_1_17, data_1_17);

%Remove the found artifacts
cfg_1_17.feedback = 'yes';
data_1_17 = ft_rejectartifact(cfg_1_17, data_1_17);

% Get rid of jumps
%The data are preprocessed (again) with the following configuration parameters,
%which are optimal for identifying jump artifacts.
cfg_18_20.artfctdef.jump.medianfilter  = 'yes';
cfg_18_20.artfctdef.jump.medianfiltord = 9;
cfg_18_20.artfctdef.jump.absdiff       = 'yes';

%  Artifacts are identified by means of thresholding the z-transformed value
%  of the preprocessed data.
cfg_18_20.artfctdef.jump.channel       = {'eeg'}; %Nx1 cell-array with selection of channels, see FT_CHANNELSEL
cfg_18_20.artfctdef.jump.cutoff        = 20; %z-value at which to threshold (default = 20)
[cfg_18_20, ~] = ft_artifact_jump(cfg_18_20, data_18_20);

%Remove the found artifacts
cfg_18_20.feedback = 'yes';
data_18_20 = ft_rejectartifact(cfg_18_20, data_18_20);

% Display Summary Stats and exclude epochs/channels
cfg_1_17=[];
cfg_1_17.showlabel='yes';
cfg_1_17.method='summary';
cfg_1_17.layout='easycapM1.lay'; %not the one we used
data_1_17=ft_rejectvisual(cfg_1_17,data_1_17);

% Display Summary Stats and exclude epochs/channels
cfg_18_20=[];
cfg_18_20.showlabel='yes';
cfg_18_20.method='summary';
cfg_18_20.layout='easycapM1.lay'; %not the one we used
data_18_20=ft_rejectvisual(cfg_18_20,data_18_20);

%Update the behavioural triggers, RT and points and save them
for i = 1:numel(data_1_17.trialinfo)
    while round(behav_1_17.triggers(i)) ~= round(data_1_17.trialinfo(i))
        behav_1_17.triggers(i) = [];
        behav_1_17.RT(i) = [];
        behav_1_17.points(i) = [];
    end
end

%Double-check that eeg and behavioural triggers are equal
if isequal(round(data_1_17.trialinfo),round(behav_1_17.triggers)) 
    disp('all good!');
else
    warning('problem with the triggers: check manually');
    keyboard;
end
% 
% %in the categorization case, this wasnt covered by the loop above &
% resulted in an error
% behav_1_17.triggers(end) = [];
% behav_1_17.RT(end) = [];
% behav_1_17.points(end) = [];

for i = 1:numel(data_18_20.trialinfo)
    while round(behav_18_20.triggers(i)) ~= round(data_18_20.trialinfo(i))
        behav_18_20.triggers(i) = [];
        behav_18_20.RT(i) = [];
        behav_18_20.points(i) = [];
    end
end

%Double-check that eeg and behavioural triggers are equal
if isequal(round(data_18_20.trialinfo),round(behav_18_20.triggers)) 
    disp('all good!');
else
    warning('problem with the triggers: check manually');
    keyboard;
end

%% Append both datasets
cfg = [];
cfg.keepsampleinfo='no'; %yes if they come from different datafiles
data = ft_appenddata(cfg, data_1_17, data_18_20);

% Append behavioural data
behav = struct;
behav.triggers = [behav_1_17.triggers;behav_18_20.triggers];
behav.RT = [behav_1_17.RT;behav_18_20.RT];
behav.points = [behav_1_17.points;behav_18_20.points];

%Transform to "timelocked" data and save the output
task_name_small = get_task_name(task);
save(sprintf('%s/%s/preprocessed_behavioural_data_%s',data_path, subname, task_name_small),'behav');
cfg.outputfile= sprintf('%s/%s/timelock_%s',data_path,subname,task_name_small); 
cfg.keeptrials='yes';
data=ft_timelockanalysis(cfg,data);
