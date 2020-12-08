function Preprocessing_PDM_full_experiment(subn,task) 
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
mainDirectoryEEG = dir(fullfile(data_path,subname,'/EEG/','*.eeg')); 
fileNameEEG = fullfile(data_path, subname,'/EEG/', mainDirectoryEEG.name); 
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
%Behavioural data
behav.triggers = [];
behav.RT = [];
behav.points = [];

for block = 1:numel(filenamesSorted)
    block_name = filenamesSorted{block};
    if contains(block_name,task_name)
        load(block_name);
        behav.triggers = [behav.triggers;data.triggers];
        behav.RT = [behav.RT; data.rt]; 
        behav.points = [behav.points; data.points]; 
    end
end

%EEG data
%Define Events
cfg=[];
cfg.dataset=fileNameEEG;
cfg.trialdef.eventtype='Stimulus';
cfg.trialdef.prestim=0.2;
cfg.trialdef.poststim=0.8;
cfg=ft_definetrial(cfg);

%Define the columns of cfg
eegtriggers = cfg.trl(:,4);
beginningEpoch = cfg.trl(:,1); %beginning of each trial relative to the beginning of the raw data
endEpoch = cfg.trl(:,2); %end of each trial; column 2-column 1 = duration of epoch
offsetTrigger = cfg.trl(:,3); %offset of the trigger with respect to the trial - defined by cfg.trialdef.prestim

numConditions = 60;
if task == 1
    task_triggers = 1:numConditions; 
elseif task == 2
    task_triggers = (1:numConditions)+100;
end

triggers_eeg = [];
trials_remaining = [];
for t = 1:numel(eegtriggers)
    if ismember(eegtriggers(t),task_triggers)
        trials_remaining = [trials_remaining;t];
        triggers_eeg = [triggers_eeg;eegtriggers(t)];
    end
end
if task == 2
    triggers_eeg = triggers_eeg-100;
end
%Remove the same trials from the other columns of cfg
beginningEpoch = beginningEpoch(trials_remaining);
endEpoch = endEpoch(trials_remaining);
offsetTrigger = offsetTrigger(trials_remaining);

%% Remove extra triggers
%remove  paperclips from behavioural data
paperclips = behav.triggers==999;
behav.triggers = behav.triggers(~paperclips);
behav.RT = behav.RT(~paperclips);
behav.points = behav.points(~paperclips);

for i = 1:numel(behav.triggers)
    while round(triggers_eeg(i)) ~= round(behav.triggers(i))
        triggers_eeg(i) = [];
        beginningEpoch(i)=[];
        endEpoch(i)=[];
        offsetTrigger(i)=[];
    end
end

if isequal(round(triggers_eeg),round(behav.triggers)) 
    disp('all good!');
else
    warning('problem with the triggers: check manually');
    keyboard;
end

%put back into the configuration file
cfg.trl = [beginningEpoch endEpoch offsetTrigger triggers_eeg];

%% Rest of preprocessing

%Filter and read
cfg.dataset=fileNameEEG;
cfg.lpfilter = 'yes';
cfg.lpfreq = 50;
cfg.demean='yes'; %apply baseline correction
cfg.baselinewindow=[-0.2,0];
data=ft_preprocessing(cfg);

%Resample data to 200hz
cfg=[];
cfg.resamplefs=200;
data=ft_resampledata(cfg,data);

% Get rid of jumps
%The data are preprocessed (again) with the following configuration parameters,
%which are optimal for identifying jump artifacts.
cfg.artfctdef.jump.medianfilter  = 'yes';
cfg.artfctdef.jump.medianfiltord = 9;
cfg.artfctdef.jump.absdiff       = 'yes';

%  Artifacts are identified by means of thresholding the z-transformed value
%  of the preprocessed data.
cfg.artfctdef.jump.channel       = {'eeg'}; %Nx1 cell-array with selection of channels, see FT_CHANNELSEL
cfg.artfctdef.jump.cutoff        = 20; %z-value at which to threshold (default = 20)
[cfg, ~] = ft_artifact_jump(cfg, data);

%Remove the found artifacts
cfg.feedback = 'yes';
data = ft_rejectartifact(cfg, data);

% Display Summary Stats and exclude epochs/channels
cfg=[];
cfg.showlabel='yes';
cfg.method='summary';
%cfg.layout='easycapM1.lay'; %not the one we used
data=ft_rejectvisual(cfg,data);

%Update the behavioural triggers, RT and points and save them
for i = 1:numel(data.trialinfo)
    while round(behav.triggers(i)) ~= round(data.trialinfo(i))
        behav.triggers(i) = [];
        behav.RT(i) = [];
        behav.points(i) = [];
    end
end

%Double-check that eeg and behavioural triggers are equal
if isequal(round(data.trialinfo),round(behav.triggers)) 
    disp('all good!');
else
    warning('problem with the triggers: check manually');
    keyboard;
    behav.triggers(end) = [];
    behav.RT(end) = [];
    behav.points(end) = [];
end
task_name_small = get_task_name(task);
save(sprintf('%s/%s/preprocessed_behavioural_data_%s',data_path, subname, task_name_small),'behav');

%Transform to "timelocked" data and save the output
cfg.outputfile= sprintf('%s/%s/timelock_%s',data_path,subname,task_name_small); 
cfg.keeptrials='yes';
data=ft_timelockanalysis(cfg,data);
