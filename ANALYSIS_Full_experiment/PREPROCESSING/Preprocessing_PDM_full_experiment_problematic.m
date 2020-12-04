function Preprocessing_PDM_full_experiment_problematic(subn,task) 
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

%The triggers for 256-260 have been converted to 0-4 and need to be
%converted back
for t = 1:numel(eegtriggers)
    if eegtriggers(t) <= 4
        eegtriggers(t) = eegtriggers(t)+256;
    end
end

%Problem: response trigger is 222, which is also the fixation trigger for
%stimulus 22. End of block trigger is 244, also the fixation trigger for
%44, but that one is easier to remo ve. Need to remove these triggers without removing the experimental trials.

%Remove the other task blocks
conditions = 1:60; 
task_triggers = conditions+task*100;
block_starts = find(eegtriggers==200);
numBlocks = numel(block_starts);
block_triggers = cell(numel(block_starts),1);

for bs = 1:numBlocks
    if bs < numBlocks
        block_triggers{bs} = eegtriggers(block_starts(bs):block_starts(bs+1)-1);
    else
        block_triggers{bs} = eegtriggers(block_starts(bs):end);
    end
end

for b = 1:numBlocks
    t = 2;
    if ~ismember(block_triggers{b}(t),task_triggers)
        block_triggers{b}(:) = NaN;
    elseif block_triggers{b}(t) == 199
        t=t+2;
        if ~ismember(block_triggers{b}(t),task_triggers)
            block_triggers{b}(:) = NaN;
        end
    end
end    
       
%Pick only the task-relevant triggers
nan_or_not = cell2mat(cellfun(@isnan, block_triggers, 'UniformOutput', false));
trials_remaining = find(nan_or_not==0);
triggers_eeg = eegtriggers(trials_remaining);
beginningEpoch = beginningEpoch(trials_remaining);
endEpoch = endEpoch(trials_remaining);
offsetTrigger = offsetTrigger(trials_remaining);

%Remove the 999, 200 and 244 trials if they are followed by 200
for t = 1:numel(triggers_eeg)
    if triggers_eeg(t) == 244
        if (t<numel(triggers_eeg) && triggers_eeg(t+1) == 200) 
            triggers_eeg(t:t+1) = 555;
            disp('replaced end of block trial');
        elseif t == numel(triggers_eeg) 
            triggers_eeg(t) = 555;
            disp('replaced end of block trial');           
        end
    elseif triggers_eeg(t) == 200 || triggers_eeg(t) == 43 || triggers_eeg(t) == 222
        triggers_eeg(t) = 555;
        disp('replaced a block start, a paperclip or a keystroke');
    end
end

beginningEpoch(triggers_eeg==555) = [];
endEpoch(triggers_eeg==555) = [];
offsetTrigger(triggers_eeg==555) = [];
trials_remaining(triggers_eeg==555)= [];
triggers_eeg(triggers_eeg==555) = [];

%Translate back into amodal triggers
triggers_eeg = triggers_eeg-task*100;

%% Remove unnecessary triggers
%remove  paperclips from behavioural data
paperclips = behav.triggers==999;
behav.triggers = behav.triggers(~paperclips);
full_behav_triggers = behav.triggers;
behav.RT = behav.RT(~paperclips);
behav.points = behav.points(~paperclips);

%remove 56 & 22 temporarily 
trials_56 = find(behav.triggers==56);
trials_22 = find(behav.triggers==22);
behav.triggers([trials_56;trials_22]) = [];

%remove any extra eeg triggers (from repeating a block for example)
for i = 1:numel(behav.triggers)
    while round(behav.triggers(i)) ~= round(triggers_eeg(i))
        triggers_eeg(i) = [];
        beginningEpoch(i)=[];
        endEpoch(i)=[];
        offsetTrigger(i)=[];
        trials_remaining(i) = [];
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
end
task_name_small = get_task_name(task);
save(sprintf('%s/%s/preprocessed_behavioural_data_%s',data_path, subname, task_name_small),'behav');

%Transform to "timelocked" data and save the output
cfg.outputfile= sprintf('%s/%s/timelock_%s',data_path,subname,task_name_small); 
cfg.keeptrials='yes';
data=ft_timelockanalysis(cfg,data);
