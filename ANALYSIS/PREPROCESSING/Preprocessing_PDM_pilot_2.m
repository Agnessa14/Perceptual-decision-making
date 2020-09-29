function Preprocessing_PDM_pilot_2(subn) 
% PREPROCESSING_PDM_pilot_2 Preprocess the EEG data for the PDM experiment for one
% participant.
%
%Steps include:
%-removing extra EEG data
%-removing the paperclip trials and response triggers
%-filtering
%-downsampling
%-getting rid of jump artifacts
%-transform into timelock.m

%% Add paths
if subn < 10
    subnum = ['0',num2str(subn)];
else
    subnum = num2str(subn);
end

addpath(genpath(fullfile('/scratch/agnek95/PDM/DATA/DATA_PILOT_2/',subnum)));
addpath('/home/agnek95/OR/TOOLBOX/fieldtrip-20190224');
ft_defaults;

%%  Specify input file and read events
mainDirectoryEEG = dir(fullfile('/scratch/agnek95/PDM/DATA/DATA_PILOT_2/',subnum,'/EEG/','*.eeg')); 
fileName = fullfile('/scratch/agnek95/PDM/DATA/DATA_PILOT_2/', subnum,'/EEG/', mainDirectoryEEG.name); 

%Define Events
cfg=[];
cfg.dataset=fileName;
cfg.trialdef.eventtype='Stimulus';
cfg.trialdef.prestim=0.2;
cfg.trialdef.poststim=0.8;
cfg=ft_definetrial(cfg);

%load behavioural data
mainDirectoryBeh = dir(fullfile('/scratch/agnek95/PDM/DATA/DATA/',subnum,'/Behavioural/','*.mat')); 

%% Fix order of blocks since they do not load in the correct order
for m = 1:numel(mainDirectoryBeh)
    filename = extractfield(mainDirectoryBeh(m),'name');
    filename = filename{1}; %in string format
    parsedName = split(filename,'_');
    mainDirectoryBeh(m).blockNum = str2double(parsedName{4});     
end

T = struct2table(mainDirectoryBeh); 
sortedT = sortrows(T, 'blockNum'); 
mainDirectoryBeh = table2struct(sortedT); 

%double-check that the order is correct: the problem would occur at block 9
% if mainDirectoryBeh(9).blockNum ~= 9
%     error('Incorrect block order: Check and fix manually')
% end

%% Remove unnecessary triggers
%Create list of triggers from the behavioural data
filenamesArray = extractfield(mainDirectoryBeh, 'name')';
triggers = [];
for n = 1:numel(filenamesArray)
    load(char(filenamesArray{n}));
    triggers = [triggers;blockData.samples];
end

%Compare eeg triggers with behavioural triggers
eegtriggers = cfg.trl(:,4);
beginningEpoch = cfg.trl(:,1); %beginning of each trial relative to the beginning of the raw data
endEpoch = cfg.trl(:,2); %end of each trial; column 2-column 1 = duration of epoch
offsetTrigger = cfg.trl(:,3); %offset of the trigger with respect to the trial - defined by cfg.trialdef.prestim

%Remove beginning of block (200), end of block (244) and response (200) triggers
removed = find(eegtriggers==200|eegtriggers==222|eegtriggers==244);
eegtriggers(eegtriggers==200|eegtriggers==222|eegtriggers==244)=[];
beginningEpoch(removed)=[];
endEpoch(removed)=[];
offsetTrigger(removed)=[];

%Change 99 to 999 to be able to compare between eegtriggers and beh. triggers
eegtriggers(eegtriggers==99) = 999;

%Remove extra EEG triggers
for i = 1:numel(triggers)
    while round(triggers(i)) ~= round(eegtriggers(i))
        eegtriggers(i) = [];
        beginningEpoch(i)=[];
        endEpoch(i)=[];
        offsetTrigger(i)=[];
    end
end

%Double-check that eeg and behavioural triggers are equal
if isequal(round(eegtriggers),round(triggers)) 
    disp('all good!');
else
    warning('problem with the triggers: check manually');
    keyboard;
end

%Remove the paperclip trials
paperclipEeg = find(eegtriggers==999);
eegtriggers(eegtriggers==999)=[];
beginningEpoch(paperclipEeg)=[];
endEpoch(paperclipEeg)=[];
offsetTrigger(paperclipEeg)=[];

%put back into the configuration file
cfg.trl = [beginningEpoch endEpoch offsetTrigger eegtriggers];
%maybe also remove those trials from cfg.event

%% Rest of preprocessing

%Filter and read
cfg.dataset=fileName;
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
cfg.artfctdef.jump.cutoff        = 20; %volts - z-value at which to threshold (default = 20)
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

%Transform to "timelocked" data and save the output
cfg.outputfile= sprintf('/scratch/agnek95/PDM/DATA/DATA_PILOT_2/%s/timelock',subnum); 
cfg.keeptrials='yes';
data=ft_timelockanalysis(cfg,data);
