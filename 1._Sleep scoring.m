clc; clear; close all;

%addpath (this + all sub-paths)
addpath(genpath('Z:\EEGlandArmandsScoringScript\EEGlabtools'),'-end')
addpath('Z:\Lab_resources\LabVIEW\functions')
%important, remove EEGLAB's octave functions
rmpath(genpath(fileparts(which('isoctave'))))

%PARAMETERS
%----------
%START SEARCH PATH (opens folder browser, starting in this path)
startPath = [''];
%SCORING FILENAME (*.set)
rFile = 'ScoredEEG.set';

channels = 1:3;

%OPTIONS (true or false)
option.reExportDat = false; %re-export EEGLAB data files (ask again)


%MAIN SCRIPT
%-----------
%INIT
fprintf('%s\n%s\n',mfilename,repmat('-',size(mfilename)))
%files
[~,rFile,rExt] = fileparts(rFile);
fileSET = [rFile,'.set'];
fileFDT = [rFile,'.fdt'];
Files = {fileSET,fileFDT};

%% SELECT DATA PATH
rPath = uigetdir(startPath,'Select Data Path');
if ~ischar(rPath)
    fprintf(2,'No Path Selected\n')
    return
end
if rPath(end)==filesep; rPath(end)=[]; end
while ~strcmpi(pwd,rPath) %care for delayed changing of path
    cd(rPath)
end
fprintf('Data Path: %s\n',rPath)


%% EXPORT EEGLAB DATA FILES
ind = cellfun(@(x)exist(fullfile(rPath,x),'file')==2,Files);
if ~all(ind) || option.reExportDat
    fprintf('\nEXPORT EEGLAB DATA\n')
    if all(ind)
        str = {'\fontsize{11}';...
            'EEGLAB data-files alread exist!';...
            sprintf('(%s and %s)',fileSET,fileFDT);'';...
            'Re-export all?';'';...
            'Abort script by closing figure without pressing any button'};
        opt.Default = 'No';
        opt.Interpreter = 'tex';
        button = questdlg(char(str),'Export data','Yes','No',opt);
    else
        button = 'Yes';
    end
    
    %EXPORT
    switch button
        case 'Yes'
            fprintf('Exporting:%s\n',sprintf(' %s',Files{:}))
            clear EEG
            EEG = eeg_emptyset; %init
            EEG.filename = fileSET;
            EEG.filepath = rPath;
            
            %LOAD LVM-FILE
            file = ls('*EEG.lvm');
            if exist(file,'file')~=2
                error('NO file found *EEG.lvm')
            end
            [data,info,bad] = lvm_import_zen(file);
            if isempty(data)
               error(info) 
            end
            EEG.data = data.Segment1.dataY(:,channels)';
            [noCHA,noSAM] = size(EEG.data);
            
            %set metadata
            EEG.srate  = data.Segment1.fs(1);
            EEG.nbchan = noCHA;
            EEG.pnts   = noSAM;
            EEG.trials = 1;
            EEG.xmin   = 0; %time vector starting with 0
            EEG.xmax   = (EEG.pnts-1)/EEG.srate;
            %check
            EEG = eeg_checkset(EEG);
            %save
            EEG = pop_saveset(EEG,'filename',fileSET);
        case 'No'
            fprintf('Nothing saved\n')
        otherwise
            fprintf(2,'Script aborted\n')
            return
    end
end

%% SLEEP SCORING
fprintf('\nSLEEP SCORING\n')
if ~exist('EEG','var') || ~isnumeric(EEG.data)
    EEG = pop_loadset(fullfile(rPath,fileSET));
end
%start scoring script
fprintf('START csc_eeg_plotter\n')
EEG = csc_eeg_plotter(EEG);

%% SAVE
file = fullfile(rPath,fileSET);
tmp = '';
while ~strcmp(file,tmp)
    [sFile,sPath] = uiputfile('.set','Confirm Saving',file);
    if ~ischar(sFile)
        break
    end
    tmp = fullfile(sPath,sFile);
    if ~strcmp(file,tmp)
        str = {'\fontsize{11}';...
            'Files Selection is only to confirm saving';...
            'Changing filename is not allowed!';'';...
            'You can just press ''Save'' for save scorings';
            'or ''Abort'' for decline saving';...
            };
        opt.WindowStyle = 'modal';
        opt.Interpreter = 'tex';
        waitfor(warndlg(char(str),'NOTE',opt));
    end
end
if ischar(sFile)
    tmp = EEG.data;
    EEG.data = fileFDT;
    EEG = pop_saveset(EEG,'filename',fileSET); 
    EEG.data = tmp;
    fprintf('Scoring Saved!\n')
else
    fprintf('[\bScoring NOT Saved!]\b\n')
end