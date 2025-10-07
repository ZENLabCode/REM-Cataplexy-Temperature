clc; clear; close all; fclose all; 

addpath('Z:\OptoLab_v4.1\function\misc') %selectionList
% addpath(genpath('Z:\OptoLab_v4.1\function')) % needed
tmp = genpath('Z:\EEGlandArmandsScoringScript\EEGlabtools\eeglab');
addpath(tmp,'-end');
tmp = regexp(tmp,';','split');
ind = contains(lower(tmp),'octavefunc');
rmpath(strjoin(tmp(ind),';'))

rFiles = ['*ScoredEEG.fdt'];
channels = {'EEG1','EEG2','EMG'}; %data channel names
Calcium = {'Ca',1,480+[-10,10]};

lastSname = [channels{end},'.mat'];

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%CHECK
if numel(channels)~=numel(unique(lower(channels)))
    error('Channels names must be unique (case in-sensitive)')
end
if numel(unique(lower(Calcium(:,1))))~=size(Calcium,1)
    error('Calcium channel names must be unique')
end

%READ FILES
if iscell(rFiles)
    ind = cellfun(@(x)exist(x,'file')~=2,rFiles);
    if any(ind)
        fprintf('[\bFiles removed (do NOT exist)!]\b\n')
        fprintf(' - %s\n',rFiles{ind})
        rFiles(ind) = [];
    end
    if numel(rFiles)==0
        fprintf(2,'No File left!\n')
        return
    end
elseif ischar(rFiles)
    if exist(rFiles,'file')==2
        rFiles = {rFiles};
    else
        [stat,list] = dos(sprintf('dir "%s" /S/B/A:-H-D',rFiles));
        if stat==1
            rFiles = {};
            fprintf(2,'No File found!\n')
            return
        else
            %split list
            list = regexprep(list,sprintf('\r'),'');
            list = regexp(list,newline,'split');
            list(cellfun(@isempty,list)) = [];
            %add path if missing (just in case is missing)
            if ~isempty(list) && isempty(fileparts(list{1}))
                list = fullfile(fileparts(rFiles),list);
            end
            %exclude dot-files (hidden mac-files)
            ind = cellfun(@(x)x(find(x==filesep,1,'last')+1)=='.',list);
            list(ind) = [];
            %select files
            rFiles = list;
            if numel(rFiles)>0
                rFiles = selectionList(rFiles,lastSname);
                if numel(rFiles)==0
                    fprintf(2,'No File selected!\n')
                    return
                end
            end
        end
    end
else
    error('Class ''%s'' not supported for variable rFiles',class(rFiles))
end

%number of ...
noFIL = numel(rFiles);
noCHA = numel(channels);
noCAL = size(Calcium,1);

%FILE LOOP
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    [rPath,rFile,rExt] = fileparts(rFiles{fil});
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath)
    
    %FILENAMES
    file.set = [rFile,'.set'];
    tmp = dir(fullfile(rPath,'*_EEG.lvm'));
    if numel(tmp)~=1
        error('check that')
    end
    file.lvm = tmp.name;
    tmp = dir(fullfile(rPath,'*_Ca.lvm'));
    if numel(tmp)==1
        file.cal = tmp.name;
    elseif numel(tmp)>1
        error('check that')
    else
        file.cal = '';
    end
    %print out
    fprintf('%s Files\n',indent)
    fields = fieldnames(file);
    for k = 1:numel(fields)
        field = fields{k};
        fprintf('%s   %s-file : %s\n',indent,field,file.(field))
    end
           
    %TIME STARTS
    clear T0
    %data (lvm-file)
    T0.dat = NaN(1,6); %init
    fid = fopen(fullfile(rPath,file.lvm),'r');
    while any(isnan(T0.dat))
        lin = strtrim(fgetl(fid));
        strs = regexp(lin,sprintf('\t'),'split');
        if strcmpi(strs{1},'Date')
            strs = regexp(strs{2},'[0-9]*','match');
            T0.dat(1:3) = cellfun(@str2double,strs);
        end
        if strcmpi(strs{1},'Time')
            strs = regexp(strs{2},'[.0-9]*','match');
            T0.dat(4:6) = cellfun(@str2double,strs);
        end
    end
    fclose(fid);
    %calcium (lvm-file)
    if ~isempty(file.cal)
        T0.cal = NaN(1,6); %init
        fid = fopen(fullfile(rPath,file.cal),'r');
        while any(isnan(T0.cal))
            lin = strtrim(fgetl(fid));
            strs = regexp(lin,sprintf('\t'),'split');
            if strcmpi(strs{1},'Date')
                strs = regexp(strs{2},'[0-9]*','match');
                T0.cal(1:3) = cellfun(@str2double,strs);
            end
            if strcmpi(strs{1},'Time')
                strs = regexp(strs{2},'[.0-9]*','match');
                T0.cal(4:6) = cellfun(@str2double,strs);
            end
        end
        fclose(fid);
    end
    %hypnogram
    T0.hyp = T0.dat;
    
    %LOAD EEG & HYPNOGRAM
    EEG = pop_loadset(fullfile(rPath,file.set));
    fs = EEG.srate;
    [noCHA,noSAM] = size(EEG.data);
    if noCHA~=numel(channels)
        fprintf(['%s [\bdefined channels N = %i, but set-file has %i ',...
            'channels]\b\n'],indent,numel(channels),noCHA);
        continue
    end
    
    %HYPNOGRAM
    fprintf('%s Hypnogram\n',indent)
    strs = {sprintf('Exported by %s.m, %s',mfilename,date);'';...
        'file  : original scoring with EEGLAB';...
        'index : If not empty, hypnogram is cut';...
        '        remaining: index(1):index(2)';...
        };
    clear H
    H.info.info = char(strs);
    H.info.file = fullfile(rPath,file.set);
    H.info.index = [];
    H.fs = fs;
    H.Hypnogram = NaN(1,noSAM);
    events = [{'init',0,NaN};EEG.csc_event_data];
    for k = 1:size(events,1)-1
        [~,t1,stage] = events{k,:};
        [~,t2,~] = events{k+1,:};
        ind1 = round(t1*fs)+1;
        ind2 = round(t2*fs);
        if ind2>noSAM && noSAM-ind2>1
            error('check that')
        end
        H.Hypnogram(ind1:ind2) = stage;
    end
    %cut
    indSTA = find(~isnan(H.Hypnogram),1,'first');
    indEND = find(~isnan(H.Hypnogram),1,'last');
    indCUT = indSTA:indEND;
    if indSTA~=1 || indEND~=noSAM
        H.info.index = [indSTA,indEND];
    end
    H.Hypnogram = H.Hypnogram(indCUT);
    T0.hyp = datevec(datenum(T0.hyp+[0 0 0 0 0 (indSTA-1)/fs]));
    duration = numel(H.Hypnogram)/fs; %[s]
    %print out / save
    fprintf('%s   tot duration : %i s (%g min)\n',indent,...
        noSAM/fs,noSAM/fs/60)
    fprintf('%s   cut duration : %i s (%g min)\n',indent,...
        duration,duration/60)
    save(fullfile(rPath,'Hypnogram.mat'),'-struct','H')
    fprintf('%s   saved : Hypnogram.mat\n',indent)
    
    %DATA CHANNELS
    fprintf('%s Channel Data\n',indent)
    %save
    for cha = 1:noCHA
        strs = {sprintf('Exported by %s.m, %s',mfilename,date);'';...
            'file    : fdt-file';...
            'index   : If not empty, data is cut';...
            '          remaining: index(1):index(2)';...
            'channel : channel index of bin-file';...
            };
        clear D
        D.info.info = char(strs);
        D.info.file = fullfile(rPath,file.set);        
        D.info.index = H.info.index;
        D.info.channel = cha;
        D.SampRate = fs;
        D.resampled_data_mV = EEG.data(cha,indCUT);
        D.stimTimes = [];
        %save
        channel = [channels{cha},'.mat'];
        save(fullfile(rPath,channel),'-struct','D')
        fprintf('%s   saved : %s\n',indent,channel)
    end
    
    %CALCIUM CHANNEL
    fprintf('%s Calcium Data\n',indent)
    if isempty(file.cal)
        fprintf('%s   [\bNot available]\b\n',indent)
    else
        %read data
        data = lvm_import_zen(fullfile(rPath,file.cal));
        fs   = unique(data.Segment1.fs);
        %time frame index
        ind = (1:duration*fs)+round(etime(T0.hyp,T0.cal)*fs);
        ind(ind>size(data.Segment1.dataY,1)) = [];  %%% in case Ca file cut differently than Hyp file
        %calcium channel loop
        for cal = 1:noCAL
            [channel,num,fRangeDemod] = Calcium{cal,:};
            strs = {sprintf('Exported by %s.m, %s',mfilename,date);'';...
                'file    : lvm-file';...
                'index   : If not empty, data is cut';...
                '          remaining: index(1):index(2)';...
                'channel : calcium channel index of bin-file';...
                'signal  : demodulated signal (bandpower)';...
                };
            clear C
            C.info.info = char(strs);
            C.info.file = fullfile(rPath,file.cal);
            C.info.index = [];
            if ind(1)~=1 || ind(end)~=noSAM
                C.info.index = [ind(1),ind(end)];
            else
                C.info.index = [];
            end
            C.info.channel = num;
            C.info.signal.unit   = 'V^2';
            C.info.signal.fRange = fRangeDemod;
            %sampling rate demodulated signal
            C.SampRate = 10; %resulting/desired sample rate
            fprintf('%s   Demodulation f-Range : %g %g\n',...
                indent,fRangeDemod)
            fprintf('%s   Sampling Rate (set)  : %g\n',indent,C.SampRate)
            %demodulated signal
            signal = data.Segment1.dataY(ind,num); %cropped, in [V]
            
                             
            n1 = fs/C.SampRate;
            n2 = floor(numel(signal)/n1); %o samples
            if n1~=round(n1) %must be integer !!!
                error('crap, this does not work')
            end
            if n2~=floor(n2) %should not be the caes
                fprintf('%s   [\bDuration only %i s (%g min)]\b\n',...
                    indent,n2/C.SampRate,n2/C.SampRate/60)
                n2 = floor(n2);
            end
            signal = reshape(signal(1:n1*n2),[n1,n2]);
            C.demodulated_signal = bandpower(signal,fs,fRangeDemod);
            %save
            save(fullfile(rPath,[channel,'.mat']),'-struct','C')
            fprintf('%s   saved : %s\n',indent,[channel,'.mat'])
        end %calcium channel loop
    end
end