% Analyze_SleepBouts
% Analyzes sleep bouts per trial and bin from txt-file.
%-----------------------------------------------------------------------
clc; clear; close all; fclose all;
addpath('Z:\OptoLab_v4.1\function\excel')

%PARAMETERS
%----------
%START PATH 

startPath=[''] ;

%ANALYZIS
durTRI=2; %trial duration [h]
trials=[1,2,3,4];  %trials to analyze
noBIN=2;  %number of bins per trial 
maxInterval=1000; %in [s]

%MOUSE/GROUP
posGroup=1;
posMouse=2;


%OPTIONS (true or false)
saveXLS=true; 
autoCellFit=true; 

%READ FILE INFOS 
file.splitSTR=sprintf('\t'); 
file.columns={... {column variiable in script, label in read file}
    'colNUM','No.';... 
    'colEPO','Epoch No.';...
    'colTIM','Time';...
    'colEPI','Episode';...
    'colCNT','Count';...
    'colDUR','Duration(sec.)';...
    };

file.dateFormat='dd/mm/yyyy HH:MM:SS';
file.dateFormat2='dd.mm.yyyy HH:MM:SS';

%MAIN PROGRAM
%------------
%SELECT FILE
if exist(startPath,'file')==2
    fname=startPath;
else
    [rFile,rPath]=uigetfile('.txt','Select A File',startPath);
    if ~ischar(rFile)
        fprintf('No File Selected!\n')
        return
    end
    fname=fullfile(rPath,rFile);
end
%get group and mouse name
[~,rFile]=fileparts(fname);
tmp=regexp(rFile,' ','split');
group=tmp{posGroup};
mouse=tmp{posMouse};

%READ FILE
fid=fopen(fname,'r');

%find label line
lab1=file.columns{1,2}; %1st label
str=strtrim(fgetl(fid));
while ~strncmpi(str,lab1,numel(lab1))
    str=strtrim(fgetl(fid));
    if isnumeric(str)
        error('Wrong File, label ''%s'', not found!',lab1)
    end
end
labelsF=regexp(str,file.splitSTR,'split'); 

for col=1:size(file.columns,1)
    [labVAR,labSTR]=file.columns{col,:};
    ind=find(strcmpi(labelsF,labSTR),1,'first');
    if isempty(ind)
        error('Label ''%'' is not available!',labSTR)
    end
    assignin('base',labVAR,ind);
end
%read data
while true
    str=fgetl(fid);
    if ~ischar(str) %end of file
        break
    end
    tmp=regexp(str,file.splitSTR,'split');
    if all(cellfun(@isempty,tmp))
        continue
    end
    if ~exist('Data','var') 
        Data=cell(0,numel(tmp));
    end
    Data=[Data;tmp];
end
fclose(fid);

%DATA VARIABLES
noDAT=size(Data,1);
vecDuration=cellfun(@str2double,Data(:,colDUR));
vecEpisode=Data(:,colEPI);
try
    timeSTA=cellfun(@(x) datenum(x,file.dateFormat),Data(:,colTIM));
catch
    timeSTA=cellfun(@(x) datenum(x,file.dateFormat2),Data(:,colTIM));
end
timeEND=datenum(datevec(timeSTA)+[zeros(noDAT,5),vecDuration]);
%durations and number of ...
T0=datevec(timeSTA(1));
answer = questdlg(...
    sprintf('Is the start time OK?  (thank Simone for this nuisance...)\n%s',datestr(T0,'HH:MM:SS')),...
	'Start Time check!','Yes','No','Yes');
if strcmpi(answer,'No')
    T0(5:6) = 0;
end
TE=datevec(timeEND(end));
durTOT=etime(TE,T0); %total duration
durTRI=durTRI*3600; %in seconds
durBIN=durTRI/noBIN;
noTRItot=ceil(durTOT/durTRI); 
if isempty(trials)
    trials=1:noTRItot;
end
noTRI=numel(trials);


%TABLE 1, duration sheet
episodes={'W','R','NR','M'};
noEPI=numel(episodes);
tmp={episodes,'uniformoutput',false};
labels=['Mouse','Group','Condition','Day','TNZGroup','Trial','Bin',...
    'REM latency','Time',...
    cellfun(@(x) sprintf('Duration %s',x),tmp{:}),...
    cellfun(@(x) sprintf('No Bouts %s',x),tmp{:}),...
    cellfun(@(x) sprintf('Mean Bout Dur %s',x),tmp{:}),...
    cellfun(@(x) sprintf('SD Bout Dur %s',x),tmp{:}),...
    ];
Table=cell(noTRI*noBIN,numel(labels)); 
Table(:,1)={mouse};
Table(:,2)={group};
%for additional analysis
D.durarion=NaN(noTRI,noBIN,noEPI);
D.boutDUR=cell(noTRI,noBIN,noEPI);
D.boutSTA=cell(noTRI,noBIN,noEPI);
D.boutEND=cell(noTRI,noBIN,noEPI);
%trial loop
for tri=1:noTRI
    trial=trials(tri);
    %bin loop
    for bin=1:noBIN
        %start and end time of trial/bin
        t=(trial-1)*durTRI+(bin-1)*durBIN; %time from T0
        TS=datenum(T0+[zeros(1,5),t]);
        TE=datenum(T0+[zeros(1,5),t+durBIN]);
        
        %FILL TABEL 1/2
        row=bin+(tri-1)*noBIN; %table row
        Table{row,strcmpi(labels,'Trial')}=trial;
        Table{row,strcmpi(labels,'Bin')}=bin;
        Table{row,strcmpi(labels,'Time')}=datestr(TS,'dd.mm.yyyy HH:MM:SS');
        %REM latency
        ind=find(timeSTA>=TS&timeSTA<TE&strcmpi(vecEpisode,'r'),1,'first');
        if isempty(ind)
            Table{row,strcmpi(labels,'REM latency')}=NaN;
        else
            Table{row,strcmpi(labels,'REM latency')}=...
                etime(datevec(timeSTA(ind)),datevec(TS));
        end
        
        %episode loop
        for epi=1:noEPI
            episode=episodes{epi};
            index = timeEND>TS & timeSTA<TE & strcmpi(vecEpisode,episode);
            
            %TIME/BOUT CORRECTION
            remLat=NaN;
            tCorr=0; %initialization
            ind=find(index,1,'first');
            if ~isempty(ind) 
                dt1=etime(datevec(timeSTA(ind)),datevec(TS));
                dt2=etime(datevec(timeEND(ind)),datevec(TS));
                if dt2<0
                    error('upps');
                end
                if dt1<0 %bout starts in previous bin
                    if abs(dt1)<dt2
                        tCorr=tCorr+dt1; %time correction
                    else
                        index(ind)=false; %exclude bout
                        tCorr=tCorr+dt2; %time correction
                    end
                end
                %time duration correction bin end
                ind=find(index,1,'last');
                dt1=etime(datevec(timeSTA(ind)),datevec(TE));
                dt2=etime(datevec(timeEND(ind)),datevec(TE));
                if dt1>0
                    error('upps');
                end
                dtEND=0;
                if dt2>0 %bout exceeds this bin
                    if abs(dt1)<dt2
                        index(ind)=false; %exclude bout
                        tCorr=tCorr-dt1; %time correction
                    else
                        tCorr=tCorr-dt2; %time correction
                    end
                end
            end
            
            %FILL TABLE 2/2
            duration=sum(vecDuration(index))+tCorr;
            boutDurations=vecDuration(index);
            Table{row,strcmpi(labels,sprintf('Duration %s',episode))}=...
                duration;
            Table{row,strcmpi(labels,sprintf('No Bouts %s',episode))}=...
                sum(index);
            Table{row,strcmpi(labels,sprintf('Mean Bout Dur %s',episode))}=...
                mean(boutDurations);
            Table{row,strcmpi(labels,sprintf('SD Bout Dur %s',episode))}=...
                std(boutDurations);
            %fill data
            D.durarion(tri,bin,epi)=duration;
            D.boutDUR(tri,bin,epi)={boutDurations};
            D.boutSTA(tri,bin,epi)={timeSTA(index)};
            D.boutEND(tri,bin,epi)={timeEND(index)};
        end %episode loop
    end %bin loop
end %trial loop
TAB{1}=[labels;Table];


%TABLE 2, summary sheet
tmp={episodes,'uniformoutput',false};
labels=['Mouse','Group','Condition','Day','TNZGroup','Bin',...
    'Mean REM latency',...
    cellfun(@(x) sprintf('Mean Duration %s',x),tmp{:}),...
    cellfun(@(x) sprintf('No Bouts %s',x),tmp{:}),...
    cellfun(@(x) sprintf('Mean Bout Dur %s',x),tmp{:}),...
    cellfun(@(x) sprintf('SD Bout Dur %s',x),tmp{:}),...
    ];
Table=cell(noBIN,numel(labels)); %initialization
Table(:,1)={mouse};
Table(:,2)={group};
%bin loop
for bin=1:noBIN
    %FILL TABEL 1/2
    row=bin; %table row
    Table{row,strcmpi(labels,'Bin')}=bin;
    
    %REM latency
    tmp=TAB{1};
    colBIN=strcmpi(tmp(1,:),'Bin');
    colLAT=strcmpi(tmp(1,:),'REM latency');
    bins=cell2mat(tmp(2:end,colBIN));
    reml=cell2mat(tmp(2:end,colLAT));
    Table{row,strcmpi(labels,'Mean REM latency')}=...
        nanmean(reml(bins==bin));
    
    %episode loop
    for epi=1:noEPI
        episode=episodes{epi};
        
        %FILL TABLE 2/2
        boutDurations=vertcat(D.boutDUR{:,bin,epi});
        Table{row,strcmpi(labels,sprintf('Mean Duration %s',episode))}=...
            sum(D.durarion(:,bin,epi))/noTRI;
        Table{row,strcmpi(labels,sprintf('No Bouts %s',episode))}=...
            numel(boutDurations);
        Table{row,strcmpi(labels,sprintf('Mean Bout Dur %s',episode))}=...
            mean(boutDurations);
        Table{row,strcmpi(labels,sprintf('SD Bout Dur %s',episode))}=...
            std(boutDurations);
    end %episode loop
end %bin loop
TAB{2}=[labels;Table];

%SAVE
%savename/sheet names
[rPath,rFile,rExt]=fileparts(fname); 
sname=fullfile(rPath,[rFile,'.xlsx']);
if exist(sname,'file')==2
    delete(sname)
end
sheets={'Duration','Summary'};
%save
warning('off','MATLAB:xlswrite:AddSheet'); 
noTAB=numel(TAB);
for k=1:noTAB
    if saveXLS
        xlswrite(sname,TAB{k},sheets{k})
        if k==noTAB
            if autoCellFit
                xls_cellFit(sname);
            end
            xls_deleteSheets(sname,sheets,-1);
            fprintf('SAVED: "%s"\n',sname)
        end
    else
        disp(TAB{k})
        if k==noTAB
            fprintf('[\bData NOT saved]\b\n')
        end
    end
end
warning('on','MATLAB:xlswrite:AddSheet')