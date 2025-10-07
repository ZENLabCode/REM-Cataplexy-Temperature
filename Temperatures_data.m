addpath('Z:\3 Markus group\Bianca\spindleDetection\functions\') %f_readHypnogram_SleepSign
addpath('Z:\OptoLab_v4.1\function\misc') %selectionList
addpath('Z:\OptoLab_v4.1\function\figures') %fig_size

clc; clear; close all
%% SELECT FILES

rPath = '';
filesHYP = fullfile(rPath,'**\*VStrn.txt');

%date formats (tries in this order)
forms = {'dd/mm/yyyy', 'dd.mm.yyyy'};


%SELCECT HYPNOGRAM FILE
tmp = dir(filesHYP);
fileHYP = fullfile({tmp.folder},{tmp.name});
tmp = cellfun(@(x)strrep(x,fileparts(fileparts(fileparts(x))),''),...
fileHYP,'uniformoutput',false);                                 

[~,ind] = selectionList(tmp,[],struct('unique',false,'sort',false));
fileHYP = fileHYP{ind};
[~,rFile] = fileparts(fileHYP);
tmp = regexp(rFile,' ','split');
name = strjoin(tmp(1:2),' ');

%READ HYPNOGRAM
[Hypnogram,hdr] = f_readHypnogram_SleepSign(fileHYP);
t0 = round(datevec(hdr.Starttime,'dd/mm/yyyy HH:MM:SS'));

t0(1:3) = 0;
n = numel(Hypnogram);
t_hyp = repmat(t0,n,1);
t_hyp(:,end) = t_hyp(:,end)+((0:n-1)');
t_hyp = datenum(t_hyp);


%SELECT SKIN FILE
rPath = fileparts(fileHYP);
tmp = dir(fullfile(rPath,'*skin.xlsx'));
tmp([tmp.isdir]) = [];
filesSKIN = fullfile({tmp.folder},{tmp.name});
fileSKIN = selectionList(filesSKIN);
fileSKIN = fileSKIN{1};
%READ SKIN FILE
[~,~,data] = xlsread(fileSKIN);
ind = find(strcmpi(data(:,1),'time'));
hdrS = data(1:ind,:);
data(1:ind,:) = [];
ind = find(cellfun(@(x)~isnumeric(x)||isnan(x),data(:,1)),1,'first');
data(ind:end,:) = [];
%define vectors
skin = cell2mat(data(:,2));
tmp = datevec(cell2mat(data(:,1)));
dt = max(diff(tmp(:,6)));
switch dt
    case 1 %everything is ok
    case 2 %interpolate
        n = numel(skin);
        skin = interp1(1:n,skin,1:0.5:n);
        skin = skin(:);
    case 0.5
        n = numel(skin);
        skin = downsample(skin,2);
        skin = skin(:);
    otherwise
        error('Check that and find a new solution')
end
n = numel(skin);
%time vector based on start time
ind = strcmpi(hdrS(:,1),'time:');
t0 = round(datevec(datenum(datevec(hdrS{ind,2})+tmp(1,:))));
t0(1:3) = 0;
t_skin = repmat(t0,n,1);
t_skin(:,end) = t_skin(:,end)+((0:n-1)');
t_skin = datenum(t_skin);

%% SELECT DEEPLABCUT FILE
rPath =  fileparts(fileHYP);
tmp   =  dir(fullfile(rPath,'*dlc.csv')); 
tmp([tmp.isdir]) = [];
filesSKINDLC = fullfile({tmp.folder},{tmp.name});
fileSKINDLC = selectionList(filesSKINDLC);
fileSKINDLC = fileSKINDLC{1};
%READ DEEPLABCUT FILE
[~,~,tmp] = xlsread(fileSKINDLC);
labelsdlc = tmp(1,:);
%remove empty labels
ind = cellfun(@(x)ischar(x)&&~isempty(x),labelsdlc);
labelsdlc(~ind) = [];
%define vectors (matrix in this case)
dlcskin = cell2mat(tmp(2:end,ind));
%interpolate
if size(dlcskin,1)<20000
    t1 = 1:size(dlcskin,1);
    t2 = 1:0.5:t1(end);
    n  = size(dlcskin,2);
    tmp = NaN(numel(t2),n);
    for k = 1:n
        tmp(:,k) = interp1(t1,dlcskin(:,k),t2,'linear');
    end
    %testing
    if false
        figure
        plot(t1,dlcskin(:,5),'o'); hold on
        plot(t2,tmp(:,5),'.');
        zoom on
        return
    end
    %assign
    dlcskin = tmp;
end


%SELECT AMB FILE
rPath = fileparts(fileHYP);
tmp = dir(fullfile(rPath,'*amb.xlsx'));
tmp([tmp.isdir]) = [];
filesAMB = fullfile({tmp.folder},{tmp.name});
fileAMB = selectionList(filesAMB);
fileAMB = fileAMB{1};
%READ AMB FILE
[~,~,data] = xlsread(fileAMB);
ind = find(strcmpi(data(:,1),'id'));
hdrA = data(1:ind,:);
data(1:ind,:) = [];
ind = find(cellfun(@(x)~isnumeric(x)||isnan(x),data(:,1)),1,'first');
if ~isempty(ind)
    data(ind:end,:) = [];
end
%define vectors
amb = cell2mat(data(:,4));
for k = 1:numel(forms)
    try
        tmp = cellfun(@(x)datevec(x,forms{k}),data(:,2),...
            'uniformoutput',false);
        err = [];
        break
    catch err
    end
end
if ~isempty(err)
    rethrow(err)
end
tmp = cat(1,tmp{:});
tmp1 = datevec(cell2mat(data(:,3)));
tmp(:,4:6) = tmp1(:,4:6);
tmp(:,1) = tmp(:,1)-tmp(1,1); 
tmp(:,2) = tmp(:,2)-tmp(1,2); 
tmp(:,3) = tmp(:,3)-tmp(1,3); 
tmp = datenum(tmp)*24*3600; %time vector in [s]
t_amb = tmp(1):1:max(tmp);
amb = interp1(tmp,amb,t_amb);
t_amb = t_amb/24/3600; %in [days]
%transpose
amb = amb(:);
t_amb = t_amb(:);


% SELECT CORE FILE
rPath = fileparts(fileHYP);
tmp = dir(fullfile(rPath,'*core.xlsx'));
tmp([tmp.isdir]) = [];
filesCORE = fullfile({tmp.folder},{tmp.name});
fileCORE = selectionList(filesCORE);
fileCORE = fileCORE{1};
%READ CORE FILE
[~,~,data] = xlsread(fileCORE);
ind = find(strcmpi(data(:,1),'Sample number'));
hdrC = data(1:ind,:);
data(1:ind,:) = [];
ind = find(cellfun(@(x)~isnumeric(x)||isnan(x),data(:,1)),1,'first');
data(ind:end,:) = [];
%define vectors
core = cell2mat(data(:,4));
for k = 1:numel(forms)
    try
        tmp = cellfun(@(x)datevec(x,forms{k}),data(:,2),...
            'uniformoutput',false);
        err = [];
        break
    catch err
    end
end
if ~isempty(err)
    rethrow(err)
end
tmp = cat(1,tmp{:});
tmp1 = datevec(cell2mat(data(:,3)));
tmp(:,4:6) = tmp1(:,4:6);
tmp(:,1) = tmp(:,1)-tmp(1,1); 
tmp(:,2) = tmp(:,2)-tmp(1,2); 
tmp(:,3) = tmp(:,3)-tmp(1,3); 
tmp = datenum(tmp)*24*3600; 
t_core = tmp(1):1:max(tmp);
core = interp1(tmp,core,t_core);
t_core = t_core/24/3600; %in [days]
%transpose
core = core(:);
t_core = t_core(:);


% SELECT BRAIN FILE
rPath = fileparts(fileHYP);
tmp = dir(fullfile(rPath,'*brain.xlsx'));
if isempty(struct2cell(tmp)) == 0
    tmp([tmp.isdir]) = [];
    filesBRA = fullfile({tmp.folder},{tmp.name});
    fileBRA = selectionList(filesBRA);
    fileBRA = fileBRA{1};
else  %IN CASE THERE'S NO TBRAIN RECORDING IN FOLDER
    rPath = 'Z:\3 Markus group\Simone\Sleep recordings TaChallenge\Mice files';
    tmp = dir(fullfile(rPath,'*brain.xlsx'));
    tmp([tmp.isdir]) = [];
    filesBRA = fullfile({tmp.folder},{tmp.name});
    fileBRA = selectionList(filesBRA);
    fileBRA = '*.xlsx';
end

% READ BRAIN FILE
[~,~,data] = xlsread(fileBRA);
brain  = cell2mat(data(2:end,strcmpi(data(1,:),'Temp')));


%% TIME SYNCH
%correct start
tmp = [t_hyp(1),t_skin(1),t_amb(1),t_core(1)];
[~,ind] = min(tmp);
switch ind(1)
    case 1
        t = t_hyp;
        n = abs(etime(datevec(t_skin(1)),datevec(t(1))));
        skin = [NaN(n,1);skin];
        dlcskin = [NaN(n,numel(labelsdlc));dlcskin];
        n = abs(etime(datevec(t_amb(1)),datevec(t(1))));
        amb = [NaN(n,1);amb];
        n = abs(etime(datevec(t_core(1)),datevec(t(1))));
        core = [NaN(n,1);core];
    case 2
        t = t_skin;
        n = abs(etime(datevec(t_hyp(1)),datevec(t(1))));
        Hypnogram = [NaN(n,1);Hypnogram];
        brain = [NaN(n,1);brain];
        n = abs(etime(datevec(t_amb(1)),datevec(t(1))));
        amb = [NaN(n,1);amb];
        n = abs(etime(datevec(t_core(1)),datevec(t(1))));
        core = [NaN(n,1);core];
    case 3
        t = t_amb;
        n = abs(etime(datevec(t_hyp(1)),datevec(t(1))));
        Hypnogram = [NaN(n,1);Hypnogram];
        brain = [NaN(n,1);brain];
        n = abs(etime(datevec(t_skin(1)),datevec(t(1))));
        skin = [NaN(n,1);skin];
        dlcskin = [NaN(n,numel(labelsdlc));dlcskin];
        n = abs(etime(datevec(t_core(1)),datevec(t(1))));
        core = [NaN(n,1);core];
    case 4
        t = t_core;
        n = abs(etime(datevec(t_hyp(1)),datevec(t(1))));
        Hypnogram = [NaN(n,1);Hypnogram];
        brain = [NaN(n,1);brain];
        n = abs(etime(datevec(t_amb(1)),datevec(t(1))));
        amb = [NaN(n,1);amb];
        n = abs(etime(datevec(t_skin(1)),datevec(t(1))));
        skin = [NaN(n,1);skin];
        dlcskin = [NaN(n,numel(labelsdlc));dlcskin];
end
%correct end
answer = questdlg('Start time?','T recording','9','10','20','10');
if strcmpi(answer,'10')
    str = 10;
    [~,ind1] = min(abs(t-datenum([0,0,0,10,0,0])));
   
    ind = (1:8*3600)+ind1-1;
elseif strcmpi(answer,'20')
    str = 20;
    [~,ind1] = min(abs(t-datenum([0,0,0,20,0,0])));
    
    ind = (1:8*3600)+ind1-1;
else
    str = 9;
    [~,ind1] = min(abs(t-datenum([0,0,0,9,0,0])));
    
    ind = (1:8*3600)+ind1-1;
end
% if numel(Hypnogram)<ind(end)
Hypnogram(numel(Hypnogram)+1:ind(end)) = NaN;
brain(numel(brain)+1:ind(end)) = NaN;
% end
Hypnogram = Hypnogram(ind);
brain = brain(ind);
if numel(skin)<ind(end)
    skin(numel(skin)+1:ind(end)) = NaN;
end


skin = skin(ind);
if size(dlcskin,1)<ind(end)
    dlcskin(size(dlcskin,1)+1:ind(end),:) = NaN;
end
dlcskin = dlcskin(ind,:);
if numel(amb)<ind(end)
    amb(numel(amb)+1:ind(end)) = NaN;
end
amb = amb(ind);
if numel(core)<ind(end)
    core(numel(core)+1:ind(end)) = NaN;
end
core = core(ind);
if numel(t)<ind(end)
    t(numel(t)+1:ind(end)) = NaN;
end
t = t(ind);
%%
%Fitting TSkin
skinNAN = skin;
skinNAN(skin<= 32) = nan;
skinfit = movmean(skinNAN,150,'omitnan');
skinNAN(skin>= 39) = nan;
dlcskinfit = movmean(dlcskin, 150, 'omitnan');
h = Hypnogram;
h(Hypnogram>4) = nan;
Hypnogram = h;


%% MAKE SUMMARY TABLE
labels = {'Clock Hour', ...
    'Sleep Stage',...
    'Core T.',...
    'Skin T.',...
    'Ambient T.',...
    'Brain T.'};
ctk = [4, 7, 9];
labelsdlc = labelsdlc(ctk);
labels = [labels, labelsdlc];   

Table  = cell(0,numel(labels));

M = cell(numel(t),numel(labels));
M(:,strcmpi(labels,'Clock Hour')) = num2cell(t);
M(:,strcmpi(labels,'Sleep Stage'))   = num2cell(Hypnogram);
M(:,strcmpi(labels,'Core T.'))   = num2cell(core);
M(:,strcmpi(labels,'Skin T.'))   = num2cell(skinfit);
M(:,strcmpi(labels,'Ambient T.'))   = num2cell(amb);
M(:,strcmpi(labels,'Brain T.'))   = num2cell(brain);
 for tb = 1:numel(labelsdlc)
     M(:,strcmpi(labels,labelsdlc(tb))) = num2cell(dlcskinfit(:,ctk(tb)));
 end
%append
Table = [Table;labels;M;cell(1,numel(labels))];

%SAVE
sFile = strcat(name, ' summaryyy');
rPath = fileparts(fileHYP);
[sFile,rPath] = uiputfile('*.xlsx','Save Table',fullfile(rPath,sFile));
if ischar(sFile)
    sname = fullfile(rPath,sFile);
    if exist(sname,'file')==2
        delete(sname)
    end
    xlswrite(sname,Table)
    fprintf('Saved: %s\n',sname)
    try
        xls_cellFit(sname)
    catch
    end
else
    fprintf('Table NOT saved\n')
end