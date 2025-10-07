% f_readHypnogram_raf
%------------------------------------------------------------------------
% Reads hypnogram from txt-files (exported from Sleepsign) and returns it
% for 1-s epochs.
% Stages:
%  1 for wake
%  2 for NREM
%  3 for REM
%  ... (add as many that exist)
%
%
% Thomas Rusterholz, 12 Dec 2019
%-------------------------------------------------------------------------
function [hypnogram,hdr] = f_readHypnogram_SleepSign(file)

%DEFINITIONS
Fields={... data fields {label in file, field data, function transform}
    'No.'           ,'no'       ,@str2double;... must come first!!!
    'Epoch No.'     ,'epochNo'  ,@str2double;...
    'Time'          ,'time'     ,{@(x)datenum(x,'dd.mm.yyyy HH:MM:SS'),...
    @(x)datenum(x,'dd/mm/yyyy HH:MM:SS')};...
    'Episode'       ,'episode'  ,@(x){x};...
    'Count'         ,'count'    ,@str2double;...
    'Duration(sec.)','duration' ,@str2double;...
    };
Stages={... as numeric stages used in OptoLab
    ... Any not defined stage will produce an error, for to be aware of
    1,'W' ,'Wake';...  Wake
    2,'NR','NREM';... NREM
    3,'R' ,'REM';...  REM
    4,'M' ,'Cataplexy';... all must be defined(think about)
    };
hypnogram=[]; %init, if error

%READ FILE
[fid,hdr.errorMessage]=fopen(file,'rt');
if fid==-1
    return
end
%HEADER
str=strtrim(fgetl(fid)); cnt=0; %count unknowns
str1=Fields{1,1}; n=numel(str1); noHDRlin=0; %number of header lined
delimiters={sprintf('\t'),','}; %find first one in this order!
for k=1:numel(delimiters)
    delimiter=delimiters{k}; %usually it's a tab, but I also found a comma
    if ismember(delimiter,str)
        break
    end
end
while numel(str)<n || ~strncmpi(str1,str,n)
    noHDRlin=noHDRlin+1;
    %create header
    if numel(str)>0
        tmp=regexp(str,delimiter,'split');
        tmp(cellfun(@isempty,tmp))=[];
        value=tmp;
        if ismember(tmp{1}(1),'0123456789')
            label='Date';
        else
            label=strrep(tmp{1},' ','');
            switch label
                case 'Stage'
                    label='Labels';
                case {'W','R','NR','M'}
                    label=Stages{strcmpi(Stages(:,2),label),3};
                otherwise
                    value=sprintf(' ,%s',value{2:end});
                    value(1:2)=[];
            end
        end
        if isempty(value)
            value='';
        end
        hdr.(label)=value;
    end
    %next line
    str=strtrim(fgetl(fid));
end
%DATA
labels=strtrim(regexp(str,delimiter,'split'));
str=fgetl(fid);
cnt=0;
while ischar(str)
    dat=regexp(str,delimiter,'split');
    cnt=cnt+1;
    for lab=1:size(Fields,1)
        [label1,label2,fun]=Fields{lab,:};
        tmp = dat(strcmpi(labels,label1));
        tmp(cellfun(@isempty,tmp))=[];
        if strcmpi(label1,'time')
            Fun=fun;
            for k=1:numel(Fun)
                fun=Fun{k};
                try
                    tmp = cellfun(@(x)fun(x),tmp);
                    break
                catch
                end
            end
        else
            tmp = cellfun(@(x)fun(x),tmp);
        end
        ind = (1:numel(tmp));
        if cnt==1
            data.(label2)(ind,1) = tmp;
        else
            if lab==1
                ind1=numel(data.(label2));  %must be same for all!
            end
            data.(label2)(ind+ind1,1) = tmp;
        end
    end
    str=fgetl(fid);
end
fclose(fid);
%sort data by time
[~,ind]=sort(data.time);
for lab=1:size(Fields,1)
    label2=Fields{lab,2};
    data.(label2)=data.(label2)(ind);
end

%CHECKS
%multiple timing
if diff(data.time)==0
    hdr.errorMessage='Multiple scorings at same time';
    return
end
%consitence sorted timing with duration
tmp=datevec(data.time);
tmp=etime(tmp(2:end,:),tmp(1:end-1,:))-data.duration(1:end-1);
if any(tmp~=0)
    str=sprintf(' %i',data.no(tmp~=0)+noHDRlin);
    hdr.errorMessage=char({...
        'Missing or overlapping scorings';...
        sprintf('Check line(s):%s',str);...
        });
    return
end
%for 1-s hypnogram --> scoring durations must be integer (in [s])
if ~isequal(data.duration,round(data.duration))
    hdr.errorMessage('Scoring durations must be multiple of 1s')
    return
end

%HYPNOGRAM, EPOCH LENGTH 1s
t=[0;cumsum(data.duration)];
hypnogram=NaN(t(end),1);
for k=1:numel(t)-1
    hypnogram(t(k)+1:t(k+1)) = ...
        Stages{strcmpi(Stages(:,2),data.episode{k}),1};
end
hdr = rmfield(hdr,'errorMessage');
hdr.Data=data;
hdr.Data.stageNUM=cell2mat(Stages(:,1));
hdr.Data.stageLAB=Stages(:,3);