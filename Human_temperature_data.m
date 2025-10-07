clc; clear; close all;
addpath('.\functions')

%PARAMETERS
%----------
%FILES (events, files to read, ...)
tmp = {'Cataplexy'};
Events = tmp{1};
%paths
rPathXLS = ''; %read path data
sPathFIG = ''; %save path
%ids & filenames
switch lower(Events)
    case 'cataplexy'
        IDS    = {''};

        Fnames = cell(numel(IDS),2); %init
            for k = 1:size(Fnames,1)
            [~,id] = fileparts(Fnames{k,2});
            id = split(id,'_');
            IDS(k) = id(2);
        end
    otherwise
        error('Events ''%s'' not implemented')
end
Fnames(:,1) = fullfile(rPathXLS,cellfun(@(x)sprintf('%s.xlsx',x),IDS(:),...
    'uniformoutput',false)); %xls-files data

%OPTIONS FIGURE AND AXES
opt.fig.save = false; %save figure or not 
opt.fig.saveFuns = {
    @(h,name)saveas(h,name);...
    @(h,name)print(h,name,'-r300','-dpng','-painters');...
    };
opt.plot.meanData = true; %plot or not, mean data across individuals
opt.plot.indData  = false; %plot or not, individual plot
%axes
tmp = '';
opt.axi.data = {... {1st row, 2nd rows, title, ylabel},
    'Temperature' ,{{'Upper Arm'},{}} ,'Upper Arm','';...
    'Temperature' ,{{'Apical'}   ,{}} ,'Apical'   ,'';...
    'Temperature' ,{{'Clavicula'},{}} ,'Clavicula','';...
    'Temperature' ,{{'Wrist'}    ,{}} ,'Wrist'    ,'';...
    'Temperature' ,{{'Wrist'}    ,{'Apical','Clavicula'}},...
        'Distal - Proximal',tmp;...
    'Temperature' ,{{'Core'}     ,{}} ,'Core'     ,tmp;...;...
    'Temperature' ,{{'Ambient'}  ,{}} ,'Ambient'  ,tmp;...;...
    };

opt.axi.unit = 'pixel'; %unit for axes size specification
opt.axi.dx = [85,0,150]; %spaces left, between and righ of axes
opt.axi.dy = [55,5,55];  %spaces bottom, between top of axes
opt.axi.length = 550; %axes length
opt.axi.height = 150; %axes height

%OPTIONS PLOT
dt = 30; %[s] time steps in xls-files (will be checked)
opt.plot.tUnit   = 'min'; %time unit in figures --> all time infos in units
opt.plot.tFactor = 1/60; %factor to [s] for corresponding tUnit
opt.plot.margin  = []; %set later, see below (in tUnit)
opt.plot.color.ids = lines(numel(IDS)); %colors per ID 
opt.plot.color.mean = 'k'; %color mean and std region
opt.plot.color.stat = 'r'; %color statistic, min/max line and region 

tmp = {'change', 'diff','zscore','absolute'};
opt.plot.fun = tmp{1};
opt.plot.margin  = [-20,10]; %plot margin

%OPTIONS EVENTS
% - clockTime : selects events within clock range 
% - timeLock  : excludes events when starting in less time after
%               prevoius one ends
%               (to prevent influence of preveous event in analysis before)
opt.event.clockTime = [9,20]; %clock range [h] (tataplexies to use)
opt.event.timeLock  = 2*10; %in [opt.plot.tUnit]

%MAIN SCRIPT
%-----------
srciptName = sprintf('%s (%s)',mfilename,Events);
fprintf('%s\n%s\n',srciptName,repmat('-',size(srciptName)));
%check
ind = ~cellfun(@(x)exist(x,'file')||isempty(x),Fnames);
if any(ind(:))
    error(sprintf('File not found:\n%s',...
        sprintf(' - %s\n',Fnames{ind})))
end

%INIT
%number of ...
noFIL = size(Fnames,1);
noAXI = size(opt.axi.data,1);
%figure properties 
Fig.position = [0, 0, sum(opt.axi.dx([1,3]))+opt.axi.length,...
    sum(opt.axi.dy([1,3]))+noAXI*opt.axi.height+(noAXI-1)*opt.axi.dy(2)];

%READ DATA
nnFIL  = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
DATA   = struct('info','','id','','t',[],'cols',{},'event',[],'Data',[]);
for fil = 1:noFIL
    [rPath,ID] = fileparts(fnameXLS);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,ID)
    if exist(fnameXLS,'file')~=2
        fprintf(2,'\b - File does not exist\n')
        continue
    end
    end
    
    %READ/PREPARE DATA
    [~,~,raw] = xlsread(fnameXLS);
    Labels = raw(1:2,:);
    Labels(~cellfun(@ischar,Labels)) = {''};
    Labels = cellfun(@strtrim,Labels,'UniformOutput',false); 
    raw(1:2,:) = [];
    noSAM = size(raw,1);
    %data matrix
    Data = NaN(noSAM,noAXI);
    for axi = 1:noAXI
        [lab1, labs2] = opt.axi.data{axi,1:2};
        ind1 = strcmpi(Labels(1,:),lab1);
        ind2 = ismember(lower(Labels(2,:)),lower(labs2{1}));
        ind  = ind1 & ind2;
        if sum(ind)~=numel(labs2{1})
            error('uups')
        end
        dat = nanmean(cell2mat(raw(:,ind)),2);
        if ~isempty(labs2{2})
            ind2 = ismember(lower(Labels(2,:)),lower(labs2{2}));
            ind  = ind1 & ind2;
            if sum(ind)~=numel(labs2{2})
                error('uups')
            end
            dat = dat-nanmean(cell2mat(raw(:,ind)),2);
        end
        Data(:,axi) = dat;
    end
    
    %TIME VECTOR
    column = strcmpi(Labels(2,:),'Date');
    assert(sum(column)==1,'Date column not found')
    dates = cell2mat(cellfun(@(x)datevec(x,'dd.mm.yyyy'),raw(:,column),...
        'UniformOutput',false));
    column = strcmpi(Labels(2,:),'Time');
    assert(sum(column)==1,'Time column not found')
    tmp = datevec(cell2mat(raw(:,column)));
    dates(:,4:6) = tmp(:,4:6);
    if dt~=mean(etime(dates(2:end,:),dates(1:end-1,:)))
        eror('Check, dt should be 30 but is %g!',dt)
    end
    %time index to select for plotting 
    margin = opt.plot.margin/opt.plot.tFactor/dt; %[s]
    ind0   = floor(margin(1)-1):ceil(margin(2)+1); %+-1 for 2 more
    switch opt.plot.fun
        case {'change','zscore','absolute'}
            t0 = ind0*dt*opt.plot.tFactor;
        case 'diff'
            t0 = (ind0(1:end-1)*dt+dt/2)*opt.plot.tFactor;
        otherwise
            error('ggrhhh')
    end
    
    %EVENTS (cataplexy)
    switch lower(Events)
        case 'cataplexy'
            column = ... find data column
                strcmpi(Labels(1,:),'Activity') & ...
                strcmpi(Labels(2,:),'Behavior');
            assert(sum(column)==1,'Activity data not found')
            tmp = strcmpi(raw(:,column),'cataplexy');
            if sum(tmp)==0
                fprintf('\b - [\bNo %s Found]\b\n',Events)
                continue
            end
            %events start/end indices
            indSTA = find( [tmp;NaN]==1 & [0;tmp]==0);
            indEND = find( [NaN;tmp]==1 & [tmp;0]==0 ) - 1;
        
    if any((indEND-indSTA)<0) 
        error('grrhhh')
    end
    %exclusions I (do this first)
    dif = (indSTA(2:end)-indEND(1:end-1))*dt*opt.plot.tFactor;
    ind = [false; dif < opt.event.timeLock];
    if any(ind)
        fprintf(['%s [\b%s excluded (N = %i): ',...
            'starting in less than %g %s after previous end]\b\n'],...
            indent,Events,sum(ind),opt.event.timeLock, opt.plot.tUnit)
        indSTA(ind) = [];
        indEND(ind) = [];
    end
    %exclusions II
    tmp = dates(indSTA,:)*[0,0,0,1,1/60,1/3600]'; %[h]
    ind = tmp<opt.event.clockTime(1) | tmp>opt.event.clockTime(2);
    if any(ind)
        fprintf(['%s [\b%s excluded (%i of %i): ',...
            'out of clock-time %05.2f - %05.2f]\b\n'],...
            indent,Events,sum(ind),numel(ind),opt.event.clockTime)
        indSTA(ind) = [];
        indEND(ind) = [];
    end
    noEVT = numel(indSTA);
    
    %DATA MATRIX
    data = NaN(numel(t0),noAXI,noEVT);
    for evt = 1:noEVT
        ind = indSTA(evt)+ind0;
        in2 = find(ind>0 & ind<noSAM);
        for axi = 1:noAXI
            dat = Data(ind(in2),axi);
            switch opt.plot.fun
                case 'absolute'
                    data(in2,axi,evt) = dat; 
                case 'change'
                    data(in2,axi,evt) = dat-dat(t0==0);
                case 'diff'
                    data(in2(1:end-1),axi,evt) = -diff(dat);
                case 'zscore'
                    zs = 1;
                    switch zs
                        case 1 %zscore stat range
                            mi = min(cellfun(@min,opt.stat.margin(:,2)));
                            ma = max(cellfun(@max,opt.stat.margin(:,2)));
                            tmp = dat(t0(in2)>=mi & t0(in2)<=ma);
                            dat = (dat-nanmean(tmp))/nanstd(tmp);
                        case 2 %zscore plot range
                            dat  = (dat-nanmean(dat))/nanstd(dat);
                        case 3 %zscore all data
                            tmp = Data(:,axi);
                            tmp = (tmp-nanmean(tmp))/nanstd(tmp);
                            dat = tmp(ind(in2));
                        otherwise
                            error('ggrhhh')
                    end
                    data(in2,axi,evt) = dat;
                otherwise
                    error('ggrhhh')
            end
        end
    end
    
    %APPEND DATA
    DATA(fil).info  = 'Data: t * cols * event (datenum)';
    DATA(fil).id    = ID;
    DATA(fil).t     = t0(:);
    DATA(fil).cols  = opt.axi.data(:,3)';
    DATA(fil).event = datenum(dates(indSTA,:));
    DATA(fil).Data  = data;
end

%% AVERAGE DATA ACROSS INDIVIDUAL EVENTS
N_meanTot = noFIL;
for fil = 1:noFIL
    data = DATA(fil);
    if fil==1
        mDATA.info = 'Data: t * cols * id';
        mDATA.t    = data.t;
        mDATA.cols = data.cols;
        mDATA.ids  = {DATA.id}';
        mDATA.Data = NaN(numel(mDATA.t),numel(mDATA.cols),...
            numel(mDATA.ids));
    end
    dat = nanmean(data.Data,3);
    mDATA.Data(:,:,fil) = dat;
    if all(isnan(dat))
        N_meanTot = N_meanTot-1;
    end
end

%SAVE TABLE
sFile = sprintf('trans%s_Human_Mean_%s',Events,opt.plot.fun);
sname = fullfile(sPathFIG,sFile);
if exist([sname,'.xlsx'],'file')==2
    delete([sname,'.xlsx'])
end
warning off
for k = 1:numel(mDATA.cols)
    dat = squeeze(mDATA.Data(:,k,:));
    tab = ['t',num2cell(mDATA.t');mDATA.ids,num2cell(dat')];
    xlswrite([sname,'.xlsx'],tab,mDATA.cols{k})
end
warning on

%% FIGURE, MEAN DATA
clc; close all; %%%
if opt.plot.meanData && N_meanTot>1
    clc; close all
    
    
    %DATA SAME FOR ALL AXES
    t    = mDATA.t;
    indR = numel(t):-1:1; 
    X  = [t; t(indR)];
    Z  = zeros(size(X));
    xl = opt.plot.margin;
    
    %FIGURE / AXES
    ax.hig = 110; ax.wid = 375;
    ax.dx  = [75,0,140]; ax.dy  = [70,7,40];
    posF = [0,0,ax.wid+sum(ax.dx),...
        noAXI*ax.hig+[1,noAXI-1,1]*ax.dy'];
    hf  = figure('unit','pixel','position',posF,'visible','off');
    if ~isequal(posF,get(hf,'position'))
        warning('Max figure size reached --> reduced size')
    end
    movegui(hf,'center'); set(hf,'visible','on');
    ha = fig_createAxes(hf,[noAXI,1],ax.dx, ax.dy,'pixel');
    set(ha,'unit','normalized')
    ht = NaN(size(ha));
    
    %axes loop
    stat = struct();
    for axi = 1:noAXI
        [~,~,labTYP,labUNI] = opt.axi.data{axi,:};
        
        %MEAN & ERROR
        Dat  = squeeze(mDATA.Data(:,axi,:));
        N    = sum(all(~isnan(Dat),1));
        mDat = nanmean(Dat,2);
        err  = nanstd(Dat,[],2)/sqrt(N); %SEM
        
        %STATISTIC
        M    = Dat'; 
        indS = 0;
        %stat 1
        t0 = t/opt.plot.tFactor; %[s]
        TC = {... time points to compare in [s]
            -0.5*60,   0*60;...
            -0.5*60, 0.5*60;...
            -1.0*60, 1.0*60;...
            -1.5*60, 1.5*60;...
            -2.5*60, 2.5*60;...
            -1.7*60, 1.7*60;...
            };
        for k = 1:size(TC,1)
            indS = indS+1;
            [t1,t2] = TC{k,:};
            [~,ind1] = min(abs(t0-t1));
            [~,ind2] = min(abs(t0-t2));
            M1 = M(:,ind1);
            M2 = M(:,ind2);
            pval = PTTEST(M1,M2);
            str  = sprintf('(%+5.1f,%+5.1f), p = %.3f',...
                t1*opt.plot.tFactor,t2*opt.plot.tFactor,pval);
            stat(axi).res(indS).label = str;
            stat(axi).res(indS).t     = t([ind1,ind2]);
            stat(axi).res(indS).pval  = pval;
            stat(axi).res(indS).larger = ...
                nanmean(M1)>nanmean(M2);
        end
        %stat 2
        indS = indS+1;
        m = nanmean(M,2);
        stat(axi).res(indS).label = ...
            'Each time point to mean';
        stat(axi).res(indS).t      = t;
        stat(axi).res(indS).pval   = NaN(size(t));
        stat(axi).res(indS).larger = false(size(t));
        for k = 1:numel(t)
            stat(axi).res(indS).pval(k)   = PTTEST(M(:,k),m);
            stat(axi).res(indS).larger(k) = ...
                nanmean(M(:,k))>nanmean(m);
        end
        stat(axi).res(indS).larger = ...
            nanmean(M,1)>nanmean(m);

        %PLOT DATA AXIS
        set(hf,'CurrentAxes',ha(axi)); %first on top
        Y    = [mDat+err; mDat(indR)-err(indR)];
        mima = [min(Y(:)),max(Y(:))];
        yl   = mima+diff(mima)*[-0.1,0.2];
                    set(gca,'ylim',yl)
        col  = ones(1,3)*0.2;
        %error + mean + event line + text
        fill(X,Y,Z,'edgecolor','none',...
            'facecolor',col,'facealpha',0.2); hold on
        plot(t,mDat,'color',col,'linewidth',2.5);
        xline(0,'-');
        %text/settings
        ht(axi) = text(min(xl),max(yl),sprintf(' N = %i',N'),...
            'horizontalalignment','left',...
            'verticalalignment','top');
        if isempty(labUNI)
            ylabel(labTYP)
        else
            ylabel({labTYP,labUNI})
        end
        set(gca,'xlim',xl,'ylim',yl,'box','off')
        
        %PLOT STATISTIK
        S   = stat(axi).res;
        statTXT = {};
        for k = 1:numel(S)
            lab = S(k).label;
            tt  = S(k).t(:);
            pp  = S(k).pval(:);
            ll  = S(k).larger(:);
            if pp(1)<0.05 %for single p-val only
                str = sprintf('{\\color{red}p = %.2g}',pp);
            else
                str = sprintf('{\\color{black}p = %.2g}',pp);
            end
            if numel(pp)==1
                if true
                    if pp<0.05
                        statTXT{end+1} = sprintf(...
                            '{\\color{red}%s}',lab);
                    else
                        statTXT{end+1} = lab;
                    end
                else %old
                    tt(tt==0) = []; %to plot
                    for q = 1:numel(tt)
                        h = xline(tt(q),':','color','k');
                        h.LineWidth  = 1;
                        h.FontWeight = 'bold';
                        h.FontSize   = 8;
                        h.LabelOrientation = 'horizontal';
                        h.LabelVerticalAlignment   = 'bottom';
                        h.LabelHorizontalAlignment = 'center';
                        if q==1
                            h.Label = str;
                        end
                    end
                end
            else
                h   = NaN(4,1);
                leg = cell(4,1);
                tmp = {...
                    ll&( pp<0.10),[1.0 0.6 0.6],'< 0.1, larger';...
                    ll&( pp<0.05),[1.0 0.0 0.0],'< 0.5, larger';...
                    ~ll&(pp<0.10),[0.6 0.6 1.0],'< 0.1, smaller';...
                    ~ll&(pp<0.05),[0.0 0.0 1.0],'< 0.5, smaller';...
                    };
                for q = 1:4
                    [ind,col,leg{q}] = tmp{q,:};
                    y      = NaN(size(pp));
                    y(ind) = max(ylim)-0.05*diff(ylim);
                    h(q)   = plot(tt,y,'.','color',col,...
                        'MarkerSize',8);
                end
                if false %axi==1
                    pos = get(gca,'position');
                    legend(h,leg,'location',...
                        'northeastoutside');
                    set(gca,'position',pos);
                end
                if axi==noAXI
                    pos = get(gca,'position');
                    hl = legend(h,leg,...
                        'orientation','horizontal',...
                        'location','bestoutside');
                    set(gca,'position',pos);
                    posL = get(hl,'position');
                    posL(1) = pos(1);
                    set(hl,'position',posL)
                end
            end
        end %stat loop
        if numel(statTXT)>0
            text(max(xlim),min(ylim),...
                [' ',strjoin(statTXT,'\n ')],...
                'horizontalalignment','left',...
                'verticalalignment','bottom')
        end  
    end %axis loop
    
    %TEXT
    %title
    set(hf,'CurrentAxes',ha(1)); %first on top
    title({sprintf('%s Transition',Events),...
        'Mean \pm SEM'})
    %xlabel/ylabel
    set(hf,'CurrentAxes',ha(end));
    xlabel(sprintf('Time [%s], dt = %g s',opt.plot.tUnit,dt))
    fig_superLabel(ha,'ylabel','Temperature [C^o]');
    %legend
    set(hf,'CurrentAxes',ha(end));
    pos = get(gca,'position');
    hl = legend(h,leg,...
        'orientation','horizontal',...
        'location','bestoutside');
    set(gca,'position',pos);
    posL = get(hl,'position');
    posL(1) = pos(1);
    set(hl,'position',posL)
    
    %SETTINGS
    linkaxes(ha,'x')
    set(ha(1:end-1),'xticklabel',[])
    zoom yon
    
    %SAVE
    if opt.fig.save
        %save figure
        for k = 1:numel(opt.fig.saveFuns)
            fun = opt.fig.saveFuns{k};
            fun(hf,sname);
        end
        fprintf('Saved: %s\n',sname)
    end 
end


%% FIGURE, INDIVIDUAL DATA
if opt.plot.indData
    %clc; close all
    t    = mDATA.t; 
    %file loop
    for fil = 1:noFIL
        %%
        %close all;
        Data = DATA(fil);
        if all(isnan(Data.Data(:)))
            fprintf(2,'%s, all excluded!\n',Data.id)
            continue
        end
        if ~isequal(Data.t,t)
            error('Fatale error, grrhhh')
        end
        noEVT  = size(Data.Data,3);
        colors = lines(noEVT); %per event
        
        %FIGURE
        hf = figure('position',Fig.position);
        movegui(hf,'center'); drawnow
        ha = fig_createAxes(hf,[noAXI,1],opt.axi.dx,opt.axi.dy,'pixel');
        set(ha,'unit','normalized')
        %axes loop
        for axi = 1:noAXI
            set(hf,'CurrentAxes',ha(axi))
            [~,~,labTYP,labUNI] = opt.axi.data{axi,:};
            data = squeeze(Data.Data(:,axi,:));
            
            %PLOT
            hp     =  NaN(noEVT+1,1); 
            legSTR = cell(noEVT+1,1);
            %event loop
            for evt = 1:noEVT 
                legSTR{evt} = datestr(Data.event(evt),'dd-mmm HH:MM:SS');
                hp(evt)     = plot(t,data(:,evt),'color',colors(evt,:),...
                    'linewidth',1.5); hold on
            end
            %mean data
            legSTR{noEVT+1} = 'Mean';
            N = noEVT-sum(all(isnan(data),1));
            if N>1
                hp(noEVT+1) = plot(t,nanmean(data,2),'k','linewidth',2.5);
            else
                hp(noEVT+1) = plot(t,NaN(size(t)),'k','linewidth',2.5);
            end
            %event line
            xline(0,'-',{Events});
            
            %SETTINGS I
            yl = [min(data(:)), max(data(:))];
            if numel(yl)==2 && all(isfinite(yl))
                set(gca,'ylim',yl+[-0.1 0.15]*diff(yl))
            else
                set(gca,'ylim',[-1,1])
            end
            
            %TEXT
            if axi==1
                title({Data.id,sprintf('Transition %s (N = %i)',...
                    Events,noEVT)})
                posA = get(gca,'position');
                legend(hp,legSTR,'location','northeastoutside')
                set(gca,'position',posA)
            end
            ylabel({labTYP,labUNI})
            if axi==noAXI
                xlabel(sprintf('Time [%s]',opt.plot.tUnit))
            end
            text(max(opt.plot.margin),max(ylim),sprintf('N = %i ',N),...
                'horizontalalignment','right','verticalalignment','top')
        end %axis loop
        
        %SETTINGS II
        linkaxes(ha,'x')
        set(ha,'xlim',opt.plot.margin,'xgrid','on','ygrid','on')
        set(ha(1:end-1),'xticklabel',[])
        zoom yon
        
        %SAVE
        sFile = sprintf('Trans%s_id_%s',Events,id,opt.plot.fun);
        if opt.fig.save
            sname = fullfile(sPathFIG,sFile);
            for k = 1:numel(opt.fig.saveFuns)
                fun = opt.fig.saveFuns{k};
                fun(hf,sname);
            end
            fprintf('Saved: %s\n',sname)
        end
    end %figure loop
end


%% PRINT OUT
if ~opt.fig.save
    fprintf('[\bFigures NOT saved!]\b\n')
end