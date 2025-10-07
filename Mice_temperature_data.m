addpath(genpath('Z:\OptoLab_v4.1\function'))
%addpath('functions')
clc; close all; clear

%PARAMETERS
%----------
%FILES

tmp   = dir('*.xlsx');
files = fullfile({tmp.folder},{tmp.name})';
files = selectionList(files);
%save path
sPath = '';


%STAGES
Stages = {... {number, label}
    1, 'Wake';...
    2, 'NREM';...
    3, 'REM' ;...
    4, 'Cataplexy';...
    };

%function to append to transition data
funs = {... {label, function}
    'absolute',@(t,T)T; ... 
    'change'  ,@(t,T)T-T(t==0);... zero at event
    };
[funLab,fun1] = funs{2,:};
%transitions cataplexy
k=1;
Trans(k).plot  = true; %find & plot transitions
Trans(k).label = 'Cataplexy';
Trans(k).trans = {{'Wake','NREM','REM'},{'Cataplexy'}}; 
Trans(k).win   = [-20,10]*60; %[s]
Trans(k).fun0  = @(t,T)T;    
Trans(k).fun1  = fun1;
Trans(k).lim.dur = 0;                %[s], min duration
Trans(k).lim.gap = abs(Trans(k).win(1)); %[s], time lock gap
Trans(k).lim.tim = NaN;              %[h], selected clock hours

Trans(k).stat.plot = true;

%PLOT OPTIONS
dt = 1; %time step in files
VarNames = {...{label,cols1, cols2} --> mean(cols1) - mean(cols2)
    'DeltaT Back '       , {'SkinT'}      , {};...
    'DeltaT Neck', {'body_center'}, {};...
    'DeltaT Right ear'  , {'ear_right'}  , {};...
    'DeltaT Tail'       , {'tail'}       , {};...
    'DeltaT Distal-Proximal' , {'tail','ear_right'}, {'SkinT','body_center'};...
    'DeltaT Core'       , {'CoreT'}      , {};...
    'DeltaT Ta'    , {'AmbientT'}   , {};...
    };
%time units to plot
plo.meanOnly  = true;
plo.time.unit = 'min';
plo.time.fac  = 1/60; 
%general options

opt.linkaxes = false; 

%smoothing
opt.smooth.do = true;
opt.smooth.funs = {... cell, functions applied in this order !!!
    @(x,y)interp1(x(~isnan(y)),y(~isnan(y)),x,'linear',nanmean(y));...
    @(x,y)passband_fourier(y,[0,0.05],1/dt);...
    @(x,y)movmean(y,200,'omitnan');...
    };
%save options
opt.save.do = true;
if opt.smooth.do
    opt.save.sFile = @(transition,id)... 
        sprintf('trans%s_smoothed_%s_%s',... 
        transition,id,funLab);
else
    opt.save.sFile = @(transition,id)... 
        sprintf('trans%s_%s_%s',... 
        transition,id,funLab);      
end
opt.save.funs = {... sname id fullfile(sPath,sFile)
    @(h,sname)print(sname,'-dpng','-r300');...
    @(h,name)saveas(h,name);...
    };

%MAIN SCRIPT
%-------------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))
%check
for k = 1:numel(Trans)
    trans = Trans(k);
    if all(isfinite(trans.lim.tim)) && trans.lim.tim(1)==trans.lim.tim(2)
        error('Wrong time limit for sleep transition')
    end
end
if ~exist(sPath,'dir')
    mkdir(sPath)
end


%for to prevent re-reading
reread = true;
if exist('RR','var') && all(ismember(fieldnames(RR),{'name','files'}))...
        && isequal(RR.files,files) && strcmp(RR.name,scriptName)
    reread = false;
else
    RR.files = files;
    RR.name  = scriptName;
end


%% READ DATA
fprintf('READ DATA\n'); tic
noFIL  = numel(files);
nnFIL  = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
if reread
    for fil = 1:noFIL
        file = files{fil};
        [~,rFile,rExt] = fileparts(file);
        id = regexp(rFile,'_','split');
        id = strjoin(id(2:end),'_');
        fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,id)
        
        %READ
        T = readtable(file,'ReadVariableNames',true);
        
        for k = 1:numel(T.Properties.VariableNames)
            name = T.Properties.VariableNames{k};
            if all(T.(name)==0)
                T.(name) = NaN(size((T.(name))));
            end
        end
        %append
        DATA(fil).id      = id;
        DATA(fil).T       = T;
        DATA(fil).hastrans = false(numel(Trans),1);
        for k = 1:numel(Trans)
            trans = Trans(k);
            ind1  = ismember(lower(Stages(:,2)),lower(trans.trans{1}));
            ind2  = ismember(lower(Stages(:,2)),lower(trans.trans{2}));
            stages1 = cell2mat(Stages(ind1,1));
            stages2 = cell2mat(Stages(ind2,1));
            DATA(fil).hastrans(k) = ...
                any(ismember(T.SleepStage,stages1)) && ...
                any(ismember(T.SleepStage,stages2));
        end
        %check
        tmp = etime(datevec(T.datenum(2)),datevec(T.datenum(1)));
        if abs(dt-tmp)/dt*100>10^-10 
            error('Time Step Error, check dt = %g',tmp)
        end
    end
    toc
    tmp = T.Properties.VariableNames;
    fprintf('Variables Names (N = %i):\n',numel(tmp))
    fprintf(' ''%s''\n',strjoin(tmp,''', '''))
else
    fprintf('[\bWas already done!]\b\n')
end


%% TRANSTIONS
%clc; close all %%%
HasTrans = [DATA.hastrans];
for kk = 1:numel(Trans)
    DataT    = DATA(HasTrans(kk,:)); %only the one having a transition
    trans    = Trans(kk);
    labelTRA = trans.label;
    fprintf('\n%s TRANSITIONS\n',upper(labelTRA));
    if numel(DataT)==0
        fprintf(2,'No data with transition\n')
        continue
    end
    clear TAB
    
    %FIND TRANSITION DATA
    % trans.plot = true;
    if trans.plot
        ind0   = floor(trans.win(1)/dt):ceil(trans.win(2)/dt);
        t0     = ind0*dt; %transition time vector [s] 
        noDAT  = numel(DataT);
        nnDAT  = numel(num2str(noDAT));
        indent = blanks(2*nnDAT+2);
        for dat = 1:noDAT
            data = DataT(dat);
            fprintf('%*i/%i: %s\n',nnDAT,dat,noDAT,data.id)
            
            %FIND TRANSITIONS
            %transition stages
            indB = ismember(lower(Stages(:,2)),lower(trans.trans{1}));
            indA = ismember(lower(Stages(:,2)),lower(trans.trans{2}));
            stageB = cell2mat(Stages(indB,1));
            stageA = cell2mat(Stages(indA,1));
            %transitions start/end
            stages = data.T.SleepStage;
            ind    = ismember(stages,stageA);
            indSTA = find([ind;0]==1 & [0;ind]==0);
            indEND = find([ind;0]==0 & [0;ind]==1)-1;
            fprintf('%s Transitions found: N = %i\n',indent,numel(indSTA));
            
            %REMOVE TRANSITIONS
            rem0 = struct('ind',false(size(indSTA)),'txt','');
            rem  = struct('tot',rem0);
            strs = fieldnames(trans.lim);
            for k = 1:numel(strs)
                str = strs{k};
                lim = trans.lim.(str);
                r   = rem0; %init
                if ~isempty(lim) && all(isfinite(lim))
                    switch str
                        case 'dur'
                            dur   = (indEND-indSTA+1)*dt;
                            r.ind = dur<lim;
                            r.txt = sprintf('duration < %g s',lim);
                        case 'gap'
                            gap   = [indSTA(1)-1;... 
                                (indSTA(2:end)-indEND(1:end-1)-1)]*dt;
                            r.ind = gap<lim;
                            r.txt = sprintf('within lock gap %g s',lim);
                        case 'tim'
                            t = datevec(data.T.datenum(indSTA)) * ...
                                [0 0 0 1 1/60 1/3600]'; %clock hour
                            if lim(1)<lim(2)
                                r.ind = ~(t>=lim(1) & t<=lim(2));
                            else
                                r.ind = t>lim(2) & t<lim(1);
                            end
                          
                            if any(ismember(lim,[0,24]))
                                r.ind(ismember(t,[0,24])) = false;
                            end
                            r.txt = sprintf(...
                                'beyond clock hours %.2f-.2f%g',lim);
                        case {'durBS','durAS'}
                            d0  = lim(1);
                            p0  = lim(2);
                            txt = sprintf('lim %s < %i%%',str,p0);
                            if strcmp(str,'durBS')
                                ind1 = t0>=d0 & t0<0;
                                indC = ismember(stages,stageB);
                            else
                                ind1 = t0>=0 & t0<lim(1);
                                indC = ismember(stages,stageA);
                            end
                            for q = 1:numel(indSTA)
                                ind2 = ind1+indSTA(k);
                                ind2(ind2<1 | ind2>numel(stages))=[];
                                p      = sum(indC(ind2))/numel(ind2)*100;
                                indR(q) = p<p0;
                            end
                            r.txt = sprintf('lim %s < %i%%',str,p0);
                        otherwise
                            error(sprintf(...
                                'trans.lim.%s not implemented',str))
                    end
                end
                %append
                rem.tot.ind = rem.tot.ind | r.ind;
                rem.(str)   = r;
            end
            %remove
            ind = rem.tot.ind;
            if sum(ind)>0
                fprintf('%s   Removed, N = %i\n',indent,sum(ind));
                indSTA(ind) = [];
                indEND(ind) = [];
                for k = 1:numel(strs)
                    str = strs{k};
                    if sum(rem.(str).ind)>0
                        fprintf('%s     %i: %s\n',indent,...
                            sum(rem.(str).ind),rem.(str).txt);
                    end
                end
                fprintf('%s   Remaining, N = %i\n',indent,numel(indSTA));
            end
            
            %SELECTED DATA / SMOOTH / TRANSITION DATA
            varNames = VarNames(:,1);
            noVAR = numel(varNames );
            noTRA = numel(indSTA);
            noSAM = size(T,1);
            %data of selected variables
            TT = NaN(noSAM,noVAR);
            for v = 1:noVAR
                [~,cols1,cols2] = VarNames{v,:};
                n1 = numel(cols1);  n2 = numel(cols2);
                d1 = NaN(noSAM,n1); d2 = NaN(noSAM,n2);
                for c = 1:n1
                    d1(:,c) = data.T.(cols1{c});
                end
                for c = 1:n2
                    d2(:,c) = data.T.(cols2{c});
                end
                if n2==0
                    TT(:,v) = nanmean(d1,2);
                else %e.g. Distal-Proximal
                    TT(:,v) = nanmean(d1,2)-nanmean(d2,2);
                end
            end
            %smooth data
            if opt.smooth.do
                for v = 1:noVAR
                    x = data.T.datenum;
                    y = TT(:,v);
                    if all(isnan(y)); continue; end
                    for f = 1:numel(opt.smooth.funs)
                        fun = opt.smooth.funs{f};
                        y = fun(x,y);
                    end
                    if all(isnan(y)); continue; end
                    TT(:,v) = y;
                end
            end
            %transition matrix
            M = NaN(noTRA,numel(ind0),noVAR); 
            for tra = 1:noTRA
                ind  = ind0+indSTA(tra);
                indM = ind>0 & ind<=numel(stages);
                ind(~indM) = [];
                for v = 1:noVAR
                    M(tra,indM,v) = trans.fun1(t0,TT(ind,v));
                end
            end
            
            %APPEND DATA
            DataT(dat).trans.info = 'Transition: trans x time x varNames';
            DataT(dat).trans.noTRA      = noTRA;
            DataT(dat).trans.datenum    = data.T.datenum(indSTA);
            DataT(dat).trans.varNames   = varNames;
            DataT(dat).trans.time       = t0;
            DataT(dat).trans.Transition = M;
            
            %MEAN DATA ACROSS TRANSITIONS
            if dat==1
                meanM = NaN(noDAT,numel(ind0),noVAR);
            end
            meanM(dat,:,:) = nanmean(M,1);
            if dat==noDAT
                data.id = 'Mean';
                data.T  = [];
                data.trans.info = 'Transition: trans x time x varNames';
                data.trans.noTRA      = noDAT;
                data.trans.datenum    = [];
                data.trans.varNames   = varNames;
                data.trans.time       = t0;
                data.trans.Transition = meanM;
                DataT(dat+1) = data;
            end
        end
        
        %PLOT TRANSITIONS
        fprintf('PLOT\n')
        if plo.meanOnly
            tmp = numel(DataT);
        else
            tmp = 1:numel(DataT);
        end
        %INIT TABLE
        clear TAB
        for dat = tmp
            %get data
            data   = DataT(dat);
            tStart = data.trans.datenum;
            t      = data.trans.time*plo.time.fac; 
            MM     = data.trans.Transition;
            NN     = squeeze(sum(~all(isnan(MM),2),1));
            labelsAXI = data.trans.varNames;
            [noTRA,noSAM,noAXI] = size(MM); 
            %check
            assert(noTRA==data.trans.noTRA   ,'ggrrrhhhh');
            assert(noSAM==numel(t)           ,'ggrrrhhhh');
            assert(noAXI==noVAR              ,'ggrrrhhhh');
            assert(noAXI==numel(labelsAXI)   ,'ggrrrhhhh');
            assert(isequal(t,t0*plo.time.fac),'ggrrrhhhh');
            
            %INIT TABLE
            if ~exist('TAB','VAR')
                labelsANI = {DataT.id};
                labelsANI(strcmpi(labelsANI,'mean')) = [];
                TAB = cell(numel(labelsANI)+1, numel(t)+1, noAXI);
            end
            
            %FIGURE AXES
            ax.hig = 110; ax.wid = 375;
            ax.dx  = [75,0,140]; ax.dy  = [70,7,40];
            posF = [0,0,ax.wid+sum([65,0,120]),...
                noAXI*ax.hig+[1,noAXI-1,1]*ax.dy'];
            hf  = figure('position',posF,'visible','off');
            movegui(hf,'center'); set(hf,'visible','on');
            ha = fig_createAxes(hf,[noAXI,1],ax.dx, ax.dy,'pixel');
            set(ha,'unit','normalized')
            %init
            ht   = NaN(noAXI,1);
            MIMA = [ones(noAXI,1)*-100,ones(noAXI,1)*100]; %init
            pval = NaN(1,noAXI);
            stat = struct();
            
            %axis loop
            for axi = 1:noAXI
                set(hf,'CurrentAxes',ha(axi));
                labelAXI = labelsAXI{axi};
                M = MM(:,:,axi);
                N = NN(axi);
                MIMA(axi,1) = min([MIMA(axi,1);M(:)]);
                MIMA(axi,2) = max([MIMA(axi,2);M(:)]);
                    
                %PLOT
                if strcmpi(data.id,'Mean')
                    %APPEND TO TABLE
                    tmp = cellfun(@(x)all(isempty(x(:))),TAB(:,:,axi));
                    assert(all(tmp(:)),'grrhhh');
                    TAB(1,1,axi)         = {'t'};
                    TAB(1,2:end,axi)     = num2cell(t);
                    TAB(2:end,1,axi)     = labelsANI(:);
                    TAB(2:end,2:end,axi) = num2cell(M);
                   
                    mea = nanmean(M,1);
                    err = nanstd(M,[],1)/sqrt(N);
                    %plot
                    xline(0,'k'); hold on
                    col  = ones(1,3)*0.2;
                    indR = numel(t):-1:1; 
                    X = [t,t(indR)]; Z = zeros(size(X));
                    Y = [mea-err,mea(indR)+err(indR)];
                    h = fill(X,Y,Z,'edgecolor','none',...
                        'facecolor',col,'FaceAlpha',0.3);
                    plot(t,mea,'color',col,'linewidth',2)
                    if axi==1
                        title({sprintf('%s Transition',labelTRA),...
                            sprintf('Mean \\pm SEM')})
                    end
                    %setting
                    mima = [min(Y(:)),max(Y(:))];
                    yl   = mima+diff(mima)*[-0.1,0.2];
                    set(gca,'ylim',yl)
                    ht(axi) = text(min(t),0,sprintf('  N = %i',N),...
                        'HorizontalAlignment','left',...
                        'VerticalAlignment'  ,'top');
                    
                    %STATISTIC
                    indS = 0;
                    %stat 1
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
                        [~,pval] = ttest(M1,M2);
                        %pval = PTTEST(M1,M2);
                        str  = sprintf('(%+5.1f,%+5.1f), p = %.3f',...
                            t1*plo.time.fac,t2*plo.time.fac,pval);
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
                        [~,tmp] = ttest(M(:,k),m);
                        stat(axi).res(indS).pval(k)   = tmp;
                        stat(axi).res(indS).larger(k) = ...
                            nanmean(M(:,k))>nanmean(m);
                    end
                    stat(axi).res(indS).larger = ...
                        nanmean(M,1)>nanmean(m);
                else
                    xline(0,'k'); hold on
                    hp  = NaN(noTRA+1,1);
                    leg = cell(noTRA+1,1);
                    for tra = 1:noTRA %transition loop
                        hp(tra)  = plot(t,M(tra,:));hold on
                        leg{tra} = datestr(tStart(tra),'HH:MM:SS');
                    end
                    mea = nanmean(M,1);
                    hp(end)  = plot(t,mea,'k','linewidth',2);
                    leg{end} = 'Mean';
                    %title/legend
                    if axi==1
                        title({sprintf('%s Transition',labelTRA),...
                            strrep(data.id,'_','\_')})
                        pos = get(gca,'position');
                        legend(hp,leg,'location','northeastoutside')
                        set(gca,'position',pos)
                    end
                    %setting
                    mima = [min(M(:)),max(M(:))];
                    if ~any(isfinite(mima)) || mima(1)==mima(2)
                        mima = [-1,1];
                    end
                    if all(mea(t==0)==0)
                        set(gca,'ylim',max(abs(mima))*[-1.1,1.1])
                    else
                        set(gca,'ylim',mima+diff(mima)*[-0.1,0.1])
                    end
                    %text
                    ht(axi) = text(max(t),0,sprintf('N = %i ',N),...
                        'HorizontalAlignment','right',...
                        'VerticalAlignment'  ,'top');
                end
                if axi==noVAR
                    xlabel(sprintf('Time [%s], dt = %g s',...
                        plo.time.unit,dt))
                end
                ylabel(strrep(labelAXI,'_','\_'))
            end %axi loop
            hy = fig_superLabel(ha,'ylabel','Temperature [C^o]');
            
            %SETTINGS
            set(ha,'xlim',t([1,end]),'box','off')
            set(ha(1:end-1),'xticklabel',[])
            if opt.linkaxes
                linkaxes(ha,'y')
                mima = [min(MIMA(:,1)),max(MIMA(:,2))];
                if all(Y(X==0)==0)
                    set(ha,'ylim',max(abs(mima))*[-1.1,1.1])
                else
                    set(ha,'ylim',mima+diff(mima)*[-0.1,0.1])
                end
            end
            linkaxes(ha,'x')
            drawnow
            for k = 1:noAXI
                pos = get(ht(k),'Position');
                pos(2) = max(get(ha(k),'ylim'));
                set(ht(k),'position',pos)
            end
            zoom yon
            
            %PLOT STATISTIC
            if dat==numel(DataT) %mean data only
                for axi = 1:noAXI
                    set(hf,'CurrentAxes',ha(axi))
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
                            else 
                                tt(tt==0) = []; 
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
                            if false 
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
                %settings 
                set(hf,'CurrentAxes',ha(axi));
                set(ha,'box','off')
                linkaxes(ha,'x')
                
            end %for mean plot
            
            %SAVE
            if opt.save.do
                sFile = opt.save.sFile(labelTRA,...
                    strrep(data.id,'_','-'));
                sname = fullfile(sPath,sFile);
                for k = 1:numel(opt.save.funs)
                    fun = opt.save.funs{k};
                    fun(hf,sname);
                    if k==numel(opt.save.funs)
                        fprintf('Saved: %s\n',sFile)
                    end
                end
            end
        end %figures
    end %trans.plot
    
    if opt.save.do
        sFile  = opt.save.sFile(labelTRA,strrep('data','_','-'));
        sname  = fullfile(sPath,[sFile,'.xlsx']);
        sheets = labelsAXI;
        if exist(sname','file')==2
            delete(sname)
        end
        warning off
        for axi = 1:noAXI
            xlswrite(sname,TAB(:,:,axi),sheets{axi});
        end
        warning on
        try xls_deleteSheets(sname,sheets,-1); catch; end
        try xls_cellFit(sname); catch; end
        
    end
    
    
end %loop transition types
if ~opt.save.do
    fprintf(2,'\nFIGURES NOT SAVED !!!\n')
end