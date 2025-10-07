clc; clear; close all;
addpath(genpath('Z:\OptoLab_v4.1\function'))

rFiles = ['*Ca.mat'];

opt.savePath = '';
opt.saveDo   = true;
opt.saveFuns = {... input for all: @(figure handle,savename)
    @(h,name)print(h,name,'-r300','-dpng');...
    ...@(h,name)print(h,name,'-r300','-dpdf');...
    };

%ANALYSIS OPTIONS
%offset
%  if true, set offset (appends offset to Ca-file)
%  else reads from Ca-file or uses zero (if not exist in Ca-file)
opt.offset = true;
%function arguments
opt.detrend.win = 2500; %[s] (for detrend SamplePoints)
opt.detrend.degree = 2; %polynomal degree
opt.movmean.win = 25;  %movmean time [s] (uses odd bins of next larger)
opt.cutAfterStage4 = false; %cuts analsys after last stage 4

%stages in hypnogram
Stages = {... {number,label}
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    4,'Cataplexy';...
    };
%stage transitions
Transitions = {... {stage from, stage to}
    'Wake','NREM';...
    'NREM','Wake';...
    'REM','Wake';...
    'NREM','REM';...
    'Wake','Cataplexy';...
    'Cataplexy', 'Wake';...
    };
margin  = 90; %[s], plots transtions +- margin (figure 2)
exclNaN = true; %exclude transitions with margin beyond data range


%PLOT PROPERTIES
props.figure = {'visible','on'}; 
props.axi.hyp = {'box','off'}; 
props.axi.cal = {'box','off','ylim',[0,1]}; 
props.axi.sta = {'box','off'}; 
props.axi.tra = {'box','off'}; 
props.equal.ylimTRA = true; 
%plot
props.plo.hyp = {};
props.plo.cal = {};
props.plo.sta = {}; 
props.plo.tra = {'linewidth',2};
%error plotting
errorTrans  = 'SEM'; 

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))
if ~opt.saveDo
    fprintf('[\bNOTE: data will not be saved (variable opt.saveDo)]\b\n')
end

%FILE-LIST
if iscell(rFiles)
    ind = cellfun(@(x)exist(x,'file')~=2,rFiles);
    if any(ind)
        fprintf('[\bFiles removed, %i of %i (do NOT exist)!]\b\n',...
            sum(ind),numel(ind))
        fprintf(' - %s\n',rFiles{ind})
        rFiles(ind) = [];
    end
    if isempty(rFiles)
        fprintf(2,'No File Left!\n')
        return
    end
    rFiles = selectionList(rFiles);
elseif ischar(rFiles)
    tmp = dir(rFiles);
    tmp([tmp(:).isdir]) = [];
    tmp = fullfile({tmp.folder},{tmp.name});
    if isempty(tmp)
        fprintf(2,'No File Found For: %s\n',rFiles)
        return
    end
    rFiles = selectionList(tmp(:),'Ca*_DFF.mat');
    if isempty(rFiles)
        fprintf(2,'No File Selected\n')
        return
    end
else
    error('Class ''%s'' not supported for variable rFiles',class(rFiles))
end
noFIL = numel(rFiles);

%INIT
if ~iscell(opt.saveFuns)
    opt.saveFuns = {opt.saveFuns};
end
%number of ...
noSTA = size(Stages,1);
noTRA = size(Transitions,1);

%PLOT PROPERTIES
tmp = cell2mat(cellfun(@(x)x(:),Stages(:,1),'UniformOutput',false));
ticks = min(tmp):max(tmp);
label(size(ticks)) = {''};
for k = 1:noSTA
    [stageNUM,stageLAB] = Stages{k,:};
    label{ismember(ticks,stageNUM)} = stageLAB;
end
props.axi.hyp = ['ytick',ticks,'yticklabel',{label},props.axi.hyp];
%check (props must not be empty)
tmp = fieldnames(props.axi);
for k = 1:numel(tmp)
    if isempty(props.axi.(tmp{k}))
        props.axi.(tmp{k}) = {'visible','on'};
    end
end
tmp = fieldnames(props.plo);
for k = 1:numel(tmp)
    if isempty(props.plo.(tmp{k}))
        props.plo.(tmp{k}) = {'visible','on'};
    end
end

%FILE LOOP
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    clear res
    [rPath,rFile,rExt] = fileparts(rFiles{fil});
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath)
    %files
    fnames.cal = fullfile(rPath,[rFile,'.mat']);
    fnames.hyp = fullfile(rPath,'Hypnogram.mat');
    if ~exist(fnames.hyp,'file')
        [~,rFile,rExt] = fileparts(fnames.hyp);
        fprintf('%s [\b%s NOT found]\b\n',indent,[rFile,rExt])
        continue
    end
    tmp = {...
        sprintf('Created by %s.m, %s',scriptName,date);'';...
        'Resulting offset is already cut from results';...
        };
    res.info.info = char(tmp);
    
    %% LOAD CA-DATA
    fprintf('%s Load Data\n',indent)
    data   = load(fnames.cal);
    fs     = data.SampRate;
    signal = data.demodulated_signal(:);
    noSAM  = numel(signal);
    fprintf('%s   Ca-signal, fs = %g Hz\n',indent,fs)
    %id for title
    id = strrep(rFile,'_','\_'); %init
    if isfield(data,'info') && isfield(data.info,'file')
        [~,id] = fileparts(data.info.file); %from recorded filename
        id = strrep(id,'_','\_');
    end
    %get offset (cut later!!!)
    if isfield(data,'offset')
        offset = data.offset;
    else
        offset = []; %init
    end
    if opt.offset
        off = fun_CaOffset(signal,'fs',fs,'offset',offset,...
            'title',regexprep(fnames.cal,{'\','_'},{'\\\\','\\_'}));
        if ~isequal(off,offset)
            offset = off;
            save(fnames.cal,'offset','-append')
            fprintf('%s     Offset: %i (appended to calcium-file)\n',...
                indent,offset)
        else
            fprintf('%s     Offset: %i (old value kept)\n',...
                indent,offset)
        end
    else
        if isempty(offset)
            offset = 0;
            fprintf('%s     [\bOffset not set, uses 0]\b\n',indent)
        else
            fprintf('%s     Offset: %i (read from calcium-file)\n',...
                indent,offset)
        end
    end
    
    %% LOAD HYPNOGRAM
    data = load(fnames.hyp);
    hypnogram = data.Hypnogram(:);
    if isfield(data,'fs')
        fsHYP = data.fs;
        fprintf('%s   Hypnogram, fs = %g Hz\n',indent,fsHYP)
    else %estimate
        fsHYP = round(fs/noSAM*numel(hypnogram));
        fprintf('%s   Hypnogram, fs = %g Hz [\b(estimated!)]\b\n',...
            indent,fsHYP)
    end
    % resample hypnogram
    if fs>fsHYP %upsample
        fac = fs/fsHYP;
        if fac==round(fac)
            hypnogram = repmat(hypnogram',fac,1);
            hypnogram = hypnogram(:);
            fprintf('%s     Upsampled by factor %i\n',indent,fac)
        else
            tmp = round(resample(hypnogram,fs,fsHYP,0));
            tmp = tmp(:);
            %time shift correction (start at 1/fs, not at 1/fsHYP)
            shift = round((1/fsHYP-1/fs)*fs); % >= 0 !!!
            if shift>0
                tmp = [repmat(tmp(1),shift,1);tmp(1:end-shift)];
            end
            hypnogram = tmp(:);
            fprintf('%s     Upsampled to %i Hz\n',indent,fs)
        end
    elseif fs<fsHYP %downsample
        %PS: upsampling signal made everything very very slow
        %    donwsampling hypnogram is quite OK (error <= 1/fs/2)
        tmp = resample(hypnogram,fs,fsHYP,0);
        tmp = tmp(:);
        %time shift correction (start at 1/fs, not at 1/fsHYP)
        shift = round((1/fs-1/fsHYP)*fs); % >= 0 !!!
        if shift>0
            tmp = [tmp(shift+1:end);repmat(tmp(end),shift,1)];
        end
        if false %test plot
            figure
            t1 = (1:numel(hypnogram))/fsHYP;
            t2 = (1:numel(tmp))/fs;
            plot(t1,hypnogram,'b','linewidth',2); hold on
            plot(t2,tmp,':r','linewidth',2);
            mima = [min(hypnogram),max(hypnogram)];
            set(gca,'ylim',mima+0.1*[-1,1]*diff(mima),'ytick',...
                mima(1):mima(2),'xlim',[0,max([t1(end),t2(end)])])
            legend('hypnogram','downsampled')
            zoom xon
            return
        end
        hypnogram = tmp(:);
        fprintf('%s     Downsampled to %i Hz\n',indent,fs)
    end
    %same size signal & hypnogram 
    if numel(hypnogram)<noSAM
        hypnogram(end+1:noSAM) = NaN;
    elseif numel(hypnogram)>noSAM
        hypnogram(noSAM+1:end) = [];
    end
    
    %CUT OFFSET
    res.info.offset = offset;
    offSIG = []; offHYP = []; %init
    if offset>0
        ind = 1:offset;
        offSIG = signal(ind);
        offHYP = hypnogram(ind);
        signal(ind)    = [];
        hypnogram(ind) = [];
        noSAM = numel(signal);
    end
    
    %% DFF
    %trend
    fprintf('%s Calculate DFF\n',indent)
    %detrend signal
    win = opt.detrend.win*fs;
    trend = signal - detrend(signal,opt.detrend.degree,... 2 is fine
        1:win:noSAM,... SamplePoints, window 2500s seems ok
        'Continuous',false); %better with false but needs to be smoothed!
    trend = smooth(trend,win/2); %smooth because continuous false
    signalD = signal-trend;
    %smooth signal
    bins = round(opt.movmean.win*fs);
    if mod(bins,2)==0
        bins = bins+1;
    end
    signalM = movmean(signalD,bins); %smooth signal
    mi = prctile(signalD,1);
    ma = prctile(signalD,99);
    %dff (like dff)
    res.fs  = fs;
    res.dff = (signalM-mi)/(ma-mi);
    res.hypnogram = hypnogram;
    
    %% MEAN ACROSS STAGES
    res.stages.label = cell(1,noSTA);
    res.stages.mean  = NaN(1,noSTA);
    res.stages.std   = NaN(1,noSTA);
    res.stages.N     = NaN(1,noSTA);
    indCUT = true(size(hypnogram));
    if opt.cutAfterStage4
        indCUT(find(hypnogram==4,1,'last')+1:end) = false;
    end
    for sta = 1:noSTA
        [stage,label] = Stages{sta,:};
        ind = ismember(hypnogram,stage) & indCUT;
        res.stages.label{sta} = label;
        res.stages.mean(sta)  = mean(res.dff(ind));
        res.stages.std(sta)   = std(res.dff(ind));
        res.stages.N(sta)     = sum(ind);
    end
    %% MEAN ACROSS TRANSITIONS
    ind0 = -round(margin*fs):round(margin*fs);
    t = ind0'/fs;
    res.trans.t = t;
    res.trans.label = cell(1,noTRA);
    res.trans.Hypno = cell(1,noTRA);
    res.trans.Data  = cell(1,noTRA);
    res.trans.mean  = NaN(numel(t),noTRA);
    res.trans.std   = NaN(numel(t),noTRA);
    res.trans.N     = NaN(1,noTRA);
    for tra = 1:noTRA
        [stage1,stage2] = Transitions{tra,:};
        nums1 = Stages{strcmpi(Stages(:,2),stage1),1};
        nums2 = Stages{strcmpi(Stages(:,2),stage2),1};
        indTRA = find(...
            ismember(hypnogram(1:end-1),nums1) & ...
            ismember(hypnogram(2:end)  ,nums2));
        %matrix of transitions
        N = numel(indTRA);
        M = NaN(numel(ind0),N); %matrix DFF
        H = NaN(numel(ind0),N); %matrix corresponding stages
        for k = 1:N
            ind = ind0+indTRA(k); %%%CHECK AGAIN
            ind1 = find(ind>=1,1,'first');    %if beyond data range
            ind2 = find(ind<=noSAM,1,'last'); %if beyond data range
            M(ind1:ind2,k) = res.dff(ind(ind1):ind(ind2));
            H(ind1:ind2,k) = res.hypnogram(ind(ind1):ind(ind2));
        end
        %exclude transitions with margin beyond data range
        if exclNaN
            ind = any(isnan(M),1);
            M(:,ind) = [];
            H(:,ind) = [];
            N = size(M,2);
        end
        
        %res
        res.trans.label{tra}  = sprintf('%s - %s',stage1,stage2);
        res.trans.Hypno{tra}  = H; 
        res.trans.Data{tra}   = M;
        res.trans.mean(:,tra) = nanmean(M,2);
        res.trans.std(:,tra)  = nanstd(M,[],2);
        res.trans.N(tra)      = N;
    end
    
    %% PLOT INIT
    fprintf('%s Plots\n',indent)
    clear hf;
    
    %% FIGURE transitions
    % clc; close all; %%%
    switch 2
        case 1
            %figure & axes
            hf(3) = figure(props.figure{:}); figLab{3} = 'transitions_90';
            n1 = floor(sqrt(noTRA)); n2 = ceil(noTRA/n1);
            ha = NaN(n2,n1);
            for k = 1:noTRA
                ha(k) = subplot(n1,n2,k);
            end
            ha = ha';
            
            %PLOT
            t = res.trans.t(:);
            mima = [inf,-inf];
            for tra = 1:noTRA
                %plot lines
                set(hf(3),'CurrentAxes',ha(tra))
                xline(0); hold on;
                y = res.trans.mean(:,tra);
                hp = plot(t,y,props.plo.tra{:});
                %plot error
                ylab = sprintf('\\DeltaF/F + %s',errorTrans);
                switch lower(errorTrans)
                    case ''
                        err = [];
                        ylab= '\DeltaF/F';
                    case 'std'
                        err = res.trans.std(:,tra);
                    case 'sem'
                        err = res.trans.std(:,tra)/sqrt(res.trans.N(tra));
                    otherwise
                        error('errorTrans = ''%s'' is NOT implemented',errorTrans)
                end
                if ~isempty(err)
                    n  = numel(t);
                    y  = [y-err ; y(n:-1:1)+err(n:-1:1)];
                    h  = fill([t;t(n:-1:1)],y,get(hp,'color'),'facealpha',0.5,...
                        'EdgeColor','none');
                    uistack(h,'bottom')
                end
                mima = [min([mima(1);y(:)]),max([mima(2);y(:)])];
                %text
                if tra==1
                    title({id,sprintf('Transition (N = %i)',res.trans.N(tra)),...
                        res.trans.label{tra}})
                else
                    title({sprintf('Transition (N = %i)',res.trans.N(tra)),...
                        res.trans.label{tra}})
                end
                xlabel('Time [s]')
                ylabel(ylab)
                %settings
                set(gca,'xlim',[-margin,margin],props.axi.tra{:})
            end
            %special axes settings
            ind = ~isnan(ha);
            linkaxes(ha(ind),'x')
            if props.equal.ylimTRA
                set(ha(ind),'ylim',mima)
       
            end
        case 2
            len = 120; hig = 140;
            dx = [50,50,80];
            dy = [50,5,70];
            posF = [10,10,noTRA*len+[1,noTRA-1,1]*dx(:),3*hig+[2 1 1]*dy'];            
            hf(3) = figure(...
                'position',posF,...
                props.figure{:});
            movegui(gcf,'center'); drawnow
            figLab{3} = 'transitions_90'; 
            ha = NaN(3,noTRA);
            yy = dy+[dy(1)+hig 0 0];
            ha(1:2,:) = fig_createAxes(gcf,[2,noTRA],dx,...
                dy+[dy(1)+hig 0 0] ,'pixel');
            ha(3,:) = fig_createAxes(gcf,[1,noTRA],dx,...
                [dy(1) 0 sum(dy)+2*hig] ,'pixel');
            set(ha,'unit','normalized')
            
            %stage label, limit ,...
            clear stage
            stage.nums = min([Stages{:,1}]):max([Stages{:,1}]);
            for k = 1:numel(stage.nums)
                ind = cellfun(@(x)x==stage.nums(k),Stages(:,1));
                if all(~ind)
                    stage.labs{k} = '';
                else
                    stage.labs(k) = Stages(ind,2);
                end
            end
            
            %PLOT
            indR = numel(res.trans.t):-1:1; %reverse index
            t  = res.trans.t(:);
            te = [t;t(indR)]; %for error plots
            mima1 = [inf,-inf]; %mean
            mima2 = [inf,-inf]; %mean +- error
            mima3 = [inf,-inf]; %all transitions
            for tra = 1:noTRA
                %GET DATA
                H = res.trans.Hypno{tra}';
                Z = res.trans.Data{tra}';
                N = res.trans.N(tra);
                y = res.trans.mean(:,tra);
                if N~=size(Z,1)
                    error('uups')
                end
                %sort & normalize
                ZN = bsxfun(@times,Z,1./max(Z,[],2)); %normalization
                [~,ind] = sort(mean(ZN,2)); %sort index
                ZN = ZN(ind,:);
                Z  = Z(ind,:);
                H  = H(ind,:);
                
                
                %error
                ylab = sprintf('\\DeltaF/F + %s',errorTrans);
                switch lower(errorTrans)
                    case ''
                        ylab= '\DeltaF/F'; %correct (without '+')
                        err = [];
                        ye  = [];
                    case 'std'
                        err = res.trans.std(:,tra);
                        ye  = [y-err ; y(indR)+err(indR)];
                    case 'sem'
                        err = res.trans.std(:,tra)/sqrt(N);
                        ye   = [y-err ; y(indR)+err(indR)];
                    otherwise
                        error('errorTrans = ''%s'' is NOT implemented',errorTrans)
                end
                %min/max
                mima1 = [min([mima1(1);y(:)]),max([mima1(2);y(:)])];
                mima2 = [min([mima2(1);ye(:)]),max([mima2(2);ye(:)])];
                mima3 = [min([mima3(1);Z(:)]),max([mima3(2);Z(:)])];
                
                %HEAT MAP
                set(hf(3),'CurrentAxes',ha(1,tra))
                %sort & normalize
                ZN = bsxfun(@times,Z,1./max(Z,[],2));
                [~,ind] = sort(mean(ZN,2));
                ZN = ZN(ind,:);
                %plot
                imagesc(t,1:N,ZN);
                axis xy
                %text
                title(res.trans.label{tra})
                ylabel({sprintf('N = %i',N),'Normalized/Sorted'})
                
                %MEAN DATA
                set(hf(3),'CurrentAxes',ha(2,tra))
                xline(0); hold on;
                hp = plot(t,y,props.plo.tra{:});
                %plot error
                if ~isempty(err)
                    h  = fill(te,ye,get(hp,'color'),'facealpha',0.5,...
                        'EdgeColor','none');
                    uistack(h,'bottom')
                end
                %text
                xlabel('Time [s]')
                ylabel(ylab)
                %settings
                set(gca,'xlim',[-margin,margin],props.axi.tra{:})
                
                %HYPNOGRAM
                set(hf(3),'CurrentAxes',ha(3,tra))
                imagesc(t,1:N,H,stage.nums([1,end]));
                axis xy
                title('Stage Map')
                xlabel('Time [s]')
                ylabel('N')
            end
            %colorbar I
            set(hf(3),'CurrentAxes',ha(1,end))
            pos = get(gca,'position');
            h = colorbar;
            set(gca,'position',pos)
            set(ha(1,:),'xticklabel',[],'clim',mima2)
            
            %colorbar II
            set(hf(3),'CurrentAxes',ha(3,end))
            pos = get(gca,'position');
            h = colorbar;
            h.Ticks = stage.nums;
            h.TickLabels = stage.labs;
            set(gca,'position',pos)
            
            %text/settings
            sgtitle({id,'Transitions'})
            linkaxes(ha,'x')
            %set(ha(1,:),'xticklabel',[],'clim',mima2)
            set(ha(1,:),'xticklabel',[],'clim',[0,1])            
            if props.equal.ylimTRA
                set(ha(2,:),'ylim',mima2)
            end
    end
    
    %% SAVE
    drawnow
    if opt.saveDo
        %save path
        if ischar(opt.savePath) && exist(opt.savePath,'dir')
            sPath = opt.savePath;
            uniqueFiles = true;
        else
            sPath = rPath;
            uniqueFiles = false;
        end
   
        %savenames
        snames = cell(noFIG+1,1);
        if ~uniqueFiles
            snames{1} = fullfile(sPath,sprintf('%s.mat',sFile));
            for k = 1:noFIG %no extension!
                snames{k+1} = fullfile(sPath,...
                    sprintf('%s_%s',sFile,figLab{k}));
            end
        else
            cnt = 0; tmp = true; %init
            while tmp
                cnt = cnt+1;
                snames{1} = fullfile(sPath,sprintf('%s_%i.mat',...
                    sFile,cnt));
                for k = 1:noFIG %no extension!
                    snames{k+1} = fullfile(sPath,sprintf('%s_%s_%i',...
                        sFile,figLab{k},cnt));
                end
                tmp = exist(snames{1},'file')==2;
            end
        end
        %save data
        sname = snames{1};
        save(sname,'-struct','res')
        [~,sFile,sExt] = fileparts(sname);
        fprintf('%s Saved Data  : %s\n',indent,[sFile,sExt])
        %save figures
        for fig = 1:noFIG
            sname = snames{fig+1};
            for k = 1:numel(opt.saveFuns)
                fun = opt.saveFuns{k};
                fun(hf(fig),sname)
            end
            [~,sFile] = fileparts(sname);
            fprintf('%s Saved Figure: %s\n',indent,sFile)
        end
        if fil<noFIL %keep only last ;-)
            close(hf)
        end
    end
end