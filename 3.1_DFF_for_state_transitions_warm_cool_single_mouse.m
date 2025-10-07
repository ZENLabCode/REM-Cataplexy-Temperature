clc; clear; close all;
addpath(genpath('Z:\OptoLab_v4.1\function'));

%PARAMETERS
%----------
%MAIN PATHS

mPath = [];
%CUT TIME & FOLDER NAME
CUT = {... {t>this, t<=this, folder name}, time in seconds! %     
%     1*1800 , 2*1800, 'hour1.60';...
%     2*1800 , 3*1800, 'hour2.30';...
%     3*1800 , 4*1800, 'hour2.60';...
%     4*1800 , 5*1800, 'hour3.30';...
%     5*1800 , 6*1800, 'hour3.60';...
%     6*1800 , 7*1800, 'hour4.30';...
%     7*1800 , 8*1800, 'hour4.60';...
%     8*1800 , 9*1800, 'hour5.30';...
%     9*1800 , 10*1800, 'hour5.60';...
%     10*1800 , 11*1800, 'hour6.30';...
%     11*1800 , 12*1800, 'hour6.60';...
%     12*1800 , 13*1800, 'hour7.30';...
%     13*1800 , 14*1800, 'hour7.60';...
%     14*1800 , 15*1800, 'hour8.30';...
%     15*1800 , 16*1800, 'hour8.60';...


    % 0*3600 , 1*3600, 'hour1';...
    % 1*3600 , 2*3600, 'hour2';...
    % 2*3600 , 3*3600, 'hour3';...
    % 3*3600 , 4*3600, 'hour4';...
    % 4*3600 , 5*3600, 'hour5';...
    % 5*3600 , 6*3600, 'hour6';...
    % 6*3600 , 7*3600, 'hour7';...
    % 7*3600 , 8*3600, 'hour8';...
    };

%STAGES
Stages = {... {number,label}
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    4,'Cataplexy';...
    };
%TRANSITIONS
Transitions = {... {stage from, stage to}
    'Wake','NREM';...
    'NREM','Wake';...
    'REM','Wake';...
    'NREM','REM';...
    'Wake','Cataplexy';...
    'Cataplexy', 'Wake';...
    };
margin  = 90; %[s], plots transtions +- margin 
exclNaN = false; %exclude transitions with margin beyond data range

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
fprintf('Main Path: %s\n',mPath)
%just in case
if strcmp(mPath(end),filesep)
    mPath(end) = [];
end
%check
if numel(unique(lower(CUT(:,3))))~=size(CUT,1)
    error('Sub-folder names in ''CUT'' must be unique')
end

%DATA PATHS
tmp = dir(fullfile(mPath,'**\Ca_DFF.mat'));
rPaths = {tmp.folder};

ind = false(size(rPaths)); 
for k = 1:numel(rPaths)
    ind(k) = any(ismember(lower(strsplit(rPaths{k},filesep)),...
        lower(CUT(:,3))));
end
rPaths(ind) = [];
if numel(rPaths)==0
    fprintf(2,'No Data Path Found!\n')
    return
end
%select paths
rPaths = selectionList(rPaths);
if numel(rPaths)==0
    fprintf(2,'No Path Selected!\n')
    return
end

%% PLOT SETTINGS
num = cell2mat(Stages(:,1));
num = min(num):max(num);
N = numel(num);
tmp = cell(1,N);
for k = 1:N
    ind = cellfun(@(x)x==num(k),Stages(:,1));
    if isempty(ind)
        tmp{k} = '';
    else
        tmp(k) = Stages(ind,2);
    end
end
mima = num([1,end]);
props.axi.hyp = ['ytick',num,'yticklabel',{tmp},...
    'ylim',mima+[-.2,.2]*diff(mima),props.axi.hyp];
clear num tmp N mima

%number of ...
noPAT = numel(rPaths);
noSTA = size(Stages,1);
noTRA = size(Transitions,1);
noCUT = size(CUT,1);

%PATH LOOP
nnPAT = numel(num2str(noPAT));
indent = blanks(2*nnPAT+2);
for pat = 1:noPAT
    rPath = rPaths{pat};
    [s0,subPath] = fileparts(rPath);
    fprintf('%*i/%i: ...%s%s\n',nnPAT,pat,noPAT,filesep,subPath)
    [~,subPath] = fileparts(rPath);
    %info string
    info = sprintf('Exported by %s.m, %s',scriptName,date);

    %SAVE PATHS
    sPaths = cell(noCUT,1);
    for k = 1:noCUT
        sPath = fullfile(s0,CUT{k,3},subPath);
        if ~exist(sPath,'dir')
            mkdir(sPath)
        end
        sPaths{k} = sPath;
    end

    
    %% Hypnogram
    rFile = 'Hypnogram.mat';
    fprintf('%s %s\n',indent,rFile)
    data  = load(fullfile(rPath,rFile));
    fsHYP = data.fs;
    t = (1:numel(data.Hypnogram))/fsHYP;
    Hypnogram = data.Hypnogram(:); %for Ca_DFF
    for k = 1:noCUT
        [t1,t2,subPI] = CUT{k,:};
        ind = t>t1 & t<=t2;
        if all(~ind)
            break
        end
        %cut
        res = data;
        res.Hypnogram(~ind) = [];
        res.tStart = t1;
        res.info   = info;
        %save
        save(fullfile(sPaths{k},rFile),'-struct','res')
    end %cut loop
    
    %% EEGs & EMGs
    files = {'EEG1','EEG2','EEG3','EEG4','EMG'};
    clear fsEEG;
    for fil = 1:numel(files)
        rFile = [files{fil},'.mat'];
        fprintf('%s %s\n',indent,rFile)
        if exist(fullfile(rPath,rFile),'file')~=2
            fprintf(2,'\b - not found\n')
            continue
        end
        data = load(fullfile(rPath,rFile));
        if ~exist('fsEEG','var')
            fsEEG = data.SampRate;
        elseif fsEEG~=data.SampRate
            error('grrrhhh')
        end
        t = (1:numel(data.resampled_data_mV))/fsEEG;
        for k = 1:noCUT
            [t1,t2,subPI] = CUT{k,:};
            ind = t>t1 & t<=t2;
            if all(~ind)
                break
            end
            %cut
            res = data;
            res.resampled_data_mV(~ind) = [];
            res.tStart    = t1;
            res.info      = info;
            if ~isempty(res.stimTimes) %care later when found
                fprintf('%s   Cutting stimTimes not yet implemented\n',...
                    indent)
                res.stimTimes = [];
            end
            %save
            save(fullfile(sPaths{k},rFile),'-struct','res')
        end %cut loop
    end %file loop
    
    
    %% CA_DFF files
    rFile = 'Ca_DFF.mat';
    fprintf('%s %s\n',indent,rFile)
    data  = load(fullfile(rPath,rFile));
    try
        fsDFF = data.stages.SampRate;
    catch
        fsDFF = data.fs;
    end
    t = (1:numel(data.dff))/fsDFF;
    %id for title
    id = strrep(rFile,'_','\_'); %init
    if isfield(data,'info') && isfield(data.info,'file')
        [~,id] = fileparts(data.info.file); %from recorded filename
        id = strrep(id,'_','\_');
    end
    %upsample hypnogram (chould basically not happen)
    if fsHYP<fsDFF
        fac = fsDFF/fsHYP;
        Hypnogram = repmat(Hypnogram',fac,1);
        Hypnogram = Hypnogram(:);
        fsHYP = fsDFF;
        fprintf('%s   Hypnogram upsampled by factor %g\n',indent,fac)
    elseif fsHYP>fsDFF
        Hypnogram = interp1((1:numel(Hypnogram))/fsHYP,Hypnogram,...
            t,'nearest');
        fprintf('%s   Hypnogram interpolated %g --> %g Hz\n',indent,...
            fsHYP,fsDFF)
        fsHYP = fsDFF;
    end
    Hypnogram(numel(data.dff)+1:end) = [];
    Hypnogram(numel(Hypnogram)+1:numel(data.dff)) = NaN;
    %cut loop
    for cut = 1:noCUT
        [t1,t2,subPI] = CUT{cut,:};
        indT = t>t1 & t<=t2;
        if all(~indT)
            break
        end
        %cut
        res = data;
        res.dff(~indT) = [];
        res.tStart = t1;
        res.info   = info;
        hypnogram  = Hypnogram(indT);
        res.hypnogram = hypnogram(:); %new, was missing

        %MEAN ACROSS STAGES
        res.stages.label = cell(1,noSTA);
        res.stages.mean  = NaN(1,noSTA);
        res.stages.std   = NaN(1,noSTA);
        res.stages.N     = NaN(1,noSTA);
        for sta = 1:noSTA
            [stage,label] = Stages{sta,:};
            ind = ismember(hypnogram,stage);
            res.stages.label{sta} = label;
            res.stages.mean(sta)  = mean(res.dff(ind));
            res.stages.std(sta)   = std(res.dff(ind));
            res.stages.N(sta)     = sum(ind);
        end
        
        %MEAN ACROSS TRANSITIONS
        ind0 = -round(margin*fsDFF):round(margin*fsDFF);
        tt = ind0'/fsDFF;
        res.trans = [];
        res.trans.t = tt;
        res.trans.label = cell(1,noTRA);
        res.trans.Hypno = cell(1,noTRA);
        res.trans.Data  = cell(1,noTRA);
        res.trans.mean  = NaN(numel(tt),noTRA);
        res.trans.std   = NaN(numel(tt),noTRA);
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
            M = NaN(numel(ind0),N);
            H = NaN(numel(ind0),N);
            for k = 1:N
                ind = ind0+indTRA(k); %%%CHECK AGAIN
                ind1 = find(ind>=1,1,'first');    %if beyond data range
                ind2 = find(ind<=numel(res.dff),1,'last');
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
        
        %% FIGURE transitions
        % Plot vatiations (2 so far)
        %   1: line plot only
        %   2: line and heat map
        % clc; try close(2); catch; end
        switch 2 %plot variation
            case 1 %transition lines plot only
                indPLO = indPLO + 1;
                hf(indPLO) = figure(props.figure{:});
                figLab{indPLO} = sprintf('transitions_%s',subPI);
                n1 = floor(sqrt(noTRA)); n2 = ceil(noTRA/n1);
                ha = NaN(n2,n1);
                for k = 1:noTRA
                    ha(k) = subplot(n1,n2,k);
                end
                ha = ha';
                
                %PLOT
                tt = res.trans.t(:);
                mima = [inf,-inf];
                for tra = 1:noTRA
                    %plot lines
                    set(hf(indPLO),'CurrentAxes',ha(tra))
                    xline(0); hold on;
                    y = res.trans.mean(:,tra);
                    hp = plot(tt,y,props.plo.tra{:});
                    %plot error
                    ylab = sprintf('\\DeltaF/F + %s',errorTrans);
                    switch lower(errorTrans)
                        case ''
                            err = [];
                            ylab= '\DeltaF/F';
                        case 'std'
                            err = res.trans.std(:,tra);
                        case 'sem'
                            err = res.trans.std(:,tra)/...
                                sqrt(res.trans.N(tra));
                        otherwise
                            error(['errorTrans = ''%s'' is NOT ',...
                                'implemented'],errorTrans)
                    end
                    if ~isempty(err)
                        n  = numel(tt);
                        y  = [y-err ; y(n:-1:1)+err(n:-1:1)];
                        h  = fill([tt;tt(n:-1:1)],y,get(hp,'color'),...
                            'facealpha',0.5,'EdgeColor','none');
                        uistack(h,'bottom')
                    end
                    mima = [min([mima(1);y(:)]),max([mima(2);y(:)])];
                    %text
                    if tra==1
                        title({sprintf('%s, %s',...
                            strrep(subPath,'_','\_'),...
                            strrep(subPI,'_','\_')),...
                            sprintf('Transition (N = %i)',...
                            res.trans.N(tra)),...
                            res.trans.label{tra}})
                    else
                        title({sprintf('Transition (N = %i)',...
                            res.trans.N(tra)),...
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
            case 2 %transition lines plot only
                indPLO = indPLO + 1;
                dx = [50,50,50];
                dy = [50,5,70];
                posF = [10,10,noTRA*120+[1,noTRA-1,1]*dx(:),2*140+sum(dy)];
                hf(indPLO) = figure('position',[10,10,...
                    noTRA*120+[1,noTRA-1,1]*dx(:) , 2*140+sum(dy)],...
                    props.figure{:});
                movegui(gcf,'center'); drawnow
                figLab{indPLO} = 'transitions';
                ha = fig_createAxes(gcf,[2,noTRA],dx,dy,'pixel');
                set(ha,'unit','normalized')
                
                %PLOT
                indR = numel(res.trans.t):-1:1; %reverse index
                tt  = res.trans.t(:);
                te = [tt;tt(indR)]; %for error plots
                mima1 = [inf,-inf]; %mean
                mima2 = [inf,-inf]; %mean +- error
                mima3 = [inf,-inf]; %all transitions
                for tra = 1:noTRA
                    %GET DATA
                    Z = res.trans.Data{tra}';
                    N = res.trans.N(tra);
                    y = res.trans.mean(:,tra);
                    if N~=size(Z,1)
                        error('uups')
                    end
                    %sort Z
                    % [~,ind] = sort(mean(Z,2));
                    [~,ind] = sort(max(Z,[],2));
                    Z = Z(ind,:);
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
                    set(hf(indPLO),'CurrentAxes',ha(1,tra))
                    imagesc(tt,1:N,Z);
                    axis xy
                    %text
                    title(res.trans.label{tra})
                    ylabel(sprintf('N = %i',N))
                    
                    %MEAN DATA
                    set(hf(indPLO),'CurrentAxes',ha(2,tra))
                    xline(0); hold on;
                    hp = plot(tt,y,props.plo.tra{:});
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
                end
                %colorbar
                set(hf(indPLO),'CurrentAxes',ha(1,end))
                pos = get(gca,'position');
                h = colorbar;
                set(gca,'position',pos)
                %text/settings
                sgtitle({id,'Transitions'})
                linkaxes(ha,'x')
                if any(~isfinite(mima1))
                    mima1 = [0,1];
                end
                set(ha(1,:),'xticklabel',[],'clim',mima1)
                %set(ha(1,:),'xticklabel',[])
                if props.equal.ylimTRA
                    set(ha(2,:),'ylim',mima2)
                end
        end
       
        %% SAVE
        %save
        sPath = sPaths{cut};
        save(fullfile(sPath,rFile),'-struct','res')
        %save figures
        for fig = 1:numel(hf)
            sname = fullfile(sPath,figLab{fig});
            print(hf(fig),sname,'-r300','-dpng')
        end
    end %cut loop
end %path loop