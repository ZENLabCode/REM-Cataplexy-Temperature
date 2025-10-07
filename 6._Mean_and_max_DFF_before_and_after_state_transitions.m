clc; clear; close all;
addpath(genpath('Z:\OptoLab_v4.1\function'))
%PARAMETERS
%----------

searchPath = [''];

filename.read = 'Ca_DFF.mat'; 
filename.save = 'Ca_DFF'; 

%OPTIONS
opt.err = 'SEM';
opt.eqYlim = false;
%heatmap color limit (1,2 or 3)
%  - 0: all transitions scaled individually, from min to max
%       PS: - other options scale all heatmaps equally
%  - 1: min/max of all transtions mean data
%  - 2: min/max of all transtions mean data +- error
%  - 3: min/max of all (not of mean data)
opt.clim = 2;

opt.remNaN = true;

%SAVE OPTIONS
%  - savePath: if empty, uses selected paths as save path.
%                 (concatenates data of sub-paths of selected paths)
%              if NOT empty, saves all data from selected paths into this
%              one. Therefore does not overwrite any file but creates
%              unique savenames (numbering)
%  - save: do or not (true of false) 
%  - saveFun: save functions, several possible
opt.savePath = '';
opt.save     = true;
opt.saveFuns = {... input for all: @(figure handle,savename)
    @(h,name)saveas(h,name);... as *.fig
    @(h,name)print(h,name,'-r300','-dpng');...
    ...@(h,name)print(h,name,'-r300','-dpdf');...
    };


%PLOT PROPERTIES
%figure
props.fig = {};
%axis
props.axiMap = {}; %heat maps
props.axiLin =  {'box','off'}; %line plots
%plot
props.ploLin = {'linewidth',2}; %line plot



%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%SELECT DATA
rFile = 'Ca_DFF.mat'; 
if ~exist(searchPath,'dir')
    error('searchPath does NOT exist:\n%s',searchPath)
end
tmp = dir(fullfile(searchPath,'**',rFile));
tmp = {tmp.folder};
tmp = unique(cellfun(@fileparts,tmp,'uniformoutput',false));
bPaths = selectionList(tmp);
if isempty(bPaths)
    fprintf('2,No Path Selected\n')
    return
end

%INIT
tmp = fieldnames(props);
for k = 1:numel(tmp)
    if isempty(props.(tmp{k}))
        props.(tmp{k}) = {'visible','on'};
    end
end

%% PATH LOOP
noPAT = numel(bPaths);
nnPAT = numel(num2str(noPAT));
indent = blanks(2*nnPAT+2);
for pat = 1:noPAT
    bPath = bPaths{pat};
    fprintf('%*i/%i: %s\n',nnPAT,pat,noPAT,bPath)
    %id for title
    [~,id] = fileparts(bPath);
    id = strrep(id,'_','\_');
    
    %find all files in bPath
    tmp = dir(fullfile(bPath,'**',rFile));
    rnames = fullfile({tmp.folder},{tmp.name});
    noHRS  = numel(rnames);
    hours  = cell(size(rnames));
    for  hrs = 1:noHRS
        [~,hours{hrs}] = fileparts(fileparts(rnames{hrs}));
    end
    fprintf('%s Concat (N = %i): %s\n',indent,noHRS,strjoin(hours,', '))
    
    %HOURS LOOP
    for hrs = 1:numel(hours)
        hour  = hours{hrs};
        rname = rnames{hrs};
        data = load(rname);
        if ~isfield(data.trans,'Data')
            error('Data is not cutted properly')
        end
        
        %INIT
        if hrs==1
            t      = data.trans.t;
            labels = data.trans.label;
            noTRA  = numel(labels);
            noTIM  = numel(t);
            Data   = cell(1,noTRA);
            Data(:) = {NaN(0,noTIM)};
        elseif ~isequal(t,data.trans.t) || ...
                ~isequal(labels,data.trans.label)
            error('uups')
        end
        
        %APPEND DATA
        for tra = 1:noTRA
            tmp = data.trans.Data{tra};
            Data{tra} = [Data{tra};tmp'];
        end
    end
    
    %FIGURE/AXES
    dx  = [50,50,50];
    dy  = [50,5,70];
    pos = [10,10,noTRA*120+[1,noTRA-1,1]*dx(:),2*140+sum(dy)];
    hf = figure('position',pos);
    movegui(gcf,'center'); drawnow
    set(hf,props.fig{:});  drawnow
    ha = fig_createAxes(gcf,[2,noTRA],dx,dy,'pixel');
    set(ha,'unit','normalized')
    if ~exist('HA','var')
        HA = NaN(2,noTRA,noPAT);
    end
    HA(:,:,pat) = ha;
    
    %PLOT
    indR = noTIM:-1:1; %reverse index
    te = [t;t(indR)]; %for error plots
    mima1 = [inf,-inf]; %mean
    mima2 = [inf,-inf]; %mean +- error
    mima3 = [inf,-inf]; %all transitions
    for tra = 1:noTRA
        %GET DATA
        Z = Data{tra};
        if opt.remNaN
            ind = isnan(sum(Z,2));
            Z(ind,:) = [];
        end
        N = size(Z,1);
        y = nanmean(Z,1);
        label = labels{tra};
        %error & label
        ylab = sprintf('\\DeltaF/F + %s',opt.err);
        switch lower(opt.err)
            case ''
                ylab= '\DeltaF/F';
                err = [];
                ye  = [];
            case 'std'
                err = nanstd(Z,[],1);
                ye  = [y-err,y(indR)+err(indR)];
            case 'sem'
                err = nanstd(Z,[],1)/sqrt(N);
                ye  = [y-err,y(indR)+err(indR)];
            otherwise
                error(['opt.err = ''%s'' is NOT ',...
                    'implemented'],opt.err)
        end
        %sort Z
        % [~,ind] = sort(mean(Z,2));
        [~,ind] = sort(max(Z,[],2));
        Z = Z(ind,:);
        %min/max
        mima1 = [min([mima1(1);y(:)]),max([mima1(2);y(:)])];
        mima2 = [min([mima2(1);ye(:)]),max([mima2(2);ye(:)])];
        mima3 = [min([mima3(1);Z(:)]),max([mima3(2);Z(:)])];
        
        %HEAT MAP
        set(hf,'CurrentAxes',ha(1,tra))
        imagesc(t,1:N,Z);
        axis xy
        %text
        title(label)
        ylabel(sprintf('N = %i',N))
        
        %MEAN DATA
        set(hf,'CurrentAxes',ha(2,tra))
        xline(0); hold on;
        hp = plot(t,y,props.ploLin{:});
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
        set(gca,'xlim',t([1,end]))
        
        
        %TABLE
        %init
        if ~exist('labelsTab','var')
            durB = abs(min(t)); %before
            durA = max(t); %after
            labelsTab = {'Transition',...
                sprintf('mean %g-s before',durB),...
                sprintf('max %g-s before',durB),...
                sprintf('mean %g-s after',durA),...
                sprintf('max %g-s after',durA),...
                'Data'};
        end
        if ~exist('Tab','var')
            Tab = cell(0,numel(labelsTab));
        end
        %fill
        ind1 = t>-durB  & t<=0;    %PS: zero belongs to before
        ind2 = t> 0     & t<=durA; 
        Tab(end+1,:) = {label ,...
            mean(y(ind1)), max(y(ind1)),...
            mean(y(ind2)), max(y(ind2)),strjoin(hours,', ')};
    end
    
    %colorbar
    set(hf,'CurrentAxes',ha(1,end))
    pos = get(gca,'position');
    hc = colorbar;
    set(gca,'position',pos)
    %text
    sgtitle({id,strrep(strjoin(hours,', '),'_','\_')})
    
    %SETTINGS
    %heatmaps
    switch opt.clim
        case 0
            hc.Ticks = hc.Limits;
            hc.TickLabels = {'min','max'};
        case 1
            cl = mima1;
        case 2
            cl = mima2;
        case 3
            cl = mima3;
        otherwise
            error('opt.clim = %i is NOT omplemented',opt.clim)
    end
    if opt.clim>0
        set(ha(1,:),'clim',cl);
        
    end
    set(ha(1,:),'xticklabel',[],props.axiMap{:})
    %line plots
    set(ha(2,:),props.axiLin{:})
    if opt.eqYlim
        if isempty(ye)
            yl = mima1 + [-1,1]*0.1*diff(mima1);
        else
            yl = mima2;
        end
        set(ha(2,:),'ylim',yl)
        
    end
    %all
    linkaxes(ha,'x')
    HA(:,:,pat) = ha; 
    
    %APPEND DATA FOR SAVING
    if ~exist('DD','var')
        DD = cell(noPAT,3);
    end
    [~,sFile] = fileparts(rFile);
    cnt   = 0;
    if ischar(opt.savePath) && exist(opt.savePath,'dir')
        sname = fullfile(opt.savePath,sFile);
        cnt = 0;
        while exist([sname,'.xslx'],'file')==2
            cnt = cnt+1; 
            sname = fullfile(opt.savePath,sprintf('%s_%i',sFile,cnt));
        end
    else
        sname = fullfile(bPath,sFile); 
    end
    DD(pat,:) = {hf,Tab,sname};
    clear Tab
end %path loop

%SAVE DATA
%% SAVE
drawnow
if opt.save
    %save path
    if ischar(opt.savePath) && exist(opt.savePath,'dir')
        sPath = opt.savePath;
        uniqueFiles = true;
        
    else
        sPath = searchPath;
        uniqueFiles = false;
    end
    %check
    noFIG = numel(hf);
    [~,sFile] = fileparts(rFile);
    sFile = sprintf('%s_DFF',sFile);
    
    %savenames
    snames = cell(noFIG+1,1);
    if ~uniqueFiles
        snames{1} = fullfile(sPath,sprintf('%s.mat',sFile));
    end
    for pat = 1:noPAT
        [hf,Tab,sname] = DD{pat,:};
        for k = 1:numel(opt.saveFuns)
            fun = opt.saveFuns{k};
            fun(hf,sname)
        end
        %table
        sname = [sname,'.xlsx'];
        if exist(sname,'file')==2
            delete(sname)
        end
        xlswrite(sname,[labelsTab;Tab])
        fprintf('Saved: %s\n',sname)
        try
            xls_cellFit(sname)
        catch
        end
    end
else
    fprintf(2,'Nothing Saved!\n')
end
