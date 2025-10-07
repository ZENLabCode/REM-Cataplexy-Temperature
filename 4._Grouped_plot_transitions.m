%     --> use first cut_DatabyHours.m
clc; clear; close all;
addpath(genpath('Z:\OptoLab_v4.1\function'));

%PARAMETERS
%----------
%SELECT DATA IDs (name of last sub-folder)
% - Opens selection list for to paths having the file 'Ca_DFF.mat'
%     PS: Selection list only shows the paths of first hour-folder in
%         Hours (e.g sub-folder 'hour1')
%         Script then concats data of corresponding hour folders
%searchPath = [... searches files here and in all it's sub-paths, animal folders in it
 searchPath = [''];


%HOUR SUB-FOLDERS TO CONCAT DATA
Hours = {... concats Hours for each cell in cell
    %{'hour5','hour6','hour7','hour8'};...
    %{'hour1','hour2','hour3','hour4'};...
  %{'hour1','hour2','hour3','hour4','hour5','hour6','hour7','hour8'};
    {'hour1.30','hour1.60','hour2.30','hour2.60','hour3.30','hour3.60','hour4.30','hour4.60','hour5.30','hour5.60','hour6.30','hour6.60','hour7.30','hour7.60','hour8.30','hour8.60'};...
    %{'hour1.30','hour1.60','hour3.30','hour3.60','hour5.30','hour5.60','hour7.30','hour7.60'};...
    %{'hour2.30','hour2.60','hour4.30','hour4.60','hour6.30','hour6.60','hour8.30','hour8.60'};...
    };

%animals = {'106','119','120','N11','L11','X11','Y11','Z11'};
% animals    = {'S11'}
recordings = cellfun(@(x)sprintf('%s_1',x),animals,'uniformoutput',false);

recordings  = animals; 
%hours   = {'hour5','hour6','hour7','hour8'};
%hours   = {'hour1','hour2','hour3','hour4'};
%hours   = {'hour1','hour2','hour3','hour4','hour5','hour6','hour7','hour8'};
hours   = {'hour1.30','hour1.60','hour2.30','hour2.60','hour3.30','hour3.60','hour4.30','hour4.60','hour5.30','hour5.60','hour6.30','hour6.60','hour7.30','hour7.60','hour8.30','hour8.60'};
%hours   = {'hour1.30','hour1.60','hour3.30','hour3.60','hour5.30','hour5.60','hour7.30','hour7.60'};...
%hours   = {'hour2.30','hour2.60','hour4.30','hour4.60','hour6.30','hour6.60','hour8.30','hour8.60'};...

%OPTIONS
%line plot errors
%  - 'SEM', 'STD' or '' (no error plot)
%  - PS: STD was huge, so SEM is probybly better
opt.err = 'SEM';
%equal ylim of all line plots, true or false
%  - PS: some trnasitions had low data compared to others. With equal ylim,
%        they looked close to horizontal lines. So false might be better
opt.eqYlim = true;
%heatmap color limit (1,2 or 3)
%  - 0: all transitions scaled individually, from min to max
%       PS: - other options scale all heatmaps equally
%  - 1: min/max of all transtions mean data
%  - 2: min/max of all transtions mean data +- error
%  - 3: min/max of all (not of mean data)
opt.clim = 2;
%remove transitions that are cropped because of the hour cut
opt.remNaN = true;


%PLOT PROPERTIES
%figure
props.fig = {};
%axis
props.axiMap = {}; %heat maps
props.axiLin =  {'box','off'}; %line plots
%plot
props.ploLin = {'linewidth',2}; %line plot

stageNum = 1:4;
stageLab = {'Wake','NREM','REM','Cataplexy'};

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))
%non-empty properties
tmp = fieldnames(props);
for k = 1:numel(tmp)
    if isempty(props.(tmp{k}))
        props.(tmp{k}) = {'visible','on'};
    end
end

%GET FILES
if ~exist(searchPath,'dir')
    error('searchPath does NOT exist:\n%s',searchPath)
end
tmp   = dir(fullfile(searchPath,'**','Ca_DFF.mat'));
tmp([tmp.isdir])= [];
Fnames = fullfile({tmp.folder},{tmp.name})';
ind = false(size(Fnames));
for k = 1:numel(ind)
    tmp = lower(strsplit(Fnames{k},'\'));
    ind(k) = ...
        any(ismember(tmp,lower(recordings))) && ...
        any(ismember(tmp,lower(hours)));
end
fnames = Fnames(ind);
noFIL = numel(fnames);

%FILE LOOP
clear res
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    fname = fnames{fil};
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,fname)
    data = load(fname);
    if ~isfield(data.trans,'Data')
        error('Data is not cutted properly')
    end
    %animal id
    tmp    = lower(strsplit(fname,'\'));
    ind    = ismember(lower(animals),tmp);
    animal = animals{ind};
    ind    = ismember(lower(hours),tmp);
    hour   = hours{ind};
    fprintf('%s %s, %s\n',indent,animal,hour)

    %INIT
    if ~exist('res','var')
        res.t        = data.trans.t;
        res.labelTRA = data.trans.label;
        res.noTRA    = numel(res.labelTRA);
        res.noTIM    = numel(res.t);
        res.Data        = cell(1,res.noTRA);
        res.Data(:)     = {NaN(0,res.noTIM)};
        res.Hypno(:)   = cell(1,res.noTRA);
        res.Hypno(:)   = {NaN(0,res.noTIM)};
        res.labelREC    = cell(1,res.noTRA);
        res.labelREC(:) = {cell(0,1)};
    elseif ~isequal(res.t, data.trans.t) || ...
            ~isequal(res.labelTRA,data.trans.label)
        error('uups')
    end

    %APPEND DATA
    for tra = 1:res.noTRA
        dat = (data.trans.Data{tra})';
        if isfield(data.trans,'Hypno')
            hyp = (data.trans.Hypno{tra})';
        else
            hyp = NaN(size(dat));
        end
        ind = (1:size(dat,1)) + numel(res.labelREC{tra});
        %append data
        tmp = res.Data{tra};
        tmp(ind,:) = dat;
        res.Data{tra} = tmp;
        %append label
        tmp = res.labelREC{tra};
        tmp(ind,1) = {animal};
        res.labelREC{tra} = tmp;
        %append stages
        tmp = res.Hypno{tra};
        tmp(ind,:) = hyp;
        res.Hypno{tra} = tmp;
    end
end

%SORT DATA
for tra = 1:res.noTRA
    lab = res.labelREC{tra};
    dat = res.Data{tra};
    
    %remove NaNs
    if opt.remNaN
        ind = isnan(sum(dat,2));
        dat(ind,:) = [];
    end
    
    %sort by data
    %[~,ind] = sort(mean(dat,2));
    [~,ind] = sort(max(dat,[],2));
    lab = lab(ind);
    dat = dat(ind,:);

    %sort by animal
    [~,ind] = sort(lab);
    lab = lab(ind);
    dat = dat(ind,:);

    %append
    res.labelREC{tra} = lab;
    res.Data{tra}    = dat;
end


%FIGURE/AXES
dx  = [50,50,70];
dy  = [50,5,70];
pos = [10,10,res.noTRA*120+[1,res.noTRA-1,1]*dx(:),3*140+sum(dy)];
hf = figure('position',pos);
movegui(gcf,'center'); drawnow
set(hf,props.fig{:});  drawnow
ha = fig_createAxes(gcf,[2,res.noTRA],dx,dy,'pixel');
set(ha,'unit','normalized')

%PLOT
indR = res.noTIM:-1:1; %reverse index
te = [res.t;res.t(indR)]; %for error plots
mima1 = [inf,-inf]; %mean
mima2 = [inf,-inf]; %mean +- error
mima3 = [inf,-inf]; %all transitions
for tra = 1:res.noTRA
    %GET DAT
    dat = res.Data{tra};
    hyp = res.Hypno{tra};
    label = res.labelTRA{tra};
    labs  = res.labelREC{tra};
    %by animal
    ani.animals = unique(labs);
    noANI    = numel(ani.animals);
    ani.ind = cell(noANI,1);
    ani.dat = NaN(noANI,size(dat,2));
    for k = 1:numel(ani.animals)
        animal = ani.animals{k};
        ind    = strcmpi(labs,animal);
        ani.ind{k} = ind;
        ani.dat(k,:) = nanmean(dat(ind,:),1);
    end

    %MEAN ACRISS ANIMALS
    N    = numel(ani.animals);
    mDat = nanmean(ani.dat,1);
    %error & label
    ylab = sprintf('\\DeltaF/F + %s',opt.err);
    switch lower(opt.err)
        case ''
            ylab= '\DeltaF/F';
            err = [];
            ye  = [];
        case 'std'
            err = nanstd(ani.dat,[],1);
            ye  = [mDat-err,mDat(indR)+err(indR)];
        case 'sem'
            err = nanstd(ani.dat,[],1)/sqrt(N);
            ye  = [mDat-err,mDat(indR)+err(indR)];
        otherwise
            error(['opt.err = ''%s'' is NOT ',...
                'implemented'],opt.err)
    end
    %min/max
    mima1 = [min([mima1(1);mDat(:)]),max([mima1(2);mDat(:)])];
    mima2 = [min([mima2(1);ye(:)]),max([mima2(2);ye(:)])];
    mima3 = [min([mima3(1);dat(:)]),max([mima3(2);dat(:)])];

    %HEAT MAP
    set(hf,'CurrentAxes',ha(1,tra))
    imagesc(res.t,1:size(dat,1),dat); hold on
    axis xy
    %text
    title(label)
    ylabel(sprintf('N = %i',size(dat,1)))
    tmp = unique(labs);
    for k = 1:noANI
        ind = find(ani.ind{k},1,'first');

        xx = max(res.t);
        yy = ind-0.5;
        text(xx,yy,sprintf(' %s',ani.animals{k}),...
            'horizontalalignment','left','verticalalignment','bottom');
        h = plot([1,1.4]*xx,[yy,yy],'k','linewidth',1,...
            'clipping','off');
    end

    %MEAN DATA
    set(hf,'CurrentAxes',ha(2,tra))
    xline(0); hold on;
    hp = plot(res.t,mDat,props.ploLin{:});
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
    set(gca,'xlim',res.t([1,end]))
    
end
%colorbar I 
set(hf,'CurrentAxes',ha(1,end))
pos = get(gca,'position');
hc = colorbar;
set(gca,'position',pos)
pos = get(hc,'Position');
pos(1) = 1-3*pos(3);
set(hc,'position',pos)

sgtitle({strjoin(animals,', '),strjoin(hours,', ')})

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
    set(ha(1,:),'clim',mima2)
elseif opt.clim==4
    set(ha(1,:),'clim',[0.2,8])
end
set(ha(1,:),'xticklabel',[],props.axiMap{:})
%line ploits
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
