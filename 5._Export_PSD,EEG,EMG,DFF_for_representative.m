clc; clear; close all;
addpath(genpath('Z:\OptoLab_v4.1\function'))
scriptName = mfilename;
scriptPath = fileparts(which(scriptName));
addpath(fullfile(scriptPath,'functions'));

%PARAMETERS
%----------
%READ PATH (select one)
tmp = dir(['']);
rPaths = {tmp.folder};
rPaths = selectionList(rPaths);

%OPTIONS
%psd windows (and sub-windows for pwelch)
spec.win = 30; %[s],
spec.subWin = 5; %PS: spec.win must be dividable by spec.subWin !!!
%limits time & frequency to plot
unit = 's'; %time unit (based on factor below)
fac = 1; %factor to seconds to get unit
lim.f = [0,20]; %freq limit for PSD (f_PSD>=lim.f(1) &  f_PSD<=lim.f(2))
lim.t = [0,28800]; %time in given units (t>lim.t(1) &  t<=lim.t(2))
%NOTE:
% if lim.f and/or lim.t are too large, table cannot be saved to xls-file


%DEFAULT SAVE PATH (askes again)
savePath = ''; %current path if empty
%savePath = 'C:\Users\i0325207\Desktop\test'; %current path if empty

%PLOT PROPERTIES
props.figure = {'position',[564 408 718 530]};
dXtick = 15; %auto if empty

%MAIN SCRIPT
%-----------
%check
tmp = spec.win/spec.subWin;
if tmp~=round(tmp)
    error('spec.win must be dividable by spec.subWin')
end


%path loop
noPAT = numel(rPaths);
nnPAT = numel(num2str(noPAT));
indent = blanks(2*nnPAT+2);
for pat = 1:noPAT
    rPath = rPaths{pat};
    [~,label] = fileparts(rPath);
    fprintf('%*i/%i: %s\n',nnPAT,pat,noPAT,rPath)
    
    %Channels of PSD
    fileEEG = fullfile(rPath,'EEG2.mat');
    fileEMG = fullfile(rPath,'EMG.mat');
    fileDFF  = fullfile(rPath,'Ca_DFF.mat');
    fileHYP  = fullfile(rPath,'Hypnogram.mat');
    
    %LOAD
    %EEG / EMG
    tmp = load(fileEEG);
    fsEEG = tmp.SampRate;
    EEG   = tmp.resampled_data_mV;
    tmp = load(fileEMG);
    EMG   = tmp.resampled_data_mV;
    if fsEEG~=tmp.SampRate
        error('Sample Rates EEG and EMG must be the same')
    end
    if numel(EMG)~=numel(EEG)
        error('something is wrong')
    end
    %DFF
    tmp = load(fileDFF);
    if isfield(tmp,'fs')
        fsDFF = tmp.fs;
    else
        fsDFF = tmp.stages.SampRate;
    end
    if isfield(tmp,'offset')
        off = tmp.offset;
    elseif isfield(tmp.info,'offset')        
        off = tmp.info.offset;
    else
        off = 0; 
    end
    if off~=0 
        warning('Re-check code for offset!')
    end
    DFF = [NaN(off*fsDFF,1);tmp.dff];
    %hypnogram
    tmp = load(fileHYP);
    hypnogram = tmp.Hypnogram; %same fs as EEG
    if numel(EEG)~=numel(hypnogram)
        error('something is wrong')
    end
    %time vectors
    t_EEG = (1:numel(EEG))/fsEEG*fac;
    t_DFF = (1:numel(DFF))/fsDFF*fac;
    
    %CUT DATA TO SELECTED TIME FRAME
    ind = t_EEG>lim.t(1) & t_EEG<=lim.t(2);
    t_EEG(~ind) = [];
    EEG(~ind) = [];
    EMG(~ind) = [];
    hypnogram(~ind) = [];
    
    
    ind = t_DFF>lim.t(1) & t_DFF<=lim.t(2);
    t_DFF(~ind) = [];
    DFF(~ind) = [];
    
    %PSD
    %reshape EEG
    rows = spec.win*fsEEG;
    cols = floor(numel(EEG)/rows);
    data = reshape(EEG(1:rows*cols),[rows,cols]);
    [PXX,f_PSD] = pwelch(data,rows/spec.subWin,0,[],fsEEG);
    %time vector
    dt = spec.win*fac;
    t_PSD = (1:size(PXX,2))*dt-dt/2 + lim.t(1);
    %limit frequences
    ind = f_PSD>=lim.f(1) &  f_PSD<=lim.f(2);
    PXX(~ind,:) = [];
    f_PSD(~ind) = [];
   
 
    %% PLOT
    %----------------------------------------------------------
    %FIGURE / AXES
    hf = figure(props.figure{:});
    ha = fig_createAxes(hf,[4,1],[0.1,0,0.1],[0.1,0.01,0.1],'normalized');
    %hypnogram axis
    ha(5) = axes('unit','normalized','position',get(ha(3),'position'));
    
    %PSD
    set(hf,'currentaxes',ha(1));
    imagesc(t_PSD,f_PSD,10*log10(PXX),[-50,-20]); 
    axis xy; colormap(jet)
    title({strrep(label,'_','\_'),'Power Spectra EEG'})
    ylabel('Frequency [Hz]');
  
    
    %EEG
    set(hf,'currentaxes',ha(2));
    plot(t_EEG,EEG)
    ylabel('EEG');
    
    %EMG
    set(hf,'currentaxes',ha(3));
    plot(t_EEG,EMG)
    ylabel('EMG');
    
    %HYPNOGRAM
    set(hf,'currentaxes',ha(5));
    plot(t_EEG,hypnogram,'k'); %same size as EEG for shure
    set(gca,'YAxisLocation','right','ylim',[0.5,8.5],'ytick',1:3,...
        'yticklabel',{'Wake','NREM','REM'},'color','none')
    
    %DFF
    set(hf,'currentaxes',ha(4));
    plot(t_DFF,DFF)
    ylabel('DFF');
    xlabel(sprintf('Time [%s]',unit))
    set(gca,'ylim',[0,1])
    
    %SETTINGS
    linkaxes(ha,'x')
    set(ha,'xlim',lim.t)
    set(ha(1),'ytick',0:5:lim.f(2))
    ind = [1:3,5];
    set(ha(ind),'xticklabel',[])
    set(ha(2:end),'box','off')
    if ~isempty(dXtick)
        set(ha,'xtick',0:dXtick:t_EEG(end))
    end
    
    
    %% SAVE
    %save figures
    sname = fullfile(sPath,sFile);
    print(hf,sname,'-dpng','-r300')
    %print out
    fprintf('%s Table & Figure Saved\n',indent)
    fprintf('%s   %s*\n',indent,sname)
    
    %maybe close figures if ...
    if noPAT>20 && pat<noPAT
        close(hf);
    end
end %path loop