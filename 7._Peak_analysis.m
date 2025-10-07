% Exports peaks detected above calculated threshold per stage bout.
clc; clear; %close all
addpath(genpath('Z:\OptoLab_v4.1'))

%PARAMETERS
%-----------
%FILES
tmp = dir(fullfile(['Ca_DFF.mat']));
files = selectionList(fullfile({tmp.folder},{tmp.name}));


%SAVENAME (xls-file, uses read path)
saveFiles = {'*.xlsx'};

%STAGES TO ANALYSE
Stages = {... {number(s), label, plot color}
    3,  'REM'        ,'b';...
    4,  'cataplexy'  ,'r';...
    };

%THERESHOLD FUNCTION 
peakMinDur = 2; %[s]
minPeakProm_readFromFile = false; %if exist
minPeakProminence = 0.02; %changed according to recordings and states


%OPTIONS
opt.plot = true; %for testing whether data is ok
opt.fs   = 10; %sampling rate in case it's not saved in file
opt.table = 1; %1 or 2 (else saves nothing)

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%remove in-existing files
files(cellfun(@(x)exist(x,'file')~=2,files)) = [];
%number of ...
noFIL = numel(files);
noSTA = size(Stages,1);
if noFIL==0
    fprintf('No File To Analyse!\n')
    return
end

%FILE LOOP
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    file  = files{fil};
    [rPath,rFile] = fileparts(file);
    id = rPath(find(rPath==filesep,1,'last')+1:end);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath)
    
    %READ DATA
    tmp = load(file);
    DFF       = tmp.dff(:);
    hypnogram = tmp.hypnogram(:);
    noSAM     = numel(DFF);
    if isfield(tmp,'fs')
        fs = tmp.fs;
    else
        try
            fs = round(1/mean(diff(tmp.trans.t)));
            fprintf(['%s [\bfs = %i Hz ',...
                '(calculated from time vector)]\b\n'],indent,fs)
        catch
            fs = opt.fs;
            fprintf('%s [\bfs = %i Hz(set using opt.fs)]\b\n',...
                indent,fs)
        end
    end
    if noSAM~=numel(hypnogram)
        tmp = load(fullfile(rPath,'Hypnogram.mat'));
        t1 = (1:noSAM)/fs;
        t2 = (1:numel(tmp.Hypnogram))/tmp.fs;
        hypnogram = interp1(t2,tmp.Hypnogram,t1,'nearest');
        hypnogram = hypnogram(:);
    end
    mea = nanmean(DFF);
    %threshold
    mad = nanmean(abs(DFF-nanmedian(DFF))); 
    ind = DFF<=median(DFF)+1.2*mad;
    thr = 1.2*nanmedian(DFF(ind));
  
    
    %Min Peak Prominence
    fname = fullfile(rPath,saveFiles{2});
    if exist(fname,'file')==2 && minPeakProm_readFromFile
        [~,~,tmp] = xlsread(fname);
        column = strcmpi(tmp(1,:),'Min Peak Prominence');
        if any(column)
            minPeakProminence = tmp{2,column};
        end
    end
    
    %FIND STAGE BOUTS
    %init
    clear Bouts
    tot.maxPKS = 0; 
    tot.boutsN = 0; 
    tot.peaksN = 0;
    tot.maxPKS2 = 0; 
    tot.boutsN2 = 0; 
    tot.peaksN2 = 0;
    
    %stage loop
    for sta = 1:noSTA
        [stage,label,~] = Stages{sta,:};
        ind = ismember(hypnogram,stage);
        indSTA = find([ind;NaN]==1 &[0;ind]==0);
        indEND = find([NaN;ind]==1 &[ind;0]==0)-1;
        noBOU = numel(indSTA);
        %check
        if noBOU~=numel(indEND)
            error('grrhhh')
        end
        %test plot, for if something seems to be weird
        if false
            plot(hypnogram,'b','marker','.'); hold on
            plot(indSTA,hypnogram(indSTA),'go')
            plot(indEND,hypnogram(indEND),'ro')
            tmp = sprintf('Hypnogram and %s-Bouts',label);
            if fig==1
                title({tmp,'(using only threshold)'})
            else
                title({tmp,'(using threshold and findpeaks)'})
            end
            ylim([0.5,4.5]); zoom xon;
            retur
        end
        
        %APPEND
        clear res
        if noBOU==0
            res.ind  = [];
            res.dur  = [];
            res.mean = [];
            res.peaks.N   = [];
            res.peaks.height  = [];
            res.peaks.area = [];
            res.peaks.sta  = [];
            res.peaks.ind  = [];
            res.peaks.end  = [];
            res.peaks.val  = [];
            res.peaks2.N   = [];
            res.peaks2.height  = [];
            res.peaks2.area = [];
            res.peaks2.sta  = [];
            res.peaks2.ind  = [];
            res.peaks2.end  = [];
            res.peaks2.val  = [];
        else
            res.ind  = [indSTA(:) indEND(:)];
            res.dur  = (diff(res.ind,[],2)+1)/fs;
            res.mean = NaN(noBOU,1);
            res.peaks(noBOU).N    = [];
            res.peaks(noBOU).height  = [];
            res.peaks(noBOU).area = [];
            res.peaks(noBOU).sta  = [];
            res.peaks(noBOU).ind  = [];
            res.peaks(noBOU).end  = [];
            res.peaks(noBOU).val  = [];
            res.peaks2(noBOU).N    = [];
            res.peaks2(noBOU).height  = [];
            res.peaks2(noBOU).area = [];
            res.peaks2(noBOU).sta  = [];
            res.peaks2(noBOU).ind  = [];
            res.peaks2(noBOU).end  = [];
            res.peaks2(noBOU).val  = [];
        end
        %fine peaks
        for bou = 1:noBOU
            ind1 = indSTA(bou);
            ind2 = indEND(bou);
            indB = false(size(DFF));
            indB(ind1:ind2) = true;
            res.mean(bou) = nanmean(DFF(ind1:ind2));
           
            %FIND PEAKS
            tmp = 1:numel(DFF);
            [pks,locs,widths,proms] = findpeaks(DFF(indB),tmp(indB),...
                'MinPeakProminence',minPeakProminence,...
                'MinPeakHeight',thr);
            %find start/end of each peak
            noPKS = numel(pks);
            pksSTA = NaN(noPKS,1);
            pksEND = NaN(noPKS,1);
            pksARE = NaN(noPKS,1);
            tmp = [find(indB,1,'first'), find(indB,1,'last')];
            for k = 1:noPKS
                pksSTA(k) = max([tmp(1),...
                    find(DFF(1:locs(k))<thr,1,'last')+1]);
                pksEND(k) = min([tmp(2),...
                    find(DFF(locs(k):end)<thr,1,'first')-2+locs(k)]);
                pksARE(k) = sum(DFF(pksSTA(k):pksEND(k))-thr)/fs;
            end
            %split overlapping peaks
            for k = 2:noPKS
                if pksSTA(k)<pksEND(k-1)
                    [~,ind] = min(DFF(locs(k-1):locs(k)));
                    ind = ind+locs(k-1)-1;
                    if pks(k-1)>=pks(k)
                        pksEND(k-1) = ind;
                        pksSTA(k)   = ind+1;
                    else
                        pksEND(k-1) = ind-1;
                        pksSTA(k)   = ind;
                    end
                end
            end
            
            %append
            res.peaks2(bou).N       = numel(pks);
            res.peaks2(bou).height  = pks-thr;
            res.peaks2(bou).area    = pksARE;
            res.peaks2(bou).sta     = pksSTA;
            res.peaks2(bou).ind     = locs;
            res.peaks2(bou).end     = pksEND;
            res.peaks2(bou).val     = pks; 
        end %bout loop
        tmp  = [res.peaks(:).N];
        tot.maxPKS = max([tot.maxPKS;tmp(:)]);
        tot.boutsN = tot.boutsN + noBOU;
        tot.peaksN = tot.peaksN + sum(tmp);
        tmp  = [res.peaks2(:).N];
        tot.maxPKS2 = max([tot.maxPKS2;tmp(:)]);
        tot.boutsN2 = tot.boutsN2 + noBOU;
        tot.peaksN2 = tot.peaksN2 + sum(tmp);
        
        %APPEND
        Bouts.(label) = res;
    end
    
    %% FIGURE
    HF = NaN(2,1);
    for fig = 1:2
        if opt.plot
            hf = figure;
            HF(fig) = hf;
            lenY    = 0.25;
            ha = [...
                axes(hf,'position',[0.1  0.12+2*lenY  0.85  1*lenY]),...
                axes(hf,'position',[0.1  0.11+0*lenY  0.85  2*lenY])];
            %time vector
            t = (1:noSAM)/fs;
            t = t/3600; tUNI = 'h'
            xl = [0,t(end)];
            
            %PLOT HYPNOGRAM
            set(hf,'currentaxes',ha(1))
            plot(t,hypnogram,'k');
            
            
            tmp = strrep(id,'_','\_');
            if fig==1
                title({tmp,'(using only threshold)'})
            else
                title({tmp,'(using threshold and findpeaks)'})
            end
            ylabel('Stages')
            tmp = max(hypnogram);
            set(gca,'xticklabel',[],'ylim',[0.5,1.2*tmp],'ytick',1:tmp);
            
            %PLOT DFF
            set(hf,'currentaxes',ha(2))
            hp     = plot(t,DFF,'k'); hold on
            legSTR = {'DFF'};
            xlabel(sprintf('Time [%s]',tUNI)); ylabel('DFF')
            %mean & threshold
            plot(xl,repmat(mea,size(xl)),'k','linewidth',1.5)
            plot(xl,repmat(thr,size(xl)),'k:','linewidth',1.5)
            %settings
            tmp = [min(DFF),max(DFF)];
            set(gca,'ylim',tmp+[-0.1,0.2]*diff(tmp));
            linkaxes(ha,'x')
            set(gca,'xlim',xl,'box','on')
            
            %PLOT BOUTS & PEAKS
            Y = [min(ylim),max(ylim),max(ylim),min(ylim)];
            Z = zeros(1,4);                               
            for sta=1:noSTA
                [stage,label,col] = Stages{sta,:};
                res   = Bouts.(label);
                for bou = 1:numel(res.dur)
                    ind = res.ind(bou,:);
                    if fig==1
                        pks = res.peaks(bou);
                    else
                        pks = res.peaks2(bou);
                    end
                    %bouts
                    hp(1+sta) = fill(t(sort([ind,ind])),Y,Z,...
                        'edgecolor','none','facecolor',col,'facealpha',0.2);
                    legSTR{1+sta} = sprintf('%s',label);
                    uistack(hp(1+sta),'bottom')
                    text(mean(t(ind)),max(ylim),num2str(bou),...
                        'color',col,...
                        'horizontalalignment','center','verticalalignment','top')
                    %peaks
                    for k = 1:numel(pks.sta)
                        ind = pks.sta(k):pks.end(k);
                        plot(t(ind),DFF(ind),col);
                        hp(2+noSTA)     = plot(t(pks.ind(k)),pks.val(k),'ko');
                        plot(repmat(t(pks.ind(k)),1,2),...
                            repmat(pks.val(k),1,2)-[0,pks.height(k)],'k')
                        legSTR{2+noSTA} = 'Peaks';
                    end
                end %bout loop
            end %stage loop
            pos = get(gca,'position');
            ind = ~cellfun(@isempty,legSTR);
            legend(hp(ind),legSTR(ind))
            set(gca,'position',pos)
            zoom xon
        end
    end
    
    %% TABLE
    TABS = cell(2,2);
    for tab = 1:2
        %function as string
        clear Tab1 Tab2 labels1 labels2
        if tab==1
            funSTR = 'using only threshold';
            n = num2cell(1:tot.maxPKS);
        else
            funSTR = 'using threshold and findpeaks';
            n = num2cell(1:tot.maxPKS2);
        end
        %labels1
        tmp = [...
            cellfun(@(x)sprintf('Peak %i Value',x),n,'uniformoutput',false);...
            cellfun(@(x)sprintf('Peak %i Height',x),n,'uniformoutput',false);...
            cellfun(@(x)sprintf('Peak %i Area',x),n,'uniformoutput',false)];
        labels1 = ['Stage','Dur [s]','Mean DFF','Peaks N',tmp(:)',...
            'Total Mean DFF','Threshold','Threshold Function',...
            'Min Peak Prominence'];
        labels2 = {'Stage','Bout','Dur [s]','Mean DFF',...
            'Peak','Value','Height','Area',...
            'Total Mean DFF','Threshold','Threshold Function',...
            'Min Peak Prominence'};
        %table
        if tab==1
            Tab1 = cell(tot.boutsN,numel(labels1));
            Tab2 = cell(tot.peaksN,numel(labels2));
        else
            Tab1 = cell(tot.boutsN2,numel(labels1));
            Tab2 = cell(tot.peaksN2,numel(labels2));
        end        
        Tab1{1,strcmpi(labels1,'Total Mean DFF')}      = mea;
        Tab1{1,strcmpi(labels1,'Threshold')}           = thr;
        Tab1{1,strcmpi(labels1,'Threshold Function')}  = funSTR;
        Tab1{1,strcmpi(labels1,'Min Peak Prominence')} = minPeakProminence;
        Tab2{1,strcmpi(labels2,'Total Mean DFF')}      = mea;
        Tab2{1,strcmpi(labels2,'Threshold')}           = thr;
        Tab2{1,strcmpi(labels2,'Threshold Function')}  = funSTR;
        Tab2{1,strcmpi(labels2,'Min Peak Prominence')} = minPeakProminence;
        %init
        row1 = 0; %init
        row2 = 0;
        for sta = 1:noSTA
            labelSTA = Stages{sta,2};
            res = Bouts.(labelSTA);
            noBOU = numel(res.dur);
            for bou = 1:numel(res.dur)
                row1 = row1+1;
                Tab1{row1,strcmpi(labels1,'Stage')}    = labelSTA;
                Tab1{row1,strcmpi(labels1,'Dur [s]')}  = res.dur(bou);
                Tab1{row1,strcmpi(labels1,'Mean DFF')} = res.mean(bou);
                %peaks
                if tab==1
                    pks = res.peaks(bou);
                else
                    pks = res.peaks2(bou);
                end
                Tab1{row1,strcmpi(labels1,'Peaks N')}  = pks.N;
                for k = 1:pks.N
                    Tab1{row1,strcmpi(labels1,sprintf('Peak %i Value',k))} = ...
                        pks.val(k);
                    Tab1{row1,strcmpi(labels1,sprintf('Peak %i Height',k))} = ...
                        pks.height(k);
                    Tab1{row1,strcmpi(labels1,sprintf('Peak %i Area',k))} = ...
                        pks.area(k);
                    %table 2
                    row2 = row2+1;
                    Tab2{row2,strcmpi(labels2,'Stage')}    = labelSTA;
                    Tab2{row2,strcmpi(labels2,'Bout')}     = bou;
                    Tab2{row2,strcmpi(labels2,'Dur [s]')}  = res.dur(bou);
                    Tab2{row2,strcmpi(labels2,'Mean DFF')} = res.mean(bou);
                    Tab2{row2,strcmpi(labels2,'Peak')}     = k;
                    Tab2{row2,strcmpi(labels2,'Value')}    = pks.val(k);
                    Tab2{row2,strcmpi(labels2,'Height')}   = pks.height(k);
                    Tab2{row2,strcmpi(labels2,'Area')}     = pks.area(k);
                end
            end
        end
        TABS(tab,:) = {[labels1;Tab1],[labels2;Tab2]};
    end
    
    %% SAVE
    %table to save
    for tab = 1:2
        switch opt.table
            case 1
                Tab = TABS{tab,1};
            case 2
                Tab = TABS{tab,2};
            otherwise
                Tab = [];
        end
        %save
        if isempty(Tab)
            fprintf('%s Nothing Saved!\n',indent)
        else
            [sPath,sFile,sExt] = fileparts(saveFiles{tab});
            if isempty(sPath)
                sPath = rPath;
            elseif ~exist(sPath,'dir')
                mkdir(sPath)
            end
            sname = fullfile(sPath,[sFile,sExt]);
            if exist(sname,'file')==2
                delete(sname)
            end
            xlswrite(sname,Tab)
            try %if available
                xls_cellFit(sname)
            catch
            end
            fprintf('%s Saved: %s\n',indent,[sFile,sExt])
            
            %SAVE FIGURES
            if opt.plot
                [sPath,sFile] = fileparts(sname);
                saveas(HF(tab),fullfile(sPath,sFile))
            end
        end
    end
end %file loop