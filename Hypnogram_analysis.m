% analyzes sleep:
    % total duration [min & %]; number of events,
    % mean episode duration [min], latency to 1st WAKE/NREM/REM [min]
    % time course [min]
clc; clear; 
% close all;
addpath('Z:\Lab_resources\LabVIEW\functions')
addpath(genpath('Z:\OptoLab_v4.1\function'))

%PARAMETERS
%----------
%READ FILES (hypnogram, mat-files)
rFiles = '*Hypnogram.mat';

%stages in hypnogram to be analyzed
Stages = {... {number,label}
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    4,'cataplexy';...
    };
noSTA = size(Stages,1);

%% PARAMETERS TO CHANGE
minduration = ;    % [s] for latency
neededuration = ;   

plotting = true; %false for testing
saving   = true; %false for testing


%% MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%FILE-LIST
rFiles = selectionList((findFiles(rFiles)),'HypnoAnalysis.mat');

%FILE LOOP
noFIL = numel(rFiles);
if noFIL==0
    fprintf(2,'No File Selected!\n')
    return
end
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);

for fil = 1:noFIL
    close all;
    rname = rFiles{fil};
    [rPath,rFile,rExt] = fileparts(rFiles{fullfile(fil)});
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath)
    slashIdx = strfind(rname, '\');
    dataname = rname(slashIdx(4)+1:slashIdx(7)-1);
    group = rname(slashIdx(3)+1:slashIdx(4)-1);
    
    %LOAD Hypnogram.mat file
    tmp = load(rname);
    Hyp = tmp.Hypnogram;
    fs = tmp.fs;
    Time = (1:numel(Hyp)/fs)/60; %[min]
    
    %TITLE for figure
    [~,id] = fileparts(rPath);
    id = strrep(id,'_','\_');
    fprintf('%s Load Data\n',indent)
    
    %Total duration
    duration = length(Hyp)/fs; %[s]
    duration2 = duration/3600; %[h]
    if duration2 == neededuration
        fprintf('%s Duration = %i h \n',indent,neededuration)
    elseif duration2 > neededuration 
        fprintf ('%s Longer than %i h! \n',indent,neededuration)
    elseif duration2 < neededuration
        fprintf(2, '%s Shorter than %i h! \n',indent,neededuration)
    end
    
    latency = NaN(1,noSTA);
    tc = NaN(neededuration-1,noSTA);
    
    for sta = 1:noSTA
        [stage,label] = Stages{sta,:};
        
        %Total duration + percentage
        dur(sta) = sum(Hyp==stage)/fs/60; %[min]
        perc(sta) = (dur(sta)/duration*100);
        
        %Find episode start-end
        tmp = [NaN;Hyp(:);NaN];
        indSTA = find(tmp(1:end-1)~=stage & tmp(2:end)==stage);
        indEND = find(tmp(1:end-1)==stage & tmp(2:end)~=stage)-1;
        
        %Episode number + mean duration
        episnr(sta) = length(indSTA);
        epis = (indEND-indSTA+1)/fs; %[s]
        if isempty(epis)
            epd(stage) = 0;
            epdstd(stage) = 0;
            tc(:,stage) = 0;
            tcstd(:,stage) = 0;
            fprintf('[\b%s No %s!]\b\n',indent,label)
        else
            epd(stage) = mean(epis); %[min]
            epdstd(stage) = std(epis);
        end
    end
    
    %Latency info
    if latency(2)<latency(1)
        fprintf(2,'%s !!! NREM (%.2f) before Wake (%.2f)!!!\n',indent,latency(2),latency(1))
    elseif latency(3)<latency(1)
        fprintf(2,'%s !!! REM (%.2f) before Wake (%.2f)!!!\n',indent,latency(3),latency(1))
    elseif latency(3)<latency(2)
        fprintf(2,'%s !!! REM (%.2f) before NREM (%.2f)!!!\n',indent,latency(3),latency(2))
    end
    
    
    %% FIGURE
    if plotting
        hf(1) = figure('WindowState','maximize');
        
        ha(1) = subplot(2,2,1);
        hold on
        bar(1,dur(1),'b')
        bar(2,dur(2),'r')
        bar(3,dur(3),'g')
        bar(4,dur(4),'y')
        title(sprintf('Total duration [min]',id));
        set(gca,'xtick',[1:4],'xticklabels',{'Wake','NREM','REM', 'cataplexy'})
        
        ha(2) = subplot(2,2,2);
        pie(perc)
        colormap([0 0 1  ; 1 0 0  ; 0 1 0 ; 1 1 0]);
        legend(Stages(5:8),'Location','eastoutside');
        title(sprintf('Total percentage',id));
        
        ha(3) = subplot(2,2,3);
        hold on
        bar(1,episnr(1),'b')
        bar(2,episnr(2),'r')
        bar(3,episnr(3),'g')
        bar(4,episnr(4),'y')
        title(sprintf('Number of episodes',id));
        set(gca,'xtick',[1:4],'xticklabels',{'Wake','NREM','REM','cataplexy'})

        ha(4) = subplot(2,2,4);
        hold on
        bar(1,epd(1),'b')
        bar(2,epd(2),'r')
        bar(3,epd(3),'g')
        bar(4,epd(4),'y')
        errorbar(1:4,epd,NaN(size(epd)),epdstd,'linestyle','none','color','k');
        title(sprintf('Mean episode duration [min]',id));
        set(gca,'xtick',[1:4],'xticklabels',{'Wake','NREM','REM','cataplexy'}) 
                
        sgtitle(sprintf('%s \n duration %i h',id,duration2),'FontWeight','bold','FontSize',20)
        
    end
    %% SAVE
    if saving
        %Figure
        snamef = fullfile(rPath,'HypnoAnalysis');
        if exist([snamef,'.png'],'file')==2
            delete([snamef,'.png'])
        end
        print(hf,snamef,'-dpng','-r300')
        fprintf('%s File PNG Saved: HypnoAnalysis.png\n',indent)
        
        if exist([snamef,'.eps'],'file')==2
            delete([snamef,'.eps'])
        end
        print(hf,snamef,'-depsc','-r300')
        fprintf('%s File EPS Saved: HypnoAnalysis.eps\n',indent)
        
        %Mat file with all variables
        snamem = fullfile(rPath,'HypnoAnalysis');
        
        HypnoAnalysis.Stages = {'Wake','NREM','REM','cataplexy'};
        
        HypnoAnalysis.Duration.TotalSeconds = duration;
        HypnoAnalysis.Duration.Minutes = dur;
        
        HypnoAnalysis.Percentage = perc;
        
        HypnoAnalysis.Episodes.Number = episnr;
        HypnoAnalysis.Episodes.MeanDurationSec = epd;
        HypnoAnalysis.Episodes.Std = epdstd;
        
        HypnoAnalysis.LatencyMinutes = latency;
        
        HypnoAnalysis.TimeCourse.Time = 1:(neededuration);
        HypnoAnalysis.TimeCourse.Minutes = tc;
        
        
        save(snamem,'-struct','HypnoAnalysis')
        fprintf('%s File MAT Saved: HypnoAnalysis2.mat \n\n',indent)
        
    else
        fprintf(2,' Nothing Saved!\n')
    end
    
end