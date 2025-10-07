clc; clear; close all; fclose all;
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
addpath('\function\excel')
addpath('\function\misc') %selectionList
addpath('\function\figures') %fig_size

%% START PATH AND SELECT FILES

startPath = '';
sumF = fullfile(startPath,'**\*Bout Temp.xlsx');
tmp = dir(sumF);
sumf = fullfile({tmp.folder},{tmp.name});
tmp = cellfun(@(x)strrep(x,fileparts(fileparts(fileparts(x))),''),...
    sumf,'uniformoutput',false);
[~,indx] = selectionList(tmp,[],struct('unique',false,'sort',false));

%% Initialize and sort
condLbs = {'Baseline', 'Warming' , 'Cooling'};
rowLbs = {'Wake', 'Nrem', 'REM', 'Cataplexy'};
colLbs = {'Core', 'Skin', 'Brain', 'Amb'};
mavgb = cell(length(rowLbs)+1, length(colLbs)+1);
mavgb(1,2:end) = colLbs;
mavgb(2:end,1) = rowLbs;                                                      %Cell array to store data for mouse
for cn = 1:numel(condLbs)
    mouse.(condLbs{cn}) = mavgb;
    final.(condLbs{cn}) = mavgb;
    mSEM.(condLbs{cn}) = mavgb;
    fSEM.(condLbs{cn}) = mavgb;
end


for an = 1:numel(indx)
    ok = true;
    if exist("mname",'var')
        rowLabels = ['Bout seconds'; num2cell((1:numel(data(:,1)))')];
        colLabels = {};
        for r = 1:numel(rowLbs)
            for c = 1:numel(colLbs)
                colLabels{end+1} = ['Mean ' rowLbs{r} ' ' colLbs{c}];
            end
        end
        exT = cell(length(rowLabels), length(colLabels)+1);
        exT(:,1) = rowLabels;
        exT(1,2:end) = colLabels;
    else
    end
    [nbr] = indx(an);
    sumfl = sumf{nbr};
    [~,rFile] = fileparts(sumfl);
    tmp = regexp(rFile,' ','split');
    name = strjoin(tmp(1:end-1),' ');
    rPath    = fileparts(sumfl);
    sheetNames= sheetnames(sumfl);
    check = contains(sheetNames(:,1), 'Warming', 'IgnoreCase', true);
    split = any(check);                                                     %Variable to distinguish Baseline and Ta Challenge


    %% Extract data columns
    if exist("mname",'var')
        mname2 = char(tmp(1));
        if ~strcmp(mname, mname2)                                            %Same mouse?
            sFile = fullfile(rPatht, sprintf...
                ('%s Bout Temp Average.xlsx',mname));
            for cn = 1:numel(condLbs)
                T = mouse.(condLbs{cn});
                Tsem = mouse.(condLbs{cn});
                assert(all(all(cellfun(@isnumeric,T(2:end, 2:end)))), 'Not numeric data')

                % Apply mean skipping labels (row 1 and col 1)
                Tsem(2:end, 2:end) = cellfun(@(x) std(x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(x), 2)), ...
                    Tsem(2:end, 2:end), 'UniformOutput', false);
                T(2:end, 2:end) = cellfun(@(x) mean(x, 2, 'omitnan') ,...
                    T(2:end, 2:end), 'UniformOutput', false);

                mouse.(condLbs{cn}) = T;
                mSEM.(condLbs{cn}) = Tsem;
                % mouse.(condLbs{cn}){rowIdx, colIdx} = mean(mouse.(condLbs{cn}){rowIdx, colIdx}, 2, 'omitnan');

                x = 1;
                for rr = 2:numel(rowLbs)+1
                    for cc = 2:numel(colLbs)+1
                        x = x + 1;
                        if ~isempty(mouse.(condLbs{cn}){rr,cc})
                            exT(2:end, x) = num2cell(mouse.(condLbs{cn}){rr,cc});
                            final.(condLbs{cn}){rr,cc} = [final.(condLbs{cn}){rr,cc}, mouse.(condLbs{cn}){rr,cc}];
                            fSEM.(condLbs{cn}){rr,cc} = [fSEM.(condLbs{cn}){rr,cc}, mSEM.(condLbs{cn}){rr,cc}];
                        else
                        end
                        
                    end
                end
                T2 = exT;
                writecell(T2, sFile, 'Sheet', [condLbs{cn},' means']);  %xlswrite
                mouse.(condLbs{cn}) = mavgb;
                mSEM.(condLbs{cn}) = mavgb;
            end
            
        else
        end
    else
    end
    if split
        data = readtable(sumfl, 'Sheet', 2);
        data2 = readtable(sumfl, 'Sheet', 4);
    else
        data = readtable(sumfl, 'Sheet', 2);
    end
    %     mAvg = {'Bout seconds'; num2cell(1:(numel(data(:,1))))'};
    %     tAvg= ['Bout seconds'; num2cell(1:(numel(data(:,1))))'];
    if any(contains(data.Properties.VariableNames, 'Cataplexy', 'IgnoreCase', true))
        t = 0;
    else
        t = 1;
    end
    for tt = 1:numel(rowLbs)-t
        stg = char(rowLbs(tt));
        for z = 1:numel(colLbs)
            tempr = char(colLbs(z));

            if split
                cn = 2;
                ind = contains(data.Properties.VariableNames, tempr)...
                    & contains(data.Properties.VariableNames, stg);
                if any(ind)
                    recW.(stg).(tempr) = data.(data.Properties.VariableNames...
                        {ind});
                else
                    recW.(stg).(tempr) = [];
                end
                %                 recW.(stg).(tempr) = data.(data.Properties.VariableNames...
                %                     {contains(data.Properties.VariableNames, tempr)...
                %                     & contains(data.Properties.VariableNames, stg)});

                ind = contains(data2.Properties.VariableNames, tempr)...
                    & contains(data2.Properties.VariableNames, stg);
                if any(ind)
                    recC.(stg).(tempr) = data2.(data2.Properties.VariableNames...
                        {ind});
                else
                    recC.(stg).(tempr) = [];
                end
                %                 recC.(stg).(tempr) = data2.(data2.Properties.VariableNames...
                %                     {contains(data2.Properties.VariableNames, tempr)...
                %                     & contains(data2.Properties.VariableNames, stg)});

            else
                cn = 1;
                ind = contains(data.Properties.VariableNames, tempr)...
                    & contains(data.Properties.VariableNames, stg);
                if any(ind)
                    recB.(stg).(tempr) = data.(data.Properties.VariableNames...
                        {ind});
                else
                    recB.(stg).(tempr) = [];
                end

                %                 recB.(stg).(tempr) = data.(data.Properties.VariableNames...
                %                     {contains(data.Properties.VariableNames, tempr)...
                %                     & contains(data.Properties.VariableNames, stg)});
            end


            % Find the row and column index
            rowIdx = find(strcmp(mouse.(condLbs{cn})(2:end,1), stg)) + 1;
            colIdx = find(strcmp(mouse.(condLbs{cn})(1,2:end), tempr)) + 1;

            % Assign value for mouse average
            if split
                mouse.(condLbs{cn}){rowIdx, colIdx} = [mouse.(condLbs{cn}){rowIdx, colIdx}, recW.(stg).(tempr)];
                mouse.(condLbs{cn+1}){rowIdx, colIdx} = [mouse.(condLbs{cn+1}){rowIdx, colIdx}, recC.(stg).(tempr)];
            else
                mouse.(condLbs{cn}){rowIdx, colIdx} = [mouse.(condLbs{cn}){rowIdx, colIdx}, recB.(stg).(tempr)];



                %     plot(numel(tS), tS(:,1), 'b', 'LineWidth', 3);                                          %sum(~isnan(tS))
                %     legSTR = {'Baseline'};
                %     xlabel('\fontsize{25}Time (s)');
                %     ylabel('\fontsize{25}Temp. (°C)');

            end
        end
    end
    mname = char(tmp(1));
    rPath2 = rPath;
    rPatht = '';
    if an == numel(indx)                                                   %Last mouse?
        sFile = fullfile(rPatht, sprintf...
            ('%s Bout Temp Average.xlsx',mname));
        for cn = 1:numel(condLbs)
            T = mouse.(condLbs{cn});
            Tsem = mouse.(condLbs{cn});
            assert(all(all(cellfun(@isnumeric,T(2:end, 2:end)))), 'Not numeric data')

            % Apply mean skipping labels (row 1 and col 1)
            Tsem(2:end, 2:end) = cellfun(@(x) std(x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(x), 2)), ...
                Tsem(2:end, 2:end), 'UniformOutput', false);

            T(2:end, 2:end) = cellfun(@(x) mean(x, 2, 'omitnan') ,...
                T(2:end, 2:end), 'UniformOutput', false);

            mouse.(condLbs{cn}) = T;
            mSEM.(condLbs{cn}) = Tsem;
            % mouse.(condLbs{cn}){rowIdx, colIdx} = mean(mouse.(condLbs{cn}){rowIdx, colIdx}, 2, 'omitnan');

            x = 1;
            for rr = 2:numel(rowLbs)+1
                for cc = 2:numel(colLbs)+1
                    x = x + 1;
                    if ~isempty(mouse.(condLbs{cn}){rr,cc})
                        exT(2:end, x) = num2cell(mouse.(condLbs{cn}){rr,cc});
                        final.(condLbs{cn}){rr,cc} = [final.(condLbs{cn}){rr,cc}, mouse.(condLbs{cn}){rr,cc}];
                    else
                    end
                end
            end
            T2 = exT;
            writecell(T2, sFile, 'Sheet', [condLbs{cn},' means']);
            
            mouse.(condLbs{cn}) = mavgb;
            mSEM.(condLbs{cn}) = mavgb;
        end
    else
    end
end


exT = cell(length(rowLabels), length(colLabels)+1);
exT(:,1) = rowLabels;
exT(1,2:end) = colLabels;
[filename, sPatht] = uiputfile('*.xlsx', 'Save As',...
    sprintf('Bout Temp Average final.xlsx'));
if isequal(filename, 0)
    disp('nope');
    return;
end
sFile = fullfile(sPatht, filename);
for cn = 1:numel(condLbs)
    T3 = final.(condLbs{cn});
    T3sem = final.(condLbs{cn});
    assert(all(all(cellfun(@isnumeric,T3(2:end, 2:end)))), 'Not numeric data')

    % Apply mean skipping labels (row 1 and col 1)
     T3sem(2:end, 2:end) = cellfun(@(x) std(x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(x), 2)), ...
                T3sem(2:end, 2:end), 'UniformOutput', false);
    T3(2:end, 2:end) = cellfun(@(x) mean(x, 2, 'omitnan') ,...
        T3(2:end, 2:end), 'UniformOutput', false);

    final.(condLbs{cn}) = T3;
    fSEM.(condLbs{cn}) = T3sem;
    % mouse.(condLbs{cn}){rowIdx, colIdx} = mean(mouse.(condLbs{cn}){rowIdx, colIdx}, 2, 'omitnan');

    x = 1;
    for rr = 2:numel(rowLbs)+1
        for cc = 2:numel(colLbs)+1
            x = x + 1;
            if ~isempty(final.(condLbs{cn}){rr,cc})
                exT(2:end, x) = num2cell(final.(condLbs{cn}){rr,cc});
            else
            end
        end
    end
    T4 = exT;
    writecell(T4, sFile, 'Sheet', [condLbs{cn},' means']);
end
col = {'black', 'r' , 'b'};
for rr = 2:numel(rowLbs)+1
        for cc = 2:numel(colLbs)+1 
for cn = 1:numel(condLbs)
time = 1:numel(final.(condLbs{cn}){rr,cc});
hf = figure;
hf.WindowState = 'maximized';
plot(time, final.(condLbs{cn}){rr,cc}, col{cn}, 'LineWidth', 3);
legSTR = colLbs(cn);
xlabel('\fontsize{25}Time')
ylabel('\fontsize{25}Temp. (°C)')
title(sprintf('\\fontsize{30}%s',(rowLbs(rr-1)),colLbs(cc-1)));
%SAVE PLOTS
oFol = [rPatht, '\plots']; % or wherever you want.
plotname = sprintf('%s', name)
fullplotname = fullfile(oFol, plotname);
%freeze paper size
set(gcf,'unit','inches','paperUnits','inches'); pause(0.5)
pos = get(gcf,'position');
set(gcf,'paperSize',pos(3:4),'paperPosition',[0 0 pos(3:4)]);
print([fullplotname,'.png'],'-dpng','-r300');
saveas(gcf,[fullplotname,'.fig'])
end
hold on
        end
end

