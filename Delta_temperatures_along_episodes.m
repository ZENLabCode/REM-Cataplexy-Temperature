clc; clear; close all; fclose all;
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
addpath('\function\excel')
addpath('\function\misc') %selectionList
addpath('\function\figures') %fig_size


%% START PATH AND SELECT FILES

startPath = '';
sumF = fullfile(startPath,'**\*summary.xlsx');
tmp = dir(sumF);
sumf = fullfile({tmp.folder},{tmp.name});
tmp = cellfun(@(x)strrep(x,fileparts(fileparts(fileparts(x))),''),...
    sumf,'uniformoutput',false);
% for every (filesparts...) the list goes back
% of 1 folder, then it substitutes the previous folders
% with a blank. Used to make the list easier to read.
[~,indx] = selectionList(tmp,[],struct('unique',false,'sort',false));



%% Load the Excel files
for an = 1:numel(indx)
    ok = true;
    [nbr] = indx(an);
    sumfl = sumf{nbr};
    [~,rFile] = fileparts(sumfl);
    tmp = regexp(rFile,' ','split');
    name = strjoin(tmp(1:end-1),' ');
    rPath    = fileparts(sumfl);
    data = readtable(sumfl);

    % Extract columns
    time = datetime(data.ClockHour, 'ConvertFrom', 'datenum');                % Time in seconds
    sleepStage = data.SleepStage;                                             % Sleep stage column
    tC = data.CoreT_;                                                         % Core Temperature column
    tA = data.AmbientT_;                                                      % Ambient Temperature column
    tS = data.SkinT_;                                                         % Skin Temperature column
    tB = data.BrainT_;                                                        % Brain Temperature column
    vars = {'tC', 'tA', 'tS', 'tB'};
    for sm = 1:numel(vars)
        if all([eval(vars{sm}) == 0] | isnan(eval(vars{sm})))
            eval([vars{sm}, '(:) = NaN;']);
        end
    end

    % Check the max temperature to decide if Baseline or Warming/Cooling
    split_into_groups = max(tA) > 28;                                         % Warming and Cooling if max temperature > 28°C, otherwise Baseline

    % Convert time to hours
    hrs = floor(hour(time)) - floor(hour(time(1))) + 1;


    %% Define bin edges
    bins = 40:5:200;
    minLength = min(bins);
    maxLength = max(bins);
    % Original 10s steps: 20s, 40s, 60s, etc.until 250s. Set here the desired binning
    disp(['Seconds bins: ', num2str(bins)]);
    binl = bins(1:end-1);

    % Initialize cell arrays for storing results
    h1R = ['Bin_Seconds'; num2cell(binl(:))];
    h1T = ['Bout seconds'; num2cell(1:(maxLength))'];
    if split_into_groups
        h2R = ['Bin_Seconds'; num2cell(binl(:))];
        h2T = ['Bout seconds'; num2cell(1:(maxLength))'];
    end
0





    %% Get unique sleep stages
    ustages = unique(sleepStage(~isnan(sleepStage)));
    stageLbs = {'Wake', 'Nrem', 'REM', 'Cataplexy'};

    % Loop over each sleep stage and analyze separately
    for stage_idx = 1:length(ustages)
        stage = ustages(stage_idx);
        stagename = stageLbs{stage_idx};

        % Mask for current sleep stage
        isstage = (sleepStage == stage);

        % Initialize episode variables
        boutdur = [];
        stagetA = [];
        stagetC = [];
        stagetB = [];
        stagetS = [];
        boutstart = [];
        currlength = 0;
        tAval = [];
        tCval = [];
        tBval = [];
        tSval = [];
        starttime = NaN;

        for i = 1:length(sleepStage)
            if isstage(i)
                if currlength == 0
                    starttime = time(i);                                      % Save start time of episode
                    idxStart = i;
                end
                currlength = currlength + 1;
                tAval = [tAval; tA(i)];
                tCval = [tCval; tC(i)];
                tBval = [tBval; tB(i)];
                tSval = [tSval; tS(i)];

            else
                if currlength > 0
                    % Save episode data
                    idxEnd = i;
                    idxEnd = min(idxStart + (maxLength-1), idxEnd);
                    boutdur = [boutdur; currlength];
                    stagetA = [stagetA; mean(tAval)];
                    stagetC = [stagetC; mean(tCval)];
                    stagetB = [stagetB; mean(tBval)];
                    stagetS = [stagetS; mean(tSval)];
                    boutstart = [boutstart; starttime];
                    idxSE(size(boutdur, 1),:) = [idxStart,idxEnd];


                    % Reset variables
                    currlength = 0;
                    tAval = [];
                    tCval = [];
                    tBval = [];
                    tSval = [];
                end
            end
        end

        % If last row was part of an episode, save it
        if currlength > 0
            boutdur = [boutdur; currlength];
            stagetA = [stagetA; mean(tAval)];
            stagetC = [stagetC; mean(tCval)];
            stagetB = [stagetB; mean(tBval)];
            stagetS = [stagetS; mean(tSval)];
            boutstart = [boutstart; starttime];
            if i == 28800
                idxEnd = 28800;
                idxEnd = min(idxStart + (maxLength -1), idxEnd);
                idxSE(size(boutdur, 1),:) = [idxStart,idxEnd];
            end
        end

        % Convert episode start times to hours
        ephours =  floor(hour(boutstart)) - floor(hour(boutstart(1))) + 1;

        % Assign episodes to Warm or Cool based on hour
        if split_into_groups
            h1episodes = ismember(ephours, [1, 3, 5, 7]) | ~split_into_groups;
            h2episodes = ismember(ephours, [2, 4, 6, 8]);
        else
            h1episodes = true(size(ephours));
        end



        alTa = NaN(size(boutdur,1), maxLength);
        alTs = NaN(size(boutdur,1), maxLength);
        alTc = NaN(size(boutdur,1), maxLength);
        alTb = NaN(size(boutdur,1), maxLength);
        for d = 1:size(boutdur,1)
            %fprintf(['aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa' ...
            %'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'])
            alTa(d,1:(idxSE(d,2)-idxSE(d,1))) = tA(idxSE(d,1):idxSE(d,2)-1,1)';
            alTs(d,1:(idxSE(d,2)-idxSE(d,1))) = tS(idxSE(d,1):idxSE(d,2)-1,1)';
            alTc(d,1:(idxSE(d,2)-idxSE(d,1))) = tC(idxSE(d,1):idxSE(d,2)-1,1)';
            alTb(d,1:(idxSE(d,2)-idxSE(d,1))) = tB(idxSE(d,1):idxSE(d,2)-1,1)';
        end

        alTa = alTa - alTa(:,1);                                               %difference from episode start
        alTs = alTs - alTs(:,1);
        alTc = alTc - alTc(:,1);
        alTb = alTb - alTb(:,1);


        meanTa = mean(alTa(h1episodes & sum(~isnan(alTa), 2) >= 40, :),...
            1, "omitnan");                                                       %average of all episodes longer than minimum set
        stdevTa = std(alTa(h1episodes & sum(~isnan(alTa), 2) >= 40, :),...
            1, "omitnan");
        meanTs = mean(alTs(h1episodes & sum(~isnan(alTs), 2) >= 40, :),...
            1, "omitnan");
        stdevTs = std(alTs(h1episodes & sum(~isnan(alTs), 2) >= 40, :),...
            1, "omitnan");
        meanTc = mean(alTc(h1episodes & sum(~isnan(alTc), 2) >= 40, :),...
            1, "omitnan");
        stdevTc = std(alTc(h1episodes & sum(~isnan(alTc), 2) >= 40, :),...
            1, "omitnan");
        meanTb = mean(alTb(h1episodes & sum(~isnan(alTb), 2) >= 40, :),...
            1, "omitnan");
        stdevTb = std(alTb(h1episodes & sum(~isnan(alTb), 2) >= 40, :),...
            1, "omitnan");

        if split_into_groups
            meanTa2 = mean(alTa(h2episodes & sum(~isnan(alTa), 2) >= 40, :),...
                1, "omitnan");
            stdevTa2 = std(alTa(h2episodes & sum(~isnan(alTa), 2) >= 40, :),...
                1, "omitnan");
            meanTs2 = mean(alTs(h2episodes & sum(~isnan(alTs), 2) >= 40, :),...
                1, "omitnan");
            stdevTs2 = std(alTs(h2episodes & sum(~isnan(alTs), 2) >= 40, :),...
                1, "omitnan");
            meanTc2 = mean(alTc(h2episodes & sum(~isnan(alTc), 2) >= 40, :),...
                1, "omitnan");
            stdevTc2 = std(alTc(h2episodes & sum(~isnan(alTc), 2) >= 40, :),...
                1, "omitnan");
            meanTb2 = mean(alTb(h2episodes & sum(~isnan(alTb), 2) >= 40, :),...
                1, "omitnan");
            stdevTb2 = std(alTb(h2episodes & sum(~isnan(alTb), 2) >= 40, :),...
                1, "omitnan");
        end

        % Function to bin temperature data
        clear ind1 ind2
        fun  =  @(ind,data) accumarray(ind(ind>0), data(ind>0), ...
            [numel(bins)-1,1], @mean, NaN);
        stagetemp = [stagetA, stagetB, stagetC, stagetS];
        %Binning sleep data
        [h1binning,~,ind1] = histcounts(boutdur(h1episodes), bins);
        h1R = [h1R, [['Episodes number ', stagename]; num2cell(h1binning(:))]];

        if split_into_groups
            [h2binning,~,ind2] = histcounts(boutdur(h2episodes), bins);
            h2R = [h2R, [['Episodes number ', stagename]; ...
                num2cell(h2binning(:))]];

        end
        %Binning temperature data. PS: ind==0 are elements not in any bins
        if numel(ind1)<=1 & max(ind1)<1
            fprintf(2, 'Something is not right in %s (stage: %s)\n', name, stagename);
            ok = false;
            continue
        end
        tempLbs = {'Ambient T. ','Brain T. ','Core T. ','Skin T. '};
        for t = 1:4
            st =  stagetemp(:, t);

            h1avgt = {meanTa, meanTb, meanTc, meanTs};

            h1tavgtemp = h1avgt{t};
            tempname = tempLbs{t};

            h1avgtemp = fun(ind1, st(h1episodes));
            %display [bins center, histogram count, mean of accumarray]
            %disp([1:15; diff(bins)/2+bins(1:end-1); h1binning; h1avgtemp'])
            h1R = [h1R, [['Average ', tempname, stagename];...
                num2cell(h1avgtemp(:))]];
            h1T = [h1T, [['Mean ', tempname, stagename];...
                num2cell(h1tavgtemp(:))]];
            if split_into_groups
                h2avgt = {meanTa2, meanTb2, meanTc2, meanTs2};
                h2tavgtemp = h2avgt{t};
                if any(ind2) > 0
                    if ~isempty(ind2) == 1
                        h2avgtemp = fun(ind2, st(h2episodes));
                        h2R = [h2R, [['Average ', tempname, stagename]; ...
                            num2cell(h2avgtemp(:))]];
                        h2T = [h2T, [['Mean ', tempname, stagename];...
                            num2cell(h2tavgtemp(:))]];
                    else
                    end
                else
                end
            end
        end
        wrmm =  'baseline';
        c = 'black';
        if split_into_groups
            wrmm = 'warm';
            c = 'r';
        else
        end
        figure
        plot(alTc(h1episodes & sum(~isnan(alTc), 2) >= 40, :)')
        hold on
        plot(meanTc, c,'LineWidth', 5);
        title([stagename, ' Tc ', wrmm]);
        hold off
        figure
        plot(alTs(h1episodes & sum(~isnan(alTs), 2) >= 40, :)')
        hold on
        plot(meanTs, c,'LineWidth', 5);
        title([stagename, ' Ts ', wrmm]);
        hold off
        figure
        plot(alTb(h1episodes & sum(~isnan(alTb), 2) >= 40, :)')
        hold on
        plot(meanTb, c,'LineWidth', 5);
        title([stagename, ' Tb ', wrmm]);
        hold off
        figure
        plot(alTa(h1episodes & sum(~isnan(alTa), 2) >= 40, :)')
        hold on
        plot(meanTa, c,'LineWidth', 5);
        title([stagename, ' Ta ', wrmm]);
        hold off
        if split_into_groups
            figure
            plot(alTc(h2episodes & sum(~isnan(alTc), 2) >= 40, :)')
            hold on
            plot(meanTc2, 'b', 'LineWidth', 5);
            title([stagename, ' Tc cool']);
            hold off
            figure
            plot(alTs(h2episodes & sum(~isnan(alTs), 2) >= 40, :)')
            hold on
            plot(meanTs2, 'b','LineWidth', 5);
            title([stagename, ' Ts cool']);
            hold off
            figure
            plot(alTb(h2episodes & sum(~isnan(alTb), 2) >= 40, :)')
            hold on
            plot(meanTb2, 'b','LineWidth', 5);
            title([stagename, ' Tb cool']);
            hold off
            figure
            plot(alTa(h2episodes & sum(~isnan(alTa), 2) >= 40, :)')
            hold on
            plot(meanTa2, 'b','LineWidth', 5);
            title([stagename, ' Ta cool']);
            hold off
        else
        end
        % Saving plots
        outF_mat = fullfile(rPath, 'delta plots mat');
        outF_jpg = fullfile(rPath, 'delta plots jpg');

        if ~exist(outF_mat, 'dir')
            mkdir(outF_mat);
        end
        if ~exist(outF_jpg, 'dir')
            mkdir(outF_jpg);
        end

        figs = findall(0, 'Type', 'figure');

        for u = 1:length(figs)
            fig = figs(u);

            % Set to full screen
            figure(fig);
            set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % Maximize
            filenm = fullfile(outF_mat, sprintf('%s_%s', name, get(get(gca, 'Title'), 'String')));
            filenj = fullfile(outF_jpg, sprintf('%s_%s', name, get(get(gca, 'Title'), 'String')));

            % Save as .fig and .jpg
            savefig(fig, [filenm, '.fig']);
            exportgraphics(gcf, [filenj, '.jpg'], 'Resolution', 300);
        end
        close all
    end

    %% Save results to an Excel file with appropriate name and sheets
    ok = true;
    if ok
        rPatht = '';
        sFile = fullfile(rPatht, sprintf('%s Bout Temp.xlsx',name));
        if split_into_groups
            writecell(h1R, sFile, 'Sheet', 'Warming bouts');
            writecell(h1T, sFile, 'Sheet', 'Warming means');
            writecell(h2R, sFile, 'Sheet', 'Cooling bouts');
            writecell(h2T, sFile, 'Sheet', 'Cooling means');
        else
            writecell(h1R, sFile, 'Sheet', 'Baseline bouts');
            writecell(h1T, sFile, 'Sheet', 'Baseline means');
        end

        if split_into_groups
            disp('Ta Challenge data');
        else
            disp('Baseline data');
        end
        disp(['Results saved to ', sFile]);
    else
    end
end
