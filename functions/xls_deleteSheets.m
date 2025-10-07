% xls_deleteSheets
%-------------------------------------------------------------------------
% Delete specific sheets from xls-file.
%
% SYNTAX:
%   xls_deleteSheets(filename,sheetNames,option)
%       - filename is the full xls-file name
%       - sheetNames is a cell-array with sheet names
%         (or char for single sheet)
%       - option is 1 or -1
%           if 1 , deletes the sheets in sheetNames
%           if -1, deletes the sheets NOT in sheetNames
%       NOTE: - filename and sheetNames are case-insensitive.
%             - alse works if the file is open. Changes will be
%               automatically saved!
%
%
% Thomas Rusterholz, Dez 1 2011
%-------------------------------------------------------------------------

function xls_deleteSheets(filename,sheetNames,option)

%INPUT CHECK
if nargin~=3
   error('Function must have 3 input variables\n') 
end
if isempty(fileparts(filename)) %full file name
    filename=fullfile(pwd,filename);
end
if ~exist(filename,'file')
    warning('File %s does not exist\n',filename)
    return
end
if ischar(sheetNames)
    sheetNames={sheetNames};
end
if ~ismember(option,[-1,1])
    error('3rd input variable must be 1 or -1\n')
end

%EXCEL/WORKBOOK, access|start/activate|open
try
    %handle to running server
    Excel=actxGetRunningServer('Excel.Application');
    serverWasRunning=true;
    %activate workbook if open
    for wb=1:Excel.Workbooks.Count
        workbookWasOpen=strcmpi(filename,...
            Excel.Workbooks.Item(wb).fullName);
        if workbookWasOpen
            Excel.Workbooks.Item(wb).Activate;
            break
        end
    end
    %open workbook if closed
    if ~workbookWasOpen
        Excel.Workbooks.Open(filename);
    end
%if server is not running
catch
    serverWasRunning=false;
    workbookWasOpen=false;
    Excel=actxserver('Excel.Application'); %start server
    Excel.Workbooks.Open(filename);        %open workbook
end

%DELETE SHEETS
Excel.DisplayAlerts=false;
noSHE=Excel.Sheets.count; %number of existing sheets
if option==1
    for she=noSHE:-1:1
        if ismember(lower(Excel.Sheets.Item(she).Name),sheetNames)
            Excel.Sheets.Item(she).Delete;
        end
    end
else
    for she=noSHE:-1:1
        if ~ismember(lower(Excel.Sheets.Item(she).Name),lower(sheetNames))
            Excel.Sheets.Item(she).Delete;
        end
    end
end
Excel.DisplayAlerts=true;
Excel.ActiveWorkbook.Save;

%close workbook/server if it was not open/running
if serverWasRunning
    if ~workbookWasOpen
        Excel.ActiveWorkbook.Close;
    end
else
    Excel.Quit;
end