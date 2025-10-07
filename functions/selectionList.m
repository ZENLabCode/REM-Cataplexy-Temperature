% selectionList
%-------------------------------------------------------------------------
% GUI for to lists data strings and return selected ones. List will be
% displayed alphabetically & unique! Basically made for to select
% automatically files found (for a re-check).
%
% Optionally (only if list is a path-list) script can check whether some
% specific files exist in given paths. It will split the list into two
% lists, one where files exist and one where files not exist. This makes it
% select paths which needs to be analyzed, and/or re-analyzed.
% NOTE: in this case, in-existing paths will be removed from list by
%       default.
%
%
% SYNTAX
%   list = selectionList(list)
%       Displays list (cell), unique and sorted with natsort (if available)
%       Returns selected list as cell.
%   list = selectionList(paths,file)
%       Displays two lists. Splits list into a list where file not exist
%       (in the same path) and a list where file exist.
%       Can be used for to check whether a script was already run or not,
%       by using file as the output file of the script.
%       NOTE: - file sholud have no path, as it will use the paths from
%               input list.
%             - Wildcard allowed, e.g. file = 'amp*.dat'
%               (check if any 'amp*.dat' exist)
%             - ignores input file if empty or no char
%       Returns a list with selections from both lists
%    list = selectionList(paths,file,opt) ... (optional)
%       opt is a structure
%       opt.title  : another titel than 'Select Files', 'Select Paths',...
%       opt.unique : unique paths/files (default is true)
%       opt.sort   : natsort paths/files (default is true) 
%       
%
%
% Thomas Rusterholz, 13 Sep 2019
%   Update 17 Aug 2020: added options unique and sort
%   Update 26 Mar 2021: also exports index, but only if there is only oe
%                       list (find other solulution for this)
%-------------------------------------------------------------------------

function [list,indList] = selectionList(list,files,opt)

%To delete GUI if failed
% delete(findall(0,'name','selectionList.m'))
% delete(findall(0,'type','figure'))

%CHECK INPUT
%list
if nargin<1 || isempty(list)
    list = {}; %return cell
    return
end
if ~iscell(list)
    error('1st input argument must be of class cell')
end
list = list(:); %column vector
%file
if nargin<2 || isempty(files) || ~ischar(files)
    noLIS = 1; %number of lists/axis
else
    noLIS = 2;
    ind = ... must be paths or files
        cellfun(@(x)exist(x,'dir'),list) | ...
        cellfun(@(x)exist(x,'file')==2,list);
    list(~ind)=[];
    if isempty(list)
        return
    end
end
%opt
if nargin<3 || ~isstruct(opt)
    opt = struct;
end
if ~isfield(opt,'unique')
    opt.unique = true;
end
if ~isfield(opt,'sort')
    opt.sort = true;
end
noTIT = 1;
if ~isfield(opt,'title')
    opt.title = ''; %default settings later
elseif iscell(opt.title)
    noTIT = numel(opt.title);
end

%PREPARE LIST
%sort
if opt.unique
    list = unique(list);
end
if opt.sort
    try %if natsort found
        tmp = cellfun(@(x)regexprep(x,{'\.','[\\/]'},{char(0),char(1)}),...
            list,'UniformOutput',false); %magic trick
        [~,ind] = natsort(tmp);
        list = list(ind);
    catch
        list = sort(list);
    end
end
%split index
if noLIS==2
    tmp = list;
    ind = cellfun(@(x)exist(x,'file')==2,tmp);
    tmp(ind) = cellfun(@fileparts,tmp(ind),'UniformOutput',false);
    tmp = fullfile(tmp,files);
    indLIS = false(size(tmp));
    for k = 1:numel(indLIS)
        indLIS(k) = isempty(dir(tmp{k}));
    end
end

%GUI POSITIONS in [pixel]
butW = 80; butH = 40; %buttons width & height
titH = 30;            %title height (auto width)  
lisW = 6*butW; lisH = 500;  %lists width (> 2*butW) & height
dd = 10;              %space for buttons, title, ...
dx = [20,20,20];      %lists x-spaces: left, between, right
dy = [butH+2*dd,0,noTIT*titH+2*dd]; %lists y-spaces: bottom, between, top
%calculate positions
pos.figure  = [200,200,noLIS*lisW+dx*[1;noLIS-1;1],lisH+sum(dy)]; %centered
pos.lists   = NaN(noLIS,4); %lists
pos.buttons = NaN(noLIS,4,2); %buttons
for k = 1:noLIS
    x=dx*[1;k-1;0]+(k-1)*lisW;
    pos.lists(k,:)     = [x ,dy(1) ,lisW ,lisH];
    pos.buttons(k,:,1) = [x         ,dd ,butW ,butH];
    pos.buttons(k,:,2) = [x+dd+butW ,dd ,butW ,butH];
end
pos.butOK = [sum(pos.lists(end,[1,3]))-butW,dd,butW,butH];
%title position, may have multipe ones
pos.title = cell(noTIT,1);
pos0 = [pos.lists(1,1),dy(1)+lisH+dd+noTIT*titH,pos.figure(3)-dy(3),titH];
for k = 1:noTIT
    pos0(2) = pos0(2)-titH;
    pos.title{k} = pos0;
end

%PATHS LISTS
if noLIS==2
    %added addPath for to deslect all.
    gui.list = {list(indLIS),list(~indLIS)};
    titleLists={...
        sprintf('''%s'' not exist in path, N = %i',files,sum(indLIS));...
        sprintf('''%s'' exist in path, N = %i',files,sum(~indLIS));...
        };
    if exist(list{1},'file')==2
        titleMain = {'Select Files'};
    else   
        titleMain = {'Select Paths'};
    end
else
    %added addPath for to deslect all
    gui.list={list};
    if exist(list{1},'file')==2
        titleLists={sprintf('Files, N = %i',numel(list))};
        titleMain = {'Select Files'};
    elseif exist(list{1},'dir')==7
        titleLists={sprintf('Paths, N = %i',numel(list))};
        titleMain = {'Select Paths'};
    else
        titleLists={sprintf('Elements, N = %i',numel(list))};
        titleMain = {'Select Elements'};
    end
    if ~isempty(opt.title)
        if ischar(opt.title)
            titleMain = {opt.title};
        else
            titleMain = opt.title;
        end
    end
end

%FIGURE
gui.h.figure=figure('position',pos.figure,...
    'name',sprintf('%s.m',mfilename),'CloseRequestFcn',@fun_close);
movegui(gui.h.figure,'center');

%PANNELS / BUTTONS / TITEL
gui.h.panels  = NaN(noLIS,1);
gui.h.lists   = NaN(noLIS,1);
gui.h.buttons = NaN(noLIS,2);
for k=1:noLIS
    %panel
    gui.h.panels(k)=uipanel('parent',gui.h.figure,...
        'unit','pixel','position',pos.lists(k,:),...
        'fontsize',12,'fontweight','bold','Title',titleLists{k});
    %list
    gui.h.lists(k)=uicontrol('parent',gui.h.panels(k),...
        'style','listbox',...
        'unit','normalized','position',[0,0,1,1],...
        'Max',2,'Min',0,... multiple selection
        'fontsize',10,'string',gui.list{k},'value',[]);
    %buttons
    gui.h.buttons(k,1) = uicontrol('parent',gui.h.figure,...
        'style','pushbutton','backgroundcolor',[1,1,0],...
        'unit','pixel','position',pos.buttons(k,:,1),...
        'fontsize',10,'string','Select All',...
        'callback',@fun_buttons);
    gui.h.buttons(k,2) = uicontrol('parent',gui.h.figure,...
        'style','pushbutton','backgroundcolor',[1,.5,.5],...
        'unit','pixel','position',pos.buttons(k,:,2),...
        'fontsize',10,'string','Deselect All',...
        'callback',@fun_buttons);
end
%ok button
gui.h.butOK=uicontrol('parent',gui.h.figure,...
    'style','pushbutton','backgroundcolor',[.2,1,.2],...
    ...'tooltip','Export 0 + 0 = 0',...
    'unit','pixel','position',pos.butOK,...
    'fontsize',12,'fontweight','bold','string','OK',...
    'callback',@fun_close);
%title
gui.h.title = cell(noTIT,1);
if noTIT==1
    fontsize = 18;
else
    fontsize = 15;
end
for k = 1:noTIT
    gui.h.title=uicontrol('parent',gui.h.figure,...
        'style','text',...
        'unit','pixel','position',pos.title{k},...
        'fontsize',fontsize,'fontweight','bold','string',titleMain{k});
end
%correction for single file-list
%if noLIS==1

%SETTINGS
%normalize unit
tmp=findall(gui.h.figure,'unit','pixel');
tmp(tmp==gui.h.figure)=[]; %would work as well
set(tmp,'unit','normalized')
%apply data
guidata(gui.h.figure,gui)
%figure always on top (if function exist)
try
    fig_alwaysOnTop(gui.h.figure,'fix')
catch
end
%wait for user interaction
waitfor(gui.h.figure);
end     

%FUCTIONS
%-------------------------
function fun_buttons(src,~)
gui=guidata(src);
[row,col] = find(ismember(gui.h.buttons,src));
switch col
    case 1
        set(gui.h.lists(row),'value',1:numel(gui.list{row}))
    case 2
        set(gui.h.lists(row),'value',[])
end
% %set export tooltip
% for k=1:numel(gui.h.lists)
%     N(k)=numel(get(gui.h.lists(k),'value'));
% end
% switch numel(N)
%     case 1
%         str=sprintf('Selected %i',N);
%     case 2
%         str=sprintf('Selected %i + %i = %i',N(1),N(2),sum(N));
% end
% set(gui.h.butOK,'tooltip',str)
end

function fun_close(src,~)
gui = guidata(src);
n   = numel(gui.list);
list={};
for k=1:n %may have two lists
    tmp  = gui.list{k};
    tmp  = tmp(get(gui.h.lists(k),'value'));
    list = [list;tmp(:)];
end
if n==1
    indList = get(gui.h.lists(k),'value')';
else
    indList = NaN; %find solution for this
end
assignin('caller','list',list)
assignin('caller','indList',indList)
delete(gui.h.figure); %PS: close does not work because it's THIS function
end