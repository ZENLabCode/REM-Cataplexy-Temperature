% fig_size
%------------------------------------------------------------------------
% Set figure position AND print/save size.
%
% Sometimes MATLAB changes print/save size of figures, especially if
% figure screen size exceeds default print size. With this function,
% saved figure will have exactly the same size as displayed on screen.
%
% SYNTAX
%   trf_figSize(h)
%     - h is figure handle (multiple possible)
%       sets print/save size to plotted figure size
%   trf_figSize(h,position,units)
%     - position (optional), figure position in the form
%       rect = [left, bottom, width, height]
%       default is actual position
%     - unit (optional), figure unit e.g. 'pixel', 'inches', ...
%       default is actual unit
%
%
% Thomas Rusterholz, July 28 2014
%------------------------------------------------------------------------

function fig_size(hf,position,unit)
%Important: figure needs to be drawn before changing anything. Drawnow
%somehow does not work for unknown reasons. So I included a pause.
drawnow; pause(0.2)

%FIGURE LOOP
for ind=1:numel(hf)
    h=hf(ind);
    
    %position and unit
    if exist('position','var') && ~isempty(position)
        pos=position;
    else
        pos=get(h,'position');
    end
    if exist('unit','var') && ~isempty(unit)
        uni=unit;
    else
        uni=get(h,'unit');
    end
    
    %set figure screen position
    set(h,'unit',uni,'position',pos)
    while pos~=get(h,'position') %sometimes too slow
        drawnow
        set(h,'position',pos);
    end
    
    %set paper size/position
    try
        set(h,'paperUnits',uni,'paperSize',pos(3:4),...
            'paperPosition',[0 0 pos(3:4)])
    catch %e.g. unit 'pixel' does not work ad paperUnit
        set(h,'unit',get(h,'paperUnits'))
        pos=get(h,'position');
        set(h,'paperSize',pos(3:4),'paperPosition',[0 0 pos(3:4)])
        set(h,'units',uni) %reset screen unit
    end
end