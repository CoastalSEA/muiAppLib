function points = gd_selectpoints(grid,npts,promptxt,isdel)
%npts,
%-------function help------------------------------------------------------
% NAME
%   gd_selectpoints.m
% PURPOSE
%   Accept figure to interactively select one or more points on a grid
% USAGE
%   points = gd_selectpoints(grid,npts,promptxt,isdel)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   npts - number of points to be selected
%   promptxt - cell array of prompts to be used for each point being defined
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - struct with x and y fields defining selected points
% NOTES
%   to trap use of 'Quit' use >> if any(isnan([points(:).x])), return; end
%   after call to gd_selectpoints in calling function 
% SEE ALSO
%   called in GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    if nargin<3, isdel = false; end
    
    figtitle = sprintf('Select points');
    figtxt = 'Use left mouse click to mark point, right click to accept';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Accept','Edit','Quit'};
    position = [0.2,0.4,0.4,0.4];
    [h_plt,h_but] = acceptfigure(figtitle,figtxt,tag,butnames,position);
    ax = gd_plotgrid(h_plt,grid);
    
    %get user to define the required points
    points = setPoints(ax,npts,promptxt);
    
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            return;             
        elseif strcmp(h_but.Tag,'Edit') 
            points = setPoints(ax,npts,promptxt);
            h_but.Tag = 'reset';
        elseif strcmp(h_but.Tag,'Quit') 
            points.x = NaN;   points.y = NaN;
            ok = 1;
        else   %user accepted
            ok = 1; 
        end        
    end

    %delete figure if isdel has been set by call.
    if isdel
        delete(h_plt.Parent)
    end
end
%%
function points = setPoints(ax,npts,promptxt)
    %use mouse to select points
    h_pnts = findobj(ax,'Tag','mypoints');
    delete(h_pnts)  %remove any existing points
    hold on
    for i=1:npts
        if length(promptxt)~=npts  && ~isempty(promptxt)
            prompt = sprintf('%s %d',promptxt{1},i);
        else
            prompt = promptxt{i};
        end
        points(i) = gd_setpoint(ax,prompt,false); %#ok<AGROW>
    end
end