function points = gd_digitisepoints(grid,promptxt,isxyz,isdel)
%
%-------function help------------------------------------------------------
% NAME
%   gd_digitisepoints.m
% PURPOSE
%   Accept figure to interactively digitise points on a grid and add
%   elevations if required
% USAGE
%   points = gd_digitisepoints(grid,promptxt,isxyz,isdel)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   promptxt - cell array of prompts to be used for each point being defined
%   isxyz - logical flag true to input z values - optional, default is false
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - struct with x and y fields defining selected points
% NOTES
%   to trap use of 'Quit' use >> if any(isnan([points(:).x])), return; end
%   after call to gd_digitisepoints in calling function 
% SEE ALSO
%   called in GDinterface.getGridLine, similar to gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    if nargin<4
        isdel = false; 
    elseif nargin<3
        isxyz = false; 
        isdel = false;
    end
 
    figtitle = sprintf('Digitise points');
    figtxt = 'Use right click to accept';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Accept','Add','Edit','Quit'};
    position = [0.2,0.4,0.4,0.4];
    [h_plt,h_but] = acceptfigure(figtitle,figtxt,tag,butnames,position);
    ax = gd_plotgrid(h_plt,grid);
    
    %h_plt.Parent.ToolBar = 'none';
    %get user to define the required points
    %points = setPoints(ax,npts,promptxt);
    count = 1;
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            points(count).x = NaN;   points(count).y = NaN; %#ok<AGROW>
            return;     
        elseif strcmp(h_but.Tag,'Add') 
            points(count) = gd_setpoint(ax,promptxt,isxyz); %#ok<AGROW>
            count = count+1;
            h_but.Tag = 'reset';
        elseif strcmp(h_but.Tag,'Edit') 
            h_pnts = findobj(ax,'Tag','mypoints');
            idx = [h_pnts(:).XData]==points(count-1).x & [h_pnts(:).YData]==points(count-1).y; 
            delete(h_pnts(idx))  %remove any existing points
            points(count-1) = []; %#ok<AGROW>
            points(count-1) = gd_setpoint(ax,promptxt,isxyz); %#ok<AGROW>
            h_but.Tag = 'reset';
        elseif strcmp(h_but.Tag,'Quit') 
            points(count).x = NaN;   points(count).y = NaN; %#ok<AGROW>
            ok = 1;
            delete(h_plt.Parent)
        else   %user accepted
            h_pnts = findobj(ax,'Tag','mypoints');
            delete(h_pnts)
            hold on
            plot(ax,[points(:).x],[points(:).y],'-ok','MarkerEdgeColor','r')
            hold off
            ok = 1; 
        end        
    end

    %delete figure if isdel has been set by call.
    if isdel
        delete(h_plt.Parent)
    end
end