function point = gd_getpoint(ax,promptxt)
%
%-------function help------------------------------------------------------
% NAME
%   gd_getpoint.m
% PURPOSE
%   interactively select a point on a plot and return the point
%   coordinates.
% USAGE
%   point = gd_getpoint(ax,promptxt);
% INPUTS
%   ax - figure axes to use to interactively select point
%   promptxt - prompt to be used for point being selected
% OUTPUTS
%   point - struct with x and y fields defining selected point
% NOTES
%   point needs to be defined using gd_setpoint so that the callback
%   provides additional information about the type of point selection
%   (left or right mouse click) in the graphical point UserData property.
% SEE ALSO
%   called in gd_digitisepoints and gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
% 
    title(ax,promptxt)
    point = [];
    but = 999;
    while but~=0
        try
            but = waitforbuttonpress;
        catch
            return;
        end
        h_pnts = findobj(ax,'Tag','mypoints');
        if isempty(h_pnts)
            warndlg('No points to defined')
            return
        end
        %callback for setpoint sets UserData to event when point selected
        if all([h_pnts.UserData]==0)      %no point set
            but = 999;                    %stay in loop
        elseif any([h_pnts.UserData]==3)  %user right clicked a point
            return;                       %exit with no point data
        end        
    end
    idx = [h_pnts.UserData]>0;
    point.x = h_pnts(idx).XData;
    point.y = h_pnts(idx).YData;
end