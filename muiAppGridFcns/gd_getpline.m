function [pline,H] = gd_getpline(ax,promptxt,ispoints)
%
%-------function help------------------------------------------------------
% NAME
%   gd_getpline.m
% PURPOSE
%   interactively select a line on a plot and return the line point
%   coordinates.
% USAGE
%   pline = getline(ax,promptxt,ispoints);
% INPUTS
%   ax - figure axes to use to interactively select point
%   promptxt - prompt to be used for point being defined
%   ispoints - true to return as an array of point structs, otherwise
%              returns an xy struct of points (optional - default is true)
% OUTPUTS
%   pline -  or a struct array of points defining selected line or a struct
%            with x and y fields, depending on value of ispoints
%   H -  handle to selected line
% NOTES
%   pline needs to be defined using gd_setlines so that the callback
%   provides additional information about the type of point selection
%   (left or right mouse click) in the graphical point UserData property.
% SEE ALSO
%   called in gd_sectionlines and see gd_setlines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
% 
    if nargin<3
        ispoints = true;
    end
    title(ax,promptxt)
    pline = [];
    but = 999;
    while but~=0
        try
            but = waitforbuttonpress;
        catch
            return;
        end

        h_lines = findobj(ax,'Tag','mylines');
        if isempty(h_lines)
            warndlg('No lines to defined')
            return
        end
        %callback for setpoint sets UserData to event when point selected
        if all([h_lines.UserData]==0)      %no point set
            but = 999;                     %stay in loop
        elseif any([h_lines.UserData]==3)  %user right clicked a point
            return;                        %exit with no point data
        end        
    end
    
    idx = [h_lines.UserData]>0;
    if isnan(h_lines(idx).XData(end))
        pline.x = h_lines(idx).XData';     %xy of lines are column vectors
        pline.y = h_lines(idx).YData';        
    else
        pline.x = [h_lines(idx).XData';NaN]; %xy of lines are column vectors
        pline.y = [h_lines(idx).YData';NaN];
    end

    if ispoints
        pline = gd_lines2points(pline);
    end
    H = h_lines(idx);                      %return handle to selected line
end