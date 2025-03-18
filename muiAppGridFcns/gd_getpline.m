function [pline,H] = gd_getpline(ax,promptxt,tagname,ispoints)
%
%-------function help------------------------------------------------------
% NAME
%   gd_getpline.m
% PURPOSE
%   interactively select a line on a plot and return the line point
%   coordinates.
% USAGE
%   [pline,H] = gd_getpline(ax,promptxt,tagname,ispoints);
% INPUTS
%   ax - figure axes to use to interactively select point
%   promptxt - prompt to be used for point being defined
%   tagname - character vector of text to be used as tag for plotted points
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
    if nargin<4
        ispoints = true;
    end
    title(ax,promptxt)
    pline = []; H = [];
    but = 999;
    while but~=0
        try
            but = waitforbuttonpress;
        catch
            return;
        end

        h_lines = findobj(ax,'Tag',tagname);
        if isempty(h_lines)
            warndlg('No lines defined')
            return
        end
        %callback for setpoint sets UserData to event when line selected
        if all([h_lines.UserData]==0)      %no line set
            but = 999;                     %stay in loop
        elseif any([h_lines.UserData]==3)  %user right clicked a line
            return;                        %exit with no line data
        end        
    end
    
    idx = [h_lines.UserData]>0;
    if sum(idx)>1, return; end             %multiple lines active

    if isempty(idx)
        return
    elseif isnan(h_lines(idx).XData(end))
        pline.x = h_lines(idx).XData';    %xy of lines are column vectors
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