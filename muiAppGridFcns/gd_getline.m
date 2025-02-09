function line = gd_getline(ax,promptxt,ispoints)
%
%-------function help------------------------------------------------------
% NAME
%   gd_getline.m
% PURPOSE
%   interactively select a line on a plot and return the line point
%   coordinates.
% USAGE
%   line = getline(ax,promptxt,ispoints);
% INPUTS
%   ax - figure axes to use to interactivvely select point
%   promptxt - prompt to be used for point being defined
%   ispoints - true to return as an array of spoint structs, otherwise
%              returns an xy struct of points (optional - default is false)
% OUTPUTS
%   line - struct with x and y fields defining selected line or a struct
%          array of points depending on value of ispoints
% NOTES
%    used with gd_setlines in functions such as gd_sectionlines
% SEE ALSO
%   called in gd_sectionlines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
% 
    if nargin<3
        ispoints = false;
    end
    title(ax,promptxt)
    line = [];
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
            but = 999;                    %stay in loop
        elseif any([h_lines.UserData]==3)  %user right clicked a point
            return;                       %exit with no point data
        end        
    end
    idx = [h_lines.UserData]>0;
    line.x = h_lines(idx).XData;
    line.y = h_lines(idx).YData;
    if ispoints
        line = gd_vec2pnt(line);
    end
end