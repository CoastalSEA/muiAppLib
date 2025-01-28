function points = gd_setpoints(ax,promptxt,isxyz)
%
%-------function help------------------------------------------------------
% NAME
%   gd_setpoints.m
% PURPOSE
%   interactively create a set of points on a plot and return the point
%   coordinates. Includes an option to enter an additional value at the
%   selected points (e.g. elevation).
% USAGE
%   points = gd_setpoints(ax,promptxt,isxyz);
% INPUTS
%   ax - figure axes to use to interactivvely select point
%   promptxt - prompt to be used for point being defined
%   isxyz - logical flag true to input z values - optional, default is false
% OUTPUTS
%   points - struct with x, y fields defining added points and z if included 
% NOTES
%   uses gd_setpoint to get each point
% SEE ALSO
%   called in gd_digitisepoints and gd_selectpoints???????????
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
    %points = struct('x',[],'y',[]);
    count = 1;
    ok = 0;
    while ok<1          
        point = gd_setpoint(ax,promptxt,isxyz);
        if isempty(point)  %user right clicks or presses return
            if ~exist('points','var'), points = []; end
            ok = 1;       
        else
            points(count) = point; %#ok<AGROW> 
            count = count+1;
        end        
    end
end
