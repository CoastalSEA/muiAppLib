function [points,Hpnts] = gd_setpoints(ax,promptxt,isxyz)
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
%   ax - figure axes to use to interactively select point
%   promptxt - prompt to be used for point being defined
%   isxyz - logical flag true to input z values - optional, default is false
% OUTPUTS
%   points - struct with x, y fields defining added points and z if included 
%   Hpnts - handle to graphical point objects
% NOTES
%   uses gd_setpoint to get each point
%    NB: returns points type and NOT a line, which has a NaN termination.
% SEE ALSO
%   called in gd_digitisepoints and gd_editlines
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
    count = 1;
    ok = 0;
    while ok<1          
        [point,H] = gd_setpoint(ax,promptxt,isxyz);
        if isempty(point)  %user right clicks or presses return
            if ~exist('points','var'), points = []; end
            ok = 1;       
        else
            points(count) = point; %#ok<AGROW> 
            Hpnts(count) = H;      %#ok<AGROW> 
            count = count+1;
        end        
    end
end
