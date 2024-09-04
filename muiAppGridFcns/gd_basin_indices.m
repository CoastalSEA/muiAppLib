function [idx,ixM,ixH] = gd_basin_indices(grid,xM,Lt)
%
%-------function help------------------------------------------------------
% NAME
%   gd_basin_indices.m
% PURPOSE
%   get the indices of the grid that fall within the basin or channel,
%   when the mouth is offset from the grid origin.
% USAGE
%   [idx,ixM,ixH] = gd_basin_indices(grid,xM,Lt)
% INPUTS
%   grid - struct of x,y,z values that define grid as well as xM and Lt,
%          distances to mouth and from mouth to tidal limit
%   xM - distance to mouth from origin of grid (min(x)) - optional
%   Lt - distance from mouth to tidal limit (m) - optional
%        if xM and/or Lt are not specified then the values in the grid struct
%        are used. If Lt=0, or if Lt not used and grid.Lt=0, defaults to
%        the end of the grid.
% OUTPUT
%   idx - indices of x-axis that comprise the channel
%   ixM - index of mouth position on x ax is
%   ixH - index of head or tidal limit position on x axis
% NOTES
%   checks direction of x-axis and orientation of channel relative to axis
%   to determine the indices within the channel taking account of any
%   offset, xM. Similar to grid orientation cases returned by gd_ax_dir
% SEE ALSO
%   used in gd_basin_hypsometry, gd_section_properties and gd_gross_properties
%   
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%
    if nargin<2
       xM = grid.xM;
       Lt = grid.Lt;
    elseif nargin<3
       Lt = grid.Lt; 
    elseif nargin==3 && isempty(xM)
       xM = grid.xM;
    end

    gdim = gd_dimensions(grid);
    xM = min(grid.x)+xM;            %xM defined as distance from min(x)
    if grid.ishead                  %head is nearest the x-origin
        xH = min(grid.x);
        if Lt>0 && Lt<xM            %check that Lt defined and within grid
            xH = xM-Lt;
        end
    else
        xH = max(grid.x);
        if Lt>0 && xM+Lt<xH         %check that Lt defined and within grid
            xH = xM+Lt;
        end
    end
    
    %now find index of mouth and head
    [~,ixM] = min(abs(grid.x-xM));  %find nearest x-index to xM
    [~,ixH] = min(abs(grid.x-xH));  %find nearest x-index to xH
     
    if grid.ishead*(gdim.xsgn>0)    %ascending x-axis
        idx = ixH:ixM;              %indices from head to mouth
    else                            %descending x-axis
        idx = ixM:ixH;              %indices from mouth to head
    end
end