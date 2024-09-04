function gdims = gd_dimensions(grid)
%-------function help------------------------------------------------------
% NAME
%   gd_dimensions.m
% PURPOSE
%   get the grid dimensions for a grid struct (as used in classes that 
%   inherit GDinterface)
% USAGE
%   gdims = gd_dimensions(grid)
% INPUTS 
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
% OUTPUTS
%   gdims - table of delx,dely,xint,yint,xmin,xmax,ymin,ymax,zmin,zmax
%       delx, dely - grid interval in the x and y dimensions
%       xsgn, ysgn - direction of y, y axes (-ve is descending)
%       xint, yint - number of grid intervals in the x and y dimensions
%       xmin,xmax,ymin,ymax,zmin,zmax - min and max x, y and z dimensions
% NOTES
%   function replaces static method GDinterface.gridDims which returns a
%   table of delx,dely,xint,yint,xmin,xmax,ymin,ymax
% SEE ALSO
%   function is also similar to setGridDimensions method in GD_GridProps
%
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%
    delx = abs(grid.x(2)-grid.x(1));          
    dely = abs(grid.y(2)-grid.y(1));
    xsgn = sign(grid.x(2)-grid.x(1));  %direction of axes (-ve descending)
    ysgn = sign(grid.y(2)-grid.y(1)); 
    xint = length(grid.x);
    yint = length(grid.y);
    xmin = min(grid.x);
    xmax = max(grid.x);
    ymin = min(grid.y);
    ymax = max(grid.y);
    zmin = min(grid.z,[],'all');
    zmax = max(grid.z,[],'all');
    %output table
    gdims = table(delx,dely,xsgn,ysgn,xint,yint,xmin,xmax,ymin,ymax,zmin,zmax);
end