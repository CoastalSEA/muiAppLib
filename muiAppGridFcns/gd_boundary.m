function points = gd_boundary(grid,figtitle,promptxt,outype,npts,isdel)
%npts,
%-------function help------------------------------------------------------
% NAME
%   gd_boundary.m
% PURPOSE
%   Accept figure to interactively select start and end points on a grid
% USAGE
%   points = gd_boundary(grid,figtitle,promptxt,outype,npts,isdel);
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   figtitle- character string used for title of figure
%   promptxt - cell array of prompts to be used for each point being defined
%              a single char string cell is apended with a count for each point,
%              whereas a cell array should have a length of npts.
%   outype - format of output - see Outputs for details
%   npts - number of points to be selected
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - outype=0: array of structs with x, y and z fields defining selected points,
%            outype=1: Nx2 or Nx3 array.
%            outype=2: struct with x, y (and z) vector fields
%            points = [] if user closes figure, or no points defined
%            outype=1:
%  points - struct with x and y fields defining selected start and end points
% NOTES
% 
% SEE ALSO
%   called in
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%

end