function [X,Y,Z] = gd_grid_line(varargin)
%
%-------function help------------------------------------------------------
% NAME
%   gd_grid_line.m
% PURPOSE
%   create a grid from scattered data points input as xyz tuples
% USAGE
%   [X,Y,Z] = gd_grid_line(xin,yin,zin,missing)
%   or
%   [X,Y,Z] = gd_grid_line(xyz,missing)
% INPUTS
%   xin,yin,zin - cartesian vectors
%   OR
%   xyz - cartesian matrix
%   missing - value used to infill missing values (e.g. 0 or NaN)
% OUTPUT
%   X,Y,Z - grid matrices
% NOTES
%   Derived from Matlab Forum code by Jack Kohoutek 2009 jack.kohoutek@gmail.com
%   http://www.mathworks.com/matlabcentral/fileexchange/4551-xyzplotter
%   Repeat xin and yin values are fine, they are not repeatedly added to the output matrices. 
%   The cell size of the resulting grid is therefore the [number of unique
%   xin values x the number of unique yin values].
%   The function xyzplotter uses @sum when callin accumarray. This sums
%   duplicate values at the same xy coordinates. Changed to @mean so that
%   duplicates are removed (assumes they have the same z value)
% SEE ALSO
%   used in GDinterface.addValleyBase
%
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%
% 
    if length(varargin)==4
        xin = varargin{1}; yin = varargin{2}; zin = varargin{3}; 
        missing = varargin{4};
        if length(xin)~=length(yin)||length(xin)~=length(zin)||length(yin)~=length(zin)
            error('Length of arguments must be equal');
        end
        if isrow(xin), xin = xin'; end  %force column vectors
        if isrow(yin), yin = yin'; end  
        if isrow(zin), zin = zin'; end  
    elseif length(varargin)==2
        xyz = varargin{1};
        if size(xyz,1)==3, xyz = xyz'; end  %force column vectors
        xin = xyz(:,1); yin = xyz(:,2); zin = xyz(:,3); 
        missing = varargin{2};
    else
        error('Incorrect number of inputs')
    end

    %find unique x and y and bin values
    xu = unique(xin);
    yu = unique(yin);

    binx = discretize(xin,xu);
    biny = discretize(yin,yu);

    [X,Y] = meshgrid(xu,yu);
    z = accumarray([binx,biny],zin,[size(X,2) size(Y,1)],@mean,missing);
    Z = z';
end