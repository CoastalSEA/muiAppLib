function vecpnts = gd_pnt2vec(points,isarray)
%
%-------function help------------------------------------------------------
% NAME
%   gd_pnt2vec.m
% PURPOSE
%   convert an array of structs with x,y (and z) fields to a [Nx2] or [Nx3] 
%   array of points, or a single stuct with vectors in the x, y (and z) fields
% USAGE
%   vecpnts = gd_pnt2vec(points,isarray);
% INPUTS
%   points - struct array of x, y and optionally z 
%   isarray - true to return array, false to return a struct (optional -
%             default is false
% OUTPUTS
%   vecpnts - struct with x, y (and z) vector fields or Nx2 or Nx3 array.
% SEE ALSO
%   used to reformat output from gd_digitisepoints and gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
% 
    if nargin<2
        isarray = false;  %return a struct of x,y vectors
    end

    x = [points(:).x]; if isrow(x), x = x'; end
    y = [points(:).y]; if isrow(y), y = y'; end

    if isarray
        vecpnts = [x,y];
    else
        vecpnts.x = x;
        vecpnts.y = y;
    end

    if isfield(points,'z')
        z = [points(:).z];  if isrow(z), z = z'; end
        if isarray
            vecpnts = [vecpnts,z];
        else
            vecpnts.z = z;
        end
    end
end