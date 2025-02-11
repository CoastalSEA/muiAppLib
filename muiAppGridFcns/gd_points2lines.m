function lines = gd_points2lines(points,outype)
%
%-------function help------------------------------------------------------
% NAME
%   gd_points2lines.m
% PURPOSE
%   convert an array of structs with x,y (and z) fields to a [Nx2] or [Nx3] 
%   array of points, or a single stuct with vectors in the x, y (and z) fields
% USAGE
%   lines = gd_points2lines(points,outype);
% INPUTS
%   points - a row struct array of x, y and optionally z 
%   outype - format of output - see Outputs for details
% OUTPUTS
%   lines - outype=0: array of structs with x, y and z fields defining selected points,
%           outype=1: Nx2 or Nx3 array.
%           outype=2: struct with x, y (and z) column vector fields
%           outype=3: table with x, y (and z) column vector fields
% NOTES
%   digitising functions work with an array of structs each one defining a
%   point. The output is converted to one of the formats defined for
%   vecpnts. This function converts the points format to any of the lines
%   formats.
% SEE ALSO
%   used to reformat output from gd_digitisepoints and gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%    
    if nargin<2 || outype==0
        %output type undefined so cannot determine what to do
        lines = points; return;                
    elseif ~isstruct(points)
        %not points array struct as input - return unchanged
        lines = points; return;         
    elseif (isstruct(points) && length(points)==1)
        %if struct of single point return unchanged, otherwise
        if outype==1      %convert to array
            lines = cell2mat(struct2cell(points)');    %struct to [Nx2] array
        elseif outype==3  %convert to table
            lines = struct2table(points);
        end
        return;
    end

    x = [points(:).x]; if isrow(x), x = x'; end
    y = [points(:).y]; if isrow(y), y = y'; end

    if outype==1
        lines = [x,y];
    elseif outype==2
        lines.x = x;
        lines.y = y;
    elseif outype==3
        lines = table(x,y,'VariableNames',{'x','y'});
    else
        warndlg('Unrecognised format type (outype should be 0-3)')
        return
    end

    if isfield(points,'z')
        z = [points(:).z];  if isrow(z), z = z'; end
        if outype==1
            lines = [lines,z];
        elseif outype==2
            lines.z = z;
        elseif outype==3
            lines = addvars(lines,z,'NewVariableNames','z');
        end
    end
end