function vecpnts = gd_pnt2vec(points,outype)
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
%   outype - format of output - see Outputs for details
% OUTPUTS
%   vecpnts - outype=0: array of structs with x, y and z fields defining selected points,
%             outype=1: Nx2 or Nx3 array.
%             outype=2: struct with x, y (and z) vector fields
%             outype=3: table with x, y (and z) vector fields
% NOTES
%   digitising functions work with an array of structs each one defining a
%   point. The output is converted to one of the formats defined for
%   vecpnts. This function converts the points format to any of the vecpnts 
%   formats.
% SEE ALSO
%   used to reformat output from gd_digitisepoints and gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
% 
    if nargin<2 || outype==0
        vecpnts = points; return;            %return unchanged        
    elseif ~isstruct(points)
        %not points array struct as input - return unchanged
        vecpnts = points; return;        
    elseif (isstruct(points) && length(points)==1)
        vecpnts = points; %return struct of vectors unchanged
        if outype==1      %convert to array
            vecpnts = cell2mat(struct2cell(points)');    %struct to [Nx2] array
        elseif outype==3  %convert to table
            vecpnts = struct2table(points);
        end
        return;
    end

    x = [points(:).x]; if isrow(x), x = x'; end
    y = [points(:).y]; if isrow(y), y = y'; end

    if outype==1
        vecpnts = [x,y];
    elseif outype==2
        vecpnts.x = x;
        vecpnts.y = y;
    elseif outype==3
        vecpnts = table(x,y,'VariableNames',{'x','y'});
    else
        warndlg('Unrecognised format type (outype should be 0-3)')
        return
    end

    if isfield(points,'z')
        z = [points(:).z];  if isrow(z), z = z'; end
        if outype==1
            vecpnts = [vecpnts,z];
        elseif outype==2
            vecpnts.z = z;
        elseif outype==3
            vecpnts = addvars(vecpnts,z,'NewVariableNames','z');
        end
    end
end