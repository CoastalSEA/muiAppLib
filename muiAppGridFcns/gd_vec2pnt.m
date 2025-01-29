function [points,outype] = gd_vec2pnt(vecpnts)
%
%-------function help------------------------------------------------------
% NAME
%   gd_vec2pnt.m
% PURPOSE
%   convert input x,y vectors in various formats to an array of structs 
%   for each point with x,y (and z) fields
% USAGE
%   vecpnts = gd_vec2pnt(points);
% INPUTS
%   vecpnts - outype=0: array of structs with x, y and z fields defining selected points,
%             outype=1: Nx2 or Nx3 array.
%             outype=2: struct with x, y (and z) vector fields
%             outype=3: table with x, y (and z) vector fields
% OUTPUTS
%   points - struct array of x, y and optionally z 
%   outype - format of input vecpnts
% NOTES
%   digitising functions work with an array of structs each one defining a
%   point. The output is converted to one of the formats defined for
%   vecpnts. This function reverses any of these formats to the points format
% SEE ALSO
%   used to reformat input in gd_editlines
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
% 
    if isstruct(vecpnts) && length(vecpnts)>1 && isfield(vecpnts(1),'x')
        points = vecpnts;             %return unchanged
        outype = 0;
        return;
    end

    points = struct('x',[],'y',[]);
    if isstruct(vecpnts)
        [points.x] = vecpnts.x;
        [points.y]= vecpnts.y;
        if isfield(vecpnts,'z')
            [points(:).z] = vecpnts.z;
        end
        outype = 2;
    elseif istable(vecpnts)
        points = assignPoints(points,vecpnts.x,'x');
        points = assignPoints(points,vecpnts.y,'y');
        if ismember('z', vecpnts.Properties.VariableNames)
            points = assignPoints(points,vecpnts.z,'z');
        end
        outype = 3;
    elseif ismatrix(vecpnts)  %NB this assumes two column vectors for x and y
        [points.x] = vecpnts(:,1);
        [points.y] = vecpnts(:,2);
        if size(vecpnts,2)==3
            [points(:).z] = vecpnts(:,3);
        end
        outype = 1;
    else
        warndlg('Unrecognised input format')
        points = []; outype = [];
        return;
    end
end
%%
function points = assignPoints(points,vec,var)
    %unpack struct or table vector and assign to points
    nrec = length(vec);
    inpv = num2cell(vec);
    [points(1,1:nrec).(var)] = inpv{:};
end