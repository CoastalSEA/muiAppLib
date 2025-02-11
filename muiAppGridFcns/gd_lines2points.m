function [points,outype] = gd_lines2points(lines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_lines2points.m
% PURPOSE
%   convert input x,y vectors in various formats to an array of structs 
%   for each point with x,y (and z) fields
% USAGE
%   [points,outype] = gd_lines2points(lines);
% INPUTS
%   lines - x, y (and z) vectors in any of the folloing formats:
%             outype=0: array of structs with x, y and z fields defining selected points,
%             outype=1: Nx2 or Nx3 array.
%             outype=2: struct with x, y (and z) column vectors
%             outype=3: table with x, y (and z) column vectors
% OUTPUTS
%   points - a row struct array of x, y and optionally z
%   outype - format of input vecpnts
% NOTES
%   digitising functions work with an array of structs each one defining a
%   point. The output is converted to one of the formats defined for
%   lines. This function reverses any of these formats to the points format
% SEE ALSO
%   used to reformat input in gd_editlines
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
% 
    if isstruct(lines) && length(lines)>1 && isfield(lines(1),'x')
        points = lines;             %return unchanged
        outype = 0;
        return;
    end

    points = struct('x',[],'y',[]);
    if isstruct(lines) || istable(lines)
        points = assignPoints(points,lines.x,'x');
        points = assignPoints(points,lines.y,'y');
        if isstruct(lines)
            outype = 2;
            if isfield(lines,'z')
                points = assignPoints(points,lines.z,'z');
            end
        else
            outype = 3;
            if ismember('z', lines.Properties.VariableNames)
                points = assignPoints(points,lines.z,'z');
            end
        end
    elseif ismatrix(lines)  %NB this assumes two column vectors for x and y
        %corrects large vectors but may fail for short vectors <=3 points
        if size(lines,2)>3, lines = lines'; end  
        points = assignPoints(points,lines(:,1),'x');
        points = assignPoints(points,lines(:,2),'y');
        if size(lines,2)==3
            points = assignPoints(points,lines(:,3),'z');
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