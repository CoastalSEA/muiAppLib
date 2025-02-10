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
    if isstruct(line) && length(line)>1 && isfield(line(1),'x')
        points = line;             %return unchanged
        outype = 0;
        return;
    end

    points = struct('x',[],'y',[]);
    if isstruct(line) || istable(line)
        points = assignPoints(points,line.x,'x');
        points = assignPoints(points,line.y,'y');
        if isstruct(line)
            outype = 2;
            if isfield(line,'z')
                points = assignPoints(points,line.z,'z');
            end
        else
            outype = 3;
            if ismember('z', line.Properties.VariableNames)
                points = assignPoints(points,line.z,'z');
            end
        end
    elseif ismatrix(line)  %NB this assumes two column vectors for x and y
        %corrects large vectors but may fail for short vectors <=3 points
        if size(line,2)>3, line = line'; end  
        points = assignPoints(points,line(:,1),'x');
        points = assignPoints(points,line(:,2),'y');
        if size(line,2)==3
            points = assignPoints(points,line(:,3),'z');
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