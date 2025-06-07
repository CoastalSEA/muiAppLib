function [clinedir,ncplines,cumlen] = gd_curvelineprops(cplines,idL)
%
%-------function help------------------------------------------------------
% NAME
%   gd_curvelineprops.m
% PURPOSE
%   for each point from idL to the end use the centre-line coordinates and 
%   direction to find the lengths and directions along the centre-line
% USAGE
%   [clinedir,ncplines,cumlen] = gd_curvelineprops(cplines,idL)
% INPUTS
%   cplines - cell array of lines with x,y struct defining line points
%   idL - index of a start point on one of the lines in cplines (1 for
%   first point in line)
% OUTPUTS
%   clinedir - mean direction of line at each point in x,y space
%   ncplines - updated cplines starting from point idL
%   cumlen - cumulative length from the defined start point
% SEE ALSO
%   used in PL_SectionLines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%    
    if nargin<3,  end
    nlines = length(cplines);
    %cumlen is local reach specific cumulative lengths from start point
    cumlen = cell(1,nlines); clinedir = cumlen; 
    nrec = length(cplines{1,1});               %number of points in first line
    %nrec = 0;
    j = 1;                                     %count of lines included
    for i=1:nlines
        lp = cplines{1,i};
        nl = length(lp);
        if idL>=nrec                           %start point not in line
            nrec = nrec+length(cplines{1,i+1});
            continue;
        elseif ~exist('dx','var')              %start point in line
            idl = idL-(nrec-nl);               %index of start point in line
            dx = diff([lp(idl:end-1).x]);      %omit trailing NaN
            dy = diff([lp(idl:end-1).y]);   
            ncplines{1,j} = lp(idl:end);       %#ok<AGROW> %crop line to start point
        else                                   %subsequent lines
            dx = diff([lp(1:end-1).x]);        %omit trailing NaN
            dy = diff([lp(1:end-1).y]);  
            ncplines{1,j} = lp;                %#ok<AGROW> %add subsequent lines
        end
        
        if ~isempty(dx)                        %trap single point at end of line
            slen = hypot(dx,dy);               %length between points
            cumlength = cumsum(slen);          %cumulative length
            cumlen{j} = [0,cumlength,NaN];     %pad to make same length as lines

            %now pad dx,dy to make angles same length as lines
            dx = [dx(1),dx,dx(end)];           %#ok<AGROW> 
            dy = [dy(1),dy,dy(end)];           %#ok<AGROW>             
            theta = atan2(dy,dx);              %direction between points            
            dirs(1) = theta(1);
            for k=2:length(theta)-1             %mean direction at point
                dirs(k) = (theta(k-1)+theta(k+1))/2; %#ok<AGROW> 
            end  
            clinedir{j} = [dirs,NaN];           %pad to make same length as lines
            dirs = [];
            j = j+1; 
        end
    end  
end    