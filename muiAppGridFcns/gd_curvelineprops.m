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
%   idL - x,y start point on one of the lines in cplines
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
    nrec = length(cplines{1,1});               %length of first line
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
            %pad to make same length as lines
            dx = [dx(1),dx,dx(end)];           %#ok<AGROW> 
            dy = [dy(1),dy,dy(end)];           %#ok<AGROW> 
    
            slen = hypot(dx,dy);               %length between points
            cumlength = cumsum(slen);          %cumulative length
            cumlen{j} = [0,cumlength(2:end-1),NaN];
            theta = atan2(dy,dx);              %direction between points
            clinedir{j}(1) = theta(1);
            for k=2:length(theta)              %mean direction at point
                clinedir{j}(k) = (theta(k-1)+theta(k))/2;
            end  
            j = j+1; 
        end
    end  
end    