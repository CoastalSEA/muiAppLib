function cplines = gd_plines2cplines(plines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_plines2cplines.m
% PURPOSE
%   convert an array of plines to a cell array of plines
% USAGE
%   cplines = gd_plines2cplines(plines);
% INPUTS
%   plines - struct array with x and y fields defining one or more plines
% OUTPUTS
%   cplines - a cell array of plines
% NOTES
%   called in gd_sectionlines
% SEE ALSO
%   inverse function gd_cplines2plines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
% 
    idN = [0,find(isnan([plines(:).x]))];  
    if isscalar(idN), cplines = {plines}; return; end  %no NaNs in plines
    for i=1:length(idN)-1
        idL =idN(i)+1:idN(i+1);    
        %cell so that number of points can vary
        cplines{1,i} = plines(idL); %#ok<AGROW> 
    end
end