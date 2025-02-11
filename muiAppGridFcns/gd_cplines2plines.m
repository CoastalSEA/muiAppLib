function plines = gd_cplines2plines(cplines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_cplines2plines.m
% PURPOSE
%   convert cell array of plines to an array of points that define plines
% USAGE
%   plines = gd_cplines2plines(cplines);
% INPUTS
%   cplines - a cell array of plines
% OUTPUTS
%   plines - struct array with x and y fields defining a set of points
% NOTES
%   called in gd_sectionlines
% SEE ALSO
%   inverse function gd_plines2cplines
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
% 
    nrec = length(cplines);
    plines = [];
    for i=1:nrec
        plines = [plines,cplines{1,i}]; %#ok<AGROW> 
    end
end