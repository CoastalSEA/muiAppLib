function [sortE,sortN,idist] = sortENdata2line(Easting,Northing)
%
%-------function help------------------------------------------------------
% NAME
%   sortENdata2line.m
% PURPOSE
%   sort the eastings and northings into an order that makes the best
%   continuous line
% USAGE
%   [sortE,sortN,idist] = sortENdata2line(Easting,Northing)
% INPUTS
%   Easting - vector of x-coordinates
%   Northings - vector of y-coordinates
% OUTPUT
%   sortE - vector of sorted x-coordinates
%   sortN - vector of sorted y-coordinates
%   idsit - indices of sorting
% NOTES
%   sorting is based on distance from the point with the minimum Easting 
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    [E1,iEmin] = min(Easting,[],'omitnan');
    N1 = Northing(iEmin);
    dist = sqrt((Easting-E1).^2+(Northing-N1).^2);
    [~,idist] = sort(dist);
    sortE = Easting(idist);
    sortN = Northing(idist);
end
        