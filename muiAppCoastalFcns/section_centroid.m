function [xc,zc,m0] = section_centroid(xi,zi)
%
%-------function help------------------------------------------------------
% NAME
%   section_centroid.m
% PURPOSE
%   Find centroid and area of a cross-section defined by points (xi,zi)
% USAGE
%   [xc,zc,m0] = section_centroid(xi,zi)
% INPUTS
%   xi - array of points of section in the x-dimension
%   zi - array of points of section in the x-dimension   
% OUTPUT
%   xc - centroid in the x-dimension 
%   zc - centroid in the z-dimension 
%   m0 - zero moment - area of section   
% NOTES
%   Uses numerical integration for centroid in both dimensions and
%   takes the average (this improves the estimate for poorly
%   sampled data sets).   
% SEE ALSO
%   Used in CT_BeachAnalysis.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%

    m0 = trapz(xi,zi);             %zero moment - area of section    
    x1 = trapz(xi,xi.*zi)/m0;
    x2 = -trapz(zi,xi.^2/2)/m0;
    xc = (x1+x2)/2;                %average of first moment in x
    z1 = trapz(xi,zi.^2/2)/m0;
    z2 = -trapz(zi,xi.*zi)/m0;
    zc = (z1+z2)/2;                %average of first moment in z
end