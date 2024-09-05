function [yp,zp,mp] = deanbeachprofile(y,zBC,z1km,ubs,isplot)
%
%-------function help------------------------------------------------------
% NAME
%   deanbeachprofile.m
% PURPOSE
%   Find the bed slope across the surf zone
%   the profile is based on a user defined slope between Hw and SWL (0mOD)
%   and a Dean profile below this level. The sediment parameter for the
%   Dean profile, A, is determined from the water depth 1km seaward of SWL.      
% USAGE
%   [yp,zp,mp] = deanbeachprofile(y,zBC,z1km,ubs,isplot)
% INPUTS
%   y    - y-coordinates to define profile (m)
%   zBC  - beach crest level (mOD) 
%   z1km - bed level 1km from shore (mOD) [OR y,z coordinates]
%   ubs  - upper beach slope used to extend profile above 
%   isplot - true if plot is required
% OUTPUTS
%   yp - profile y -co-ordinates
%   zp - profile z -co-ordinates
%   mp - beach slope at each point on profile
%
% Author: Ian Townend
% CoastalSEA (c) April 2020
%----------------------------------------------------------------------
%
    if length(z1km)>1
        y1km = z1km(1);       %y and z have been specified
        z1km = z1km(2);
    else
        y1km = 1000;          %default distance
    end
    
    y0 = zBC*ubs;      
    A = (0-z1km)/y1km^(2/3); 
    h = A*y.^(2/3);
    m = 1./(2/3*A*y.^(-1/3));
    %now set up arrays for full profile
    yp = [0,y0+y]';     
    zp = [zBC,-h]';
    mp = [ubs,ubs,m(2:end)]';
    if isplot
      plotprofile(yp,zp,mp)
    end
end
%%
function plotprofile(yp,zp,mp)
    %plot beach profile
    figure('Name','Profile','Tag','PlotFig');
    subplot(2,1,1)
    plot(yp,zp)
    title('Profile')
    subplot(2,1,2)
    plot(yp,mp)
    title('Slope 1:bs')
end