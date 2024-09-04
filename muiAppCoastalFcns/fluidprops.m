function [rhow,visc] = fluidprops(S,T)
%
%-------header-------------------------------------------------------------
% NAME
%   fluidprops.m
% PURPOSE
%   Function to calculate fluid density and kinematic viscosity based on 
%   salinity and temperature.
% USAGE
%    [rhow,visc] = fluidprops(S,T)
% INPUTS
%   S - salinity as a weight ratio e.g. g/kg or ppt
%   T - temperature (deg C)
% RESULTS
%   rhow - density of sea water (kg/m3)
%   visc - kinematic viscosity (m2/s)   
% NOTES
%   van Rijn L C, 1993, Principles of sediment transport in rivers, estuaries 
%   and coastal seas, Aqua Publications, Amsterdam - see chapter 3.
%   Whitehouse R, Soulsby R, Roberts W, Mitchener H, Dynamics of estuarine 
%   muds: a manual for practical applications, Thomas Telford, London, 2000 
%   see chapter 2.
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%

    %Chlorinity
    Cl = 0.554*(S-0.03); 
    %Density of sea water
    rhow = 1000+1.455*Cl-0.0065*(T-4+0.4*Cl)^2;
    %kinematic viscosity
    t15 = T-15;
    visc = (1.14-0.031*t15+0.00068*t15^2)*10^-6;
end
