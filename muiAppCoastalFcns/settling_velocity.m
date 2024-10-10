function ws = settling_velocity(d50,g,rhow,rhos,visc,rhoc)
%
%-------function help------------------------------------------------------
% NAME
%   settling_velocity.m
% PURPOSE
%   calculate the settling velocity using Soulsby equation in the Mud
%   Manual 
% INPUT
%   d50  = median sediment grain size (m)
%   rhoc = density of sediment(kg/m^3) - default is 0.
%   g = acceleration due to gravity (m/s2)
%   rhow = density of water (kg/m3)
%   rhos = density of sediment (kg/m3)
%   visc = kinematic viscosity of water (m2/s)
%   [the above parameters can also be passed as a struct]
% OUTPUT
%   ws - settling velocity (m/s)
% NOTES
%   Whitehouse R J S, Soulsby R L, Roberts W and Mitchener H, 2000, 
%   Dynamics of estuarine muds, Thomas Telford, London.
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%--------------------------------------------------------------------------
%
    if nargin==1 && isstruct(d50)
        [d50,rhow,rhos,rhoc,visc,g] = getInput(d50);  %unpack struct input
    elseif nargin<6
        rhoc=0;
    end
    %
    conc  = rhoc/rhos; %volume concentration (-)
    %
    % calculate settling velocity.  Mud Manual, eqn 5.7 neglecting flocculation
    Ds   = d50*(g*(rhos/rhow-1)/visc^2)^(1/3);
    ws   = visc/d50*(sqrt(107.3+1.049*((1-conc)^4.7)*Ds^3)-10.36);
end
%%
function [d50,rhow,rhos,rhoc,visc,g] = getInput(inp)
    %unpack struct of input parameters
    d50 = inp.d50;
    rhow = inp.rhow;
    rhos = inp.rhos;
    rhoc = inp.rhoc;
    visc = inp.visc;    
    g =  inp.g;
end