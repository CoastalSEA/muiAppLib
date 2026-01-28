function [R2,setup,swash,gR2] = runup(bs,H,T)
%
%-------function help------------------------------------------------------
% NAME
%   runup.m
% PURPOSE
%   wave runup magnitude
% USAGE
%   [R2,zR] = runup(beta,H,T)
% INPUTS
%   beta- beach foreshore slope (1:bs)
%   H   - offshore significant wave height (m)
%   T   - offshore wave period (s) - paper is not clear whether this is peak or mean
% RESULTS
%   R2 - 2% runup (m)
%   setup - setup component of runup (m)
%   swash - swash component of runup (m)
%   gR2 - 2% runup using the eqn for gravel beaches by Poate etal (m)
% NOTES
%   based on equation proposed by Stockdon etal, Coastal Eng. 2006.
%   gravel beach equation using Poate etal, Coastal Eng.2016
%
% Author: Ian Townend
% CoastalSEA (c)June 2018
%--------------------------------------------------------------------------
%
g = 9.81;       %acceleration due to gravity
beta = 1./bs;   %beach slope as tan(alpha) or 1/bs
HoLo = H.*(g*T.^2/2/pi());
setup = 0.35*beta.*sqrt(HoLo);
swash = sqrt(HoLo.*(0.56*beta.^2+0.004))/2;
R2 = 1.1*(setup+swash);                     %R2 Stockdon et al

gR2 = 0.33*sqrt(beta).*T.*H;                %R2 Poate et al (gravel)

