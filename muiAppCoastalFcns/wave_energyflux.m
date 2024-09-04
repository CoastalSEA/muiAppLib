function Ef = wave_energyflux(Hs,Tp,wdep,rhow,g)
%
%-------function help------------------------------------------------------
% NAME
%   wave_energyflux.m
% PURPOSE
%   function to calculate the wave energy flux using linear wave theory
% USAGE
%   Ef = wave_energyflux(swl,Hs,Tp,bs,bl)
% INPUTS
%   Hs  - incident significant wave height (m)
%   Tp  - peak wave period (s)
%   wdep - water depth (m)
%   rhow- water density(kg/m3)
%   g - acceleration due to gravity(m/s2)
% RESULTS
%   Ef - wave energy flux(J/ms)
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%--------------------------------------------------------------------------
%                    
en = rhow*g*Hs.^2/8;                %wave energy per unit wave crest(J/m2)
cel = celerity(Tp,wdep);            %wave celerity (uses Hunt eqn)
pidL0 = 4*pi()*wdep./(cel.*Tp);     
cg = cel/2.*(1+pidL0./sinh(pidL0)); %group celerity (linear wave theory)
Ef = en.*cg;                        %energy flux(J/ms)