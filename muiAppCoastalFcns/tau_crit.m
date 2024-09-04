function [tau,me] = tau_crit(rhob,d50,pcm,visc,rhow,rhos)
%
%-------function help------------------------------------------------------
% NAME
% tau_crit.m
% PURPOSE
%   Calculate the critical erosion shear stress and erosion rate for sand, 
%   mud or mixed sediments
% USAGE
%   [tau me] = tau_crit(rhob,d50,visc,rhow,rhos,pcm)
% INPUTS
%   rhob - bulk density of bed (kg/m^3)
%   d50  - median grain size (m)
%   pcm  - percentage mud content (%)
%   visc - viscosity of water (m2/s)
%   rhow - density of water (default = 1025 kg/m^3)
%   rhos - density of sediment (default = 2650 kg/m^3)     
% OUTPUTS
%   tau - critical erosion shear stress (Pa)
%   me  - erosion rate (m/s)
% NOTES 
%   Set pcm=0 to use for sand, pcm=1 to use for mud. Otherwise define %mud
%
% Author: Ian Townend
% CoastalSEA (c)June 2016
%--------------------------------------------------------------------------
%
g     = 9.81;         %acceleration due to gravity(m/s^2)
%
% calculate the critical erosion threshold
if pcm == 0
    Ds  = d50*(g*(rhos/rhow-1)/visc^2)^(1/3);
    theta_cr = 0.3/(1+1.2*Ds)+0.055*(1-exp(-0.02*Ds));
    tau = g*(rhos-rhow)*d50*theta_cr;
elseif pcm == 1
    rhod = rhos*(rhob-rhow)/(rhos-rhow); %dry density(kg/m^3)
    tau = 0.185*g*(rhod/rhow)^0.75;
else
    %sand tau_crit
    Ds  = d50*(g*(rhos/rhow-1)/visc^2)^(1/3);
    theta_cr = 0.3/(1+1.2*Ds)+0.055*(1-exp(-0.02*Ds));%eq(77) in Sands Manual
    tau_s = g*(rhos-rhow)*d50*theta_cr;  %shear stress for sand sediments  
    %mud tau_crit
    rhod = rhos*(rhob-rhow)/(rhos-rhow); %dry density(kg/m^3)
    tau_m = 0.185*(rhod/rhow)^0.75;      %shear stress for mud sediments
    %mixed tau_mx value (NB constant should vary with sediment properties)
    tau_mx = 5*tau_s; 
    pcm_mx = 0.2;  %percentage at which tau_e max occurs
    if pcm <= pcm_mx
        tau = linterp(tau_s,tau_mx,pcm/pcm_mx);
    else
        tau = linterp(tau_mx,tau_m,(pcm-pcm_mx)/(1-pcm_mx));
    end
end
%
% calculate the erosion rate
me_s = 0.001; %typical value for sand Zeng(2013) (m/s)
me_m = 0.002; %typical value for mud, Mud manual p67 (m/s)
me = linterp(me_s,me_m,pcm);
%
%--------------------------------------------------------------------------
%
function yi = linterp(ymin, ymax, xi)
    x = [0 1];
    y = [ymin ymax];
    yi= interp1(x,y,xi);

