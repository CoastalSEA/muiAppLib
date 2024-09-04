function var = sedprops(rhob,rhos,rhow,prop,S,T,sd)
%
%-------function help------------------------------------------------------
% NAME
% sedprops.m
% PURPOSE
%   Function to return one of a range of sediment properties based on selection
%   defined in 'prop'
% USAGE
%   var = sedprops(rhob,rhos,rhow,prop,S,T,sd)
% INPUTS
%   rhob - bulk density or wet density (kg/m3)
%   rhos - density of sediment (kg/m3)
%   rhow - density of sea water (kg/m3)
%   prop - selected propery to be returned as 'var'
%             SpecificDensity (-)
%             DryDensity (kg/m3)
%             MassSedConc (kg/m3), sediment concentration by mass (saturated)     
%             VolSedConc (-), concentration by volume
%             FluidSedMixDensity (kg/m3), fluid-sediment mixture density
%             Porosity (-)
%             VoidRatio (-)
%             SolidContent (-)
%             WaterContent (-)
%             UnitWeightSed (N/m3)
%             UnitSubWeightSed (N/m3)
%             DynamicViscosity (Ns/m2)
%             WaterSedDynVisc, dynamic viscosity coefficient for water-sediment (m2/s)
%             WaterSedKinVisc, kinematic viscosity coefficient for water-sediment (Ns/m2)
%   S - salinity %weight ratio e.g. g/kg
%   T - temperature (deg C)
%   sd - degree of saturation (1=saturated, 0=dry)
% OUTPUTS
%   var - value of selected property
% SEE ALSO
%   called from sediment_properties and uses fluidprops.
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
    if (nargin<7 || isempty(sd)) && strcmp(prop,'WaterContent')
        sd = 1;  %default degree of saturation if not defined
    else
        sd = 0;
    end
    
    if nargin>5 && isempty(rhow) && isempty(S)
        S = 32;  %defaulat salinity if not defined
        T = 10;  %default temperature if not defined
        [rhow,visc] = fluidprops(S,T);
    elseif nargin>5 && isempty(S)
        S = 32;  %defaulat salinity if not defined
        T = 10;  %default temperature if not defined
        [~,visc] = fluidprops(S,T);
    elseif nargin<5
        S = 32;  %defaulat salinity if not defined
        T = 10;  %default temperature if not defined
        [~,visc] = fluidprops(S,T);
    else
        [rhow,visc] = fluidprops(S,T);
    end
    
    g = 9.81;                             %acceleration due to gravity (m/s2)
    s = rhos/rhow;                        %specific density (-)
    rhod = (rhob-rhow)/(rhos-rhow)*rhos;  %dry density (kg/m3)
    
    switch prop
        case 'SpecificDensity'
            var = s;
        case 'DryDensity'
            var = rhod;
        case 'MassSedConc'
            var = rhod;
        case 'VolSedConc'
            var = rhod/rhos;
        case 'FluidSedMixDensity'
            cv = rhod/rhos;
            var = rhow*(1+(s-1)*cv);
        case 'Porosity'
            var = (rhos-rhod)/rhos;
        case 'VoidRatio'
            p = (rhos-rhod)/rhos;
            var = p/(1-p);
        case 'SolidContent'
            var = s*rhod/(rhos+(s-1)*rhod);       
        case 'WaterContent'
            p = (rhos-rhod)/rhos;
            e = p/(1-p);
            var = e*sd*rhow/rhos;
        case 'UnitWeightSed'
            var = g*rhos;
        case 'UnitSubWeightSed'
            var = g*(rhos-rhow);
        case 'DynamicViscosity'
            var = rhow*visc;
        case 'WaterSedDynVisc'
            cv = rhod/rhos;
            lambda = ((0.74/cv)^(1/3)-1)^-1;
            eta = rhow*visc;
            var = eta*(1+lambda)*(1+0.5*lambda);
        case 'WaterSedKinVisc'
            cv = rhod/rhos;
            lambda = ((0.74/cv)^(1/3)-1)^-1;
            eta = rhow*visc;
            etam = eta*(1+lambda)*(1+0.5*lambda);
            rhom = rhow*(1+(s-1)*cv);
            var = etam/rhom; 
    end
end
       
    