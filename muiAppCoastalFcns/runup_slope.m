function Rslope = runup_slope(H0,Tp,ubs,C,d50,S,phi,g,visc)
%
%-------function help------------------------------------------------------
% NAME
%   runup_slope.m
% PURPOSE
%   Calculate runup beach slope using Reis A H and Gama C, 2010 
% USAGE
%   slope = runup_slope(H,T,ubs,C,d50,S,phi,g,visc)
% INPUTS
%   H0  - wave height (m)  NB should be deepwater value
%   Tp  - wave period (s)
%   ubs - upper beach slope (1:ubs) [tested using average around mtl]
%   C   - Setup coefficient (-)
%   d50 - sediment grain size (m)
%   S   - sphericity of sediment (-)
%   phi - beach porosity (-)
%   g   - gravity (m/s^2)
%   visc- viscosity (m2/s)
% OUTPUTS
%   Rslope - slope of the beach in the runup zone: rs = 1/tan(theta)
% NOTES
%   Reis A H and Gama C, 2010, Sand size versus beachface slope — An 
%   explanation based on the Constructal Law. Geomorphology, 114 (3), 276-283
%   Eq(20) is implicit because iribarren number depends on the slope
%   Eq(15a) is a = Iri0.*H0*S^2*phi^3*g/(150*visc*(1-phi^2))
%   Eq(20) is beta = 1/Rslope = (a/2).^(2/3).*B.^(-4/3)*d50^(4/3)
%   modified here to be explicit in beta
%   Eq(15a) becomes a = H0*S^2*phi^3*g/(150*visc*(1-phi^2))/sqrt(H0/L0)
%   Eq(20) becomes beta = (a/2).^2.*B.^-4*d50^4
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    n = 0.015+0.2*(d50-0.0003);         %Manning's n
    [Iri0,~] = iribarren(H0,Tp,ubs,g);  %deepwater iribarren number
    a = Iri0.*H0*S^2*phi^3*g/(150*visc*(1-phi^2));
    [~,setup,~] = runup(ubs,H0,Tp);     %wave setup
    h = C*setup;                        %estimate of flow depth    
    B = (h.^(5/2))/n;                   %Eq(17)
    beta = (a/2).^(2/3).*B.^(-4/3)*d50^(4/3);
    Rslope = 1./beta;                   %slope is 1/tan(theta)
end