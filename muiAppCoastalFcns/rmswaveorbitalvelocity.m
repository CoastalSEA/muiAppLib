function Urms = rmswaveorbitalvelocity(Hs,Tp,dep,g)
%
%-------function help------------------------------------------------------
% NAME
%   rmswaveorbitalvelocity.m
% PURPOSE
%   Calculate the root mean square wave orbital velocity 
% USAGE
%   Urms = rmswaveorbitalvelocity(Hs,Tp,dep,g)
% INPUTS
%   Hs - significant wave height (m)
%   Tp - peak period (s)
%   dep  - water depth (m)
%   g    - acceleration due to gravity(m/s2)
% OUTPUT
%   Urms - root mean square orbital wave velocity (m/s)
% NOTES
%   Uses method given in HR Wallingford report SR 76 by Soulsby and Smallman 
% SEE ALSO
%   wave_friction.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    Tz = 0.7775*Tp;  %assumes a JONSWAP spectrum with gamma=3.3
    
    Tn = sqrt(dep/g);  %natural scaling period based on water depth
    
    %using eq(20)
    t = Tn./Tz;
    A = (6500+(0.56+15.54*t).^6).^(1/6);
    Urms = 0.25*Hs./(Tn.*(1+A.*t.^2).^3);
end