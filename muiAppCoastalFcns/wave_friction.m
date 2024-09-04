function fw = wave_friction(Hs,Tp,dep,d50,g,visc)
%
%-------function help------------------------------------------------------
% NAME
%   wave_friction.m
% PURPOSE
%   compute the wave friction factor for rough and smooth turbulent
%   conditions and smooth laminar conditions
% USAGE
%   Cd = drag_coefficient()
% INPUTS
%   Hs - significant wave height (m)
%   Tp - peak period (s)
%   dep  - water depth (m)
%   d50  - grain size d50 (m)
%   g    - acceleration due to gravity(m/s2)
%   visc - viscosity of water (m2/s)
% OUTPUT
%   fw - wave friction factor
% NOTES
%   Based on Soulsby, Dynamics of Marine Sands, 
% SEE ALSO
%   rmswaveorbitalvelocity.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    %root mean square wave orbital velocity;
    Urms = rmswaveorbitalvelocity(Hs,Tp,dep,g); %Soulsby and Smallman, 1986
    Uw = sqrt(2)*Urms;                          %equivalent monochromatic wave
    A = Uw.*Tp/2/pi();                          %semi-orbital excursion, eq.58
    Rw = Uw.*A/visc;                            %Reynolds no., eq.58
%     ks = 2.5*d50;                             %Nikuradse grain roughness, eq.24
    zo = d50/12;                                %bed roughness length (m)   
    fwr = 1.39*(A./zo).^-0.52;                  %rough turbulent friction, eq.62
    fwl = 2*Rw.^-0.5;                           %smooth laminar friction, eq.63 
    fws = 0.0521*Rw.^-0.187;                    %smooth turbulent friction, eq.63
    fws(Rw<=5e5) = fwl(Rw<=5e5);
    if isrow(fwr), fws = fws'; fwr = fwr'; end
    fw = max([fws,fwr],[],2);