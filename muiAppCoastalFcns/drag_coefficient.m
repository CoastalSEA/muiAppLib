function coeffs = drag_coefficient(U,dep,d50,visc)
%
%-------function help------------------------------------------------------
% NAME
%   drag_coefficient.m
% PURPOSE
%   Compute the drag cofficient, Cd, under steady flow. also returns the
%   bed roughness length and otehr equivalent friction coefficients
% USAGE
%   Cd = drag_coefficient(U,h,d50,visc)
% INPUTS
%   U - depth averaged current speed (m/s)
%   dep - water depth (m)
%   d50  - grain size d50 (m)
%   visc - viscosity of water (m2/s)
% OUTPUT
%   coeffs - struct containing:
%       ks - Nikuradse roughness length (m);
%       zo - bed roughness length (m)
%       Cd - drag coefficient
%       f  - Darcy-Weisbach resistance coefficient
%       Ch - Chezy coefficient
%       n  - Manning coefficient
% NOTES
%   Based on Soulsby, Dynamics of Marine Sands, 
% SEE ALSO
%   settling_velocity.m, wave_friction.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
ks = 2.5*d50;  %Nikuradse roughness length (m), eq.24
ustar = U/7.*(d50./dep).^(1/7); %friction velocity (m/s), eq.34 
zo = ks/30*(1-exp(-ustar*ks/27/visc))+visc/9./ustar;  %bed roughness length (m), eq.23a
Cd = (0.4./(1+log(zo./dep))).^2;      %drag coefficient, eq.36
%eq31 give following conversions to other friction factors
f = 8*Cd;                       %Darcy-Weisbach resistance coefficient
Ch = sqrt(9.81./Cd);            %Chezy coefficient
n = sqrt(dep.^(1/3).*Cd/9.81);  %Manning coefficient
coeffs = struct('ks',ks,'zo',zo,'Cd',Cd,'f',f,'Ch',Ch,'n',n);