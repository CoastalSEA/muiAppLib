function Qx = xshore_bailard(Hs,Tp,Dir,theta,bs,d50,g,rhow,rhos,visc)
%
%-------function help------------------------------------------------------
% NAME
%   xshore_bailard.m
% PURPOSE
%   Computes the cross-shore transport for given wave and beach conditions
% USAGE
%   cst = xshore_bailard(Hs,Tp,Dir,theta,bs,d50,g,rhow,rhos,visc)
% INPUTS
%   Hs   - significant wave height (m)
%   Tp   - peak period (s)
%   Dir  - wave direction (degTN)
%   theta- angle of shoreline fron north (degrees TN)
%   bs   - beach slope (m=1:bs)
%   d50  - grain size d50 (m)
%   g    - acceleration due to gravity(m/s2)
%   rhow - density of water (default = 1025 kg/m^3)
%   rhos - density of sediment (default = 2650 kg/m^3)
%   visc - viscosity of water (m^2/s)
% OUTPUT
%   Qx  - cross-shore volumetric transport rate (m^3/s)
% NOTES
%   Based on Bailard & Inman, 1981, Bailard, 1981,1982
% NEEDS UPDATING SEE
%   Gallagher E L, Elgar S and Guza R T, 1998, JGR
% SEE ALSO
%   littoraldrift.m, drag_coefficient.m, wave_friction.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
%assumed constants (see Bailard, 1981
etab = 0.21*100;             %bed load efficiency with dimension correction
etas = 0.025*10;             %suspended load efficiency with dimension corrction
tanphi = 0.625;              %assumed internal angle of friction for sand of 32 deg

% empirical velocity relationships from Bailard, 1982 (units of m/s)
psi1 = 0.00303-0.00144*Hs;   %wave velocity skewness parameter
psi2 = 0.00603-0.00510*Hs;   %wave velocituy skewness parameter
delu = -0.00157*Hs;          %normalised onshore current
um   = 0.319+0.403*Hs;       %near bed oscillatory velocity
u3st = 0.00548+0.000733*Hs;  %normalised velocity of u3
u5st = 0.015+0.00346*Hs;     %normalsised velocity of u5
Cf = 0.005;                  %Bailard's assumed friction value
% empirical velocity relationships from Bailard, 1982 (units of cm/s)
% Hs = Hs*100;
% psi1 = 0.303-0.00144*Hs;   %wave velocity skewness parameter
% psi2 = 0.603-0.00510*Hs;   %wave velocituy skewness parameter
% delu = 0.45-0.00157*Hs;          %normalised onshore current
% um   = 31.9+0.403*Hs;       %near bed oscillatory velocity
% u3st = 0.548+0.000733*Hs;  %normalised velocity of u3
% u5st = 1.5+0.00346*Hs;     %normalsised velocity of u5
% Cf = 0.005;                  %Bailard's assumed friction value

% alternaive formulations considered
% coeffs = drag_coefficient(1,dep,d50,visc); %indicative drag coefficient
% Cf = coeffs.Cd;                            %lower than Cf based on 1m/s
% fw = wave_friction(Hs,Tp,dep,d50,g,visc);  %higher than Cf
% for the test values Urms was lower than um
% um = rmswaveorbitalvelocity(Hs,Tp,1,g); %Soulsby and Smallman, 1986
% fprintf(' Cd=%f, fw=%f, Urms=%f, um=%f\n',[Cd;fw';Urms;um])

ws = settling_velocity(d50,g,rhow,rhos,visc)*2;
% ws = ws*100; %cm/s
% ws=4;
beta = 1./bs; %assumes direction has no effect
% [alpi,quad] = getalp(Dir,theta);
% alpi = 20*pi()/180;

%Bailard, 1982, TN 1649 eqn 23
term1 = rhow*Cf.*um.^3;
term2 = etab/tanphi.*(psi1+3/2*delu-tan(beta)/tanphi.*u3st);
term3 = um/ws*etas.*(psi2+delu.*u3st-um/ws.*etas.*u5st.*tan(beta));
ix = term1.*(term2+term3);        %immersed weight transport rate (N/m)
Qx = ix/(rhos-rhow)/g/0.6;        %volumetric transport rate (m3/m/s)

