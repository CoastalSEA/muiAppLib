function Qx = xshore_bailard(Hs,Tp,Dir,dep,theta,bs,d50,g,rhow,rhos,visc)
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
%   dep  - water depth (m)
%   theta- angle of shoreline fron north (degrees TN)
%   bs   - beach slope (m=1:bs)
%   d50  - grain size d50 (m)
%   g    - acceleration due to gravity(m/s2)
%   rhow - density of water (default = 1025 kg/m^3)
%   rhos - density of sediment (default = 2650 kg/m^3)
%   visc - viscosity of water (m^2/s)
% OUTPUT
%   Qx  - cross-shore volumetric transport rate (m^3/s) and is positive in 
%         the onshore direction. Returns a struct with the results for
%         Bailard Eq13 and Eq23 as two fields.
% NOTES
%   Based on Bailard & Inman, 1981, Bailard, 1981,1982 (see TN1649,'82 for Eqs used)
%   Gallagher E L, Elgar S and Guza R T, 1998, JGR
%   Thornton E B, Humiston R T and Birkemeier W, 1996, JGR: Oceans
% SEE ALSO
%   littoraldrift.m, drag_coefficient.m, wave_friction.m
% EXAMPLE - comparison with the results given by Bailard
%   % Bailard, 1982 data set Volumes m3/m/day, Tables 1 and 2
%   Hs = [0.55,0.65,0.59,1.01,0.99,1.4,0.91,0.62,0.73];
%   bs = 50; %beach face is 1:20 reducing to 1:50 in surf zone
%   d50 = 0.17/1000;
%   Qx = xshore_bailard(Hs,[],[],[],[],bs,d50,g,rhow,rhos,visc);
%   Vx = Qx*3600*24;      %convert m3/m/s to m3/m/d
% 
%   % check plots against Bailard's data - gives a similar order but not the same.
%   Vp = [0.72,0.02,0.52,-2.41,-3.74,-3.9,0.36,0.79,0.34];
%   Vm = [0.21,-8.2,9.3,1.08,1.08,-13.4,0.05,0.78,-1.99];
%   figure('Tag','PlotFig');
%   scatter(Vm,Vp)
%   hold on
%   scatter(Vm,Vx)
%   plot(xlim, 0*[1 1],'--k');
%   plot(0*[1 1], ylim,'--k');
%   hold off
%   xlabel('Measured volume change (m^3/m/day)')
%   ylabel('Predicted volume change (m^3/m/day)')
%   title('Comparison of Bailard''s measured and predicted volume change')
%   legend('Bailard model','Code');
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
%assumed constants (see Bailard, 1982)
etab = 0.1;                  %bed load efficiency
etas = 0.02;                 %suspended load efficiency
tanphi = 0.625;              %assumed internal angle of friction for sand of 32 deg

% empirical velocity relationships from Bailard, 1982 (units of cm/s)
Hs = Hs*100;
psi1 = 0.303-0.00144*Hs;          %wave velocity skewness parameter
psi2 = 0.603-0.00510*Hs;          %wave velocituy skewness parameter
delu = 0-0.00157*Hs;              %normalised onshore current (based on Fig 8 and not Eq.24)
um   = 31.9+0.403*Hs;             %near bed oscillatory velocity (cm/s)
u3st = 0.548+0.000733*Hs;         %normalised velocity of u3
u5st = 1.5+0.00346*Hs;            %normalised velocity of u5
Cf = 0.005;                       %Bailard's assumed friction value

beta = 1./bs; %assumes direction has no effect

% % alternative formulations considered
% um = rmswaveorbitalvelocity(Hs,Tp,dep,g);   %Soulsby and Smallman, 1986
% coeffs = drag_coefficient(um,dep,d50,visc); %indicative drag coefficient
% Cf = coeffs.Cd;                             %lower than Cf based on 1m/s
% fw = wave_friction(Hs,Tp,dep,d50,g,visc);   %higher than Cf
% % for the test values Urms was lower than um
% % fprintf(' Cd=%f, fw=%f, Urms=%f, um=%f\n',[Cd;fw';Urms;um])
% 
% [alpi,~] = getalp(Dir,theta);
% % alpi = 20*pi()/180;
% um = um.*cos(alpi);

ws = settling_velocity(d50,g,rhow,rhos,visc);
ws = ws*100; %cm/s
%ws = 4;                            %Bailard's assumed value in cm/s

%Bailard, 1982, TN 1649 eqn 13 

%Bailard, 1982, TN 1649 eqn 23 when delu <<1, delv <<1, and alp2 and alp3 are small
rho = rhow/1000;
term1 = rho*Cf.*um.^3;
term2 = etab/tanphi.*(psi1+3/2*delu-tan(beta)/tanphi.*u3st);
term3 = um/ws*etas.*(psi2+delu.*u3st-um/ws.*etas.*u5st.*tan(beta));
ix = term1.*(term2+term3);        %immersed weight transport rate (dyne/cm)
Qx = ix/1e3/(rhos-rhow)/g/0.6;    %volumetric transport rate (m3/m/s)



