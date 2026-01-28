function [Q,gR2,zgR2] = otopBeach(swl,Hs,Tp,Ho,beach)
%
%-------function help------------------------------------------------------
% NAME
%   otopBeach.m
% PURPOSE
%   wave overtopping on gravel beaches
% USAGE
%   Q = otopBeach(bs,Ho,Tp);
% INPUTS
%   swl  - still water level (m above datum)
%   Hs   - significant wave height (m)
%   Tp   - offshore peak wave period (s)
%   Ho   - "effective" deepwater offshore wave height
%   beach - struct with following fields (* defined in ctWaveParameters):
%               BeachToeLevel  - toe elevation (mOD)
%               BeachCrestLevel  - crest elevation (mOD) *
%               UpperBeachSlope - upper beach slope (1:ubs) *        
% RESULTS
%   Q - overtopping discharge in m3/s/m-run of beach
%   gR2 - 2% runup (m)
%   zgR2 - Runup elevation (gR2+swl)
% NOTES
%   based on equation proposed by Stokes et al, Coastal Eng., 2021, with
%   the runup determined using gravel beach equation using Poate etal, Coastal Eng., 2016
%
% Author: Ian Townend
% CoastalSEA (c)Jan 2026
%--------------------------------------------------------------------------
%
g = 9.81;
tl = beach.BeachToeLevel;         %beach toe level (m above datum)
cl = beach.BeachCrestLevel;       %beach crest level (mOD)
ubs = beach.UpperBeachSlope;      %beach slope (1:ubs), ie tan(alpha)
beta = 1/ubs;
z1km = beach.BedLevelat1km;       %bed level 1km out from SWL (mOD)
dep = swl-tl;                     %water deoth at toe (m)

%find the foreshore slope in a water depth that is 0.78*Hs deeper than toe
fs = profileslope(dep+0.78*Hs,swl,z1km,ubs); %first argument is depth
hbflag = 3;      %specific case in hb_break for overtopping, which uses Tz
Tz = Tp*0.7775;  %assumes JONSWAP with gamma=3.3
Hb = hb_break(dep,fs,Tz,g,hbflag);
Hsb = min(Hb,Hs,'includenan');
gR2 = 0.33*sqrt(beta).*Tp.*Ho; %R2% Poate et al (gravel)
zgR2 = swl+gR2;                %elevation of runup
Rs = (zgR2-cl)./Hsb;           %relative freeboard
Qs = 8*10^-5*exp(3.1*Rs);      %normalised overtopping discharge
Qo = Qs.*sqrt(g*Hsb.^3);       %overtopping discharge (m3/s/m-run of beach)     
Q = Qo.*(Qo>1e-6);             %remove very small values