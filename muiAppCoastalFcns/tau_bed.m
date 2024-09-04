function [tau,Cd] = tau_bed(d,d50,visc,rhow,Uc,H,Tp,phid)
%
%-------header-------------------------------------------------------------
% NAME
%   taubed.m
% PURPOSE
%   Bed shear stress under combined wave-current action
% USAGE
%   tcr = taubed(d,visc,rhow,Uc,H,Tp,phi)
% INPUTS
% d     = depth (m)
% d50   = median grain size (m)
% visc  = kinematic viscosity (m^2s^-1)
% rhow  = density of water (kg/m^3)
% Uc    = depth averaged current speed (m.s^-1)
% H     = wave height(m) - usually rms value
% Tp    = peak wave period (s)
% phid  = angle between flow and waves (deg)
% RESULTS
%   tau - struct containing:
%       tauc - current induced bed shear stress (N/m2)
%       tauw - wave induced bed shear stress amplitude (N/m2)
%       taum - mean wave-current induced bed shear stress (N/m2)
%       taux - maximum wave-current induced bed shear stress (N/m2)
%       taur - root-mean-square wave-current induced bed shear stress (N/m2)
% NOTES
%   Based on the Soulsby algorithm, TR137, 2005
%
%   It can also be used for rippled sand provided that a suitable ripple 
%   scale is used instead of d50. For ripples of height, D and wavelength L, 
%   an equivalent d50 grainsize of approximately d50 = 12*D^2/L can be used. 
%   For rippled beds the resulting shear-stresses are total stresses 
%   (including form drag of the ripples). For flat beds the shear stresses 
%   represent the skin-friction.
%
%   Regarding the choice of wave parameters, the theoretical expressions have
%   largely been derived for monochromatic sea states.  Consequently Soulsby 
%   (Marine Sands, p69) suggests that Hrms and Tp are the most appropriate 
%   wave parameters to use.
%
%   Assumes wave parameters H,Tp are valid (eg within acceptable steepness range).
%   Inputs can be vectors of Hs,Tp,phid,d,d50. All vectors must be same length
% SEE ALSO
% Calls celerity.m. Used in totaltransport_model
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%
% calculate basic parameters
%
Uc   = Uc + (Uc==0)*eps;         %adjust Uc=0 to Uc=eps to avoid divide by zero
phi  = phid*pi/180;              %convert phi from degrees to radians
Rec  = Uc.*d./visc;              %current Reynolds number
%
Tp   = Tp + (Tp==0)*eps;         %adjust Tp=0 to Tp=eps to avoid divide by zero
kp   = 2*pi./celerity(Tp,d)./Tp; %wave number
uwv  = pi*H./Tp./sinh(kp.*d);    %orbital wave velocity at bed
uwv  = uwv + (uwv==0)*eps;       %adjust uwv=0 to uwv=eps to avoid divide by zero
Aw   = uwv.*Tp/2/pi;             %wave semi-orbital excursion
Rew  = uwv.*Aw./visc;            %wave Reynolds number
%
a = 0.0001615; b = 6; c = -0.08; %defined empirical fit parameters
Cds  = a*exp(b*(Rec).^c);        %current smooth turbulent drag coefficient
fws  = 0.0521*Rew.^-0.187;       %wave friction factor
%
zo   = d50/12;                   %bed roughness length
zo   = zo + (zo==0)*eps;         %adjust zo=0 to zo=eps to avoid divide by zero
Cdr  = (0.4./(log(d./zo)-1)).^2; %current rough turbulent drag coefficient
fwr  = 1.39*(Aw./zo).^-0.52;     %wave friction factor
%
taucl = 3*rhow.*visc.*Uc./d;     %laminar current shear stress
tauwl = rhow.*Rew.^-0.5.*uwv.^2; %laminar wave shear stress
%
tauws = 0.5*rhow.*fws.*uwv.^2;   %smooth turbulent wave shear stress
tauwr = 0.5*rhow.*fwr.*uwv.^2;   %rough turbulent wave shear stress
tauw  = tauwl.*(Rew <= 1.5*10^5);%check for laminar flow 
tauw  = tauw + max(tauws, tauwr).*(Rew > 1.5*10^5);  %wave induced bed shear stress
%
taucs = rhow.*Cds.*Uc.^2;        %smooth turbulent current shear stress
taucr = rhow.*Cdr.*Uc.^2;        %rough turbulent current shear stress
tauc  = taucl.*(Rec <= 2000);    %check for laminar flow                     
tauc  = tauc + max(taucs, taucr).*(Rec > 2000); %current induced bed shear stress
%
[Cdms,Cdxs] = Cdmsxs;            %smooth wave-current drag coefficients
[Cdmr,Cdxr] = Cdmrxr;            %rough wave-current drag coefficients
%
taums = rhow.*Cdms.*Uc.^2;       %smooth mean shear stress
tauxs = rhow.*Cdxs.*Uc.^2;       %smooth maximum shear stress
taumr = rhow.*Cdmr.*Uc.^2;       %rough mean shear stress
tauxr = rhow.*Cdxr.*Uc.^2;       %rough maximum shear stress
%
% determine flow regime
taum  = zeros(size(d));
taux  = zeros(size(d));
%
idx = find(Uc==0 & H==0);      %no flow case
taum(idx) = 0;
taux(idx) = 0;
tauw(idx) = 0;
%
idx = find(Uc>0 & H==0);       %current flow only case
taum(idx) = tauc(idx);
taux(idx) = taum(idx);
tauw(idx) = 0;
%
idx = find(Uc==0 & H>0);       %wave only case
taum(idx) = 0;
taux(idx) = tauw(idx);
%
Reccr = 2000+(5.92*10^5.*Rew).^0.35;
Rewcr = 1.5*10^5;
idx = find(Rec<=Reccr & Rew<=Rewcr);             %laminar wave-current case
taum(idx) = taucl(idx);
taux(idx) = (taum(idx)+tauw(idx).*abs(cos(phi(idx)))).^2;
taux(idx) = (taux(idx)+ (tauw(idx).*abs(sin(phi(idx)))).^2).^0.5;
%
idx = find((Rec>Reccr | Rew>Rewcr) & tauxr<=tauxs);%smooth turbulent wave-current case
taum(idx) = taums(idx);
taux(idx) = tauxs(idx);
Cdm(idx) = Cdms(idx); Cdx(idx) = Cdxs(idx);
idx = find((Rec>Reccr | Rew>Rewcr) & tauxr>tauxs); %rough turbulent wave-current case
taum(idx) = taumr(idx);
taux(idx) = tauxr(idx);
Cdx(idx) = Cdmr(idx); Cdx(idx) = Cdxr(idx);
%
% Calculate the root-mean-square shear stress
taur = (taum.^2+0.5*tauw.^2).^0.5;

tau = struct('tauc',tauc,'tauw',tauw,'taum',taum,'taux',taux,'taur',taur);
Cd = struct('Cdm',Cdm,'Cdx',Cdx); 
%
%---nested functions ------------------------------------------------------
%
    function [Cdms,Cdxs] = Cdmsxs   %obtain Cdm and Cdmax for smooth turbulent case
        as   = 0.22;  %empirical coefficient
        T1   = 9*as.*Rew.*(fws/2).^0.5.*(Cds.^2.*(Uc./uwv).^4+(fws/2).^2).^0.25;
        T2   = Rec./Rew.*uwv./Uc./as.*(2./fws).^0.5;
        T3   = (Cds.^2+(fws/2).^2.*(uwv./Uc).^4).^0.25;
        A1   = T3.*(log(T2)-1)/2./log(T1);
        A2   = 0.4*T3./log(T1);
        Cdms = ((A1.^2+A2).^0.5-A1).^2;
        Cdxs = (Cdms + T3.*uwv./Uc.*(fws/2).^0.5.*abs(cos(phi))).^2;
        Cdxs = (Cdxs + (T3.*uwv./Uc.*(fws/2).^0.5.*abs(sin(phi))).^2).^0.5;
    end
%
    function [Cdmr,Cdxr] = Cdmrxr   %obtain Cdm and Cdmax for rough turbulent case
        ar   = 0.26;  %empirical coefficient
        T1   = max(ar*(fwr/2).^0.5.*Aw./zo,12);
        T2   = d./T1./zo;
        T3   = (Cdr.^2+(fwr/2).^2.*(uwv./Uc).^4).^0.25;
        A1   = T3.*(log(T2)-1)/2./log(T1);
        A2   = 0.4*T3./log(T1);
        Cdmr = ((A1.^2+A2).^0.5-A1).^2;
        Cdxr = (Cdmr + T3.*uwv./Uc.*(fwr/2).^0.5.*abs(cos(phi))).^2;
        Cdxr = (Cdxr + (T3.*uwv./Uc.*(fwr/2).^0.5.*abs(sin(phi))).^2).^0.5;
    end
%
%--------------------------------------------------------------------------
%
end
    

