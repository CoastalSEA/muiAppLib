function [U,Cd] = ucrit(d,d50,rhow,taucr,Hs,Tp,wcflag)
%
%-------function help------------------------------------------------------
% NAME
%   ucrit.m
% PURPOSE
%   compute the critical flow velocity for a given critical shear stress and
%   wave conditions in the wave-currrent case
% USAGE
%       [U,Cd] = ucrit(d,d50,rhow,taucr,Hs,Tp,wcflag)
% INPUTS
%   Parameters for current only case
%   d     - depth of water, m
%   d50   - median sediment grain size, m
%   rhow  - density of water, kgm^-3
%   taucr - critical shear stress, Pa
%   Additional parameters for wave-current case
%   Hs    - significant wave height (m)
%   Tp    - peak wave period (s)
%   wcflg - flag for current or wave-current estimate (default is 0)
% OUTPUT
%   U threshold flow velocity, taking account of waves if included  
%   Cd drag coefficient for defined conditions
% NOTES
%   Based on the inverse of the Soulsby algorithm, TR137, 2005
%   Assumes waves and currents are aligned so that phi=0
%   Use tma_spectrum.m to get Hs and Tp if only wind data is available
% SEE ALSO
%   Calls celerity.m and lambertw.m. All other funtions are internal
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%--------------------------------------------------------------------------
% check inputs
if nargin<5
    Hs = 0; Tp = 0; wcflag = 0;
elseif nargin<7
    wcflag = 0;
end
%
intflag = check_vector_lengths(d,Hs,Tp);
if isempty(intflag)
    warndlg('Vector inputs are a different length');
    U =[]; Cd=[];
    return;
end
%
visc = 1.34e-6;      %kinematic viscosity, m^2s^-1
Uthreshold = 0.05;   %depth above which currents are calculated, m
Hthreshold = 0.1;    %depth above which waves are taken into account, m
%
d = d.*(d>=Uthreshold);
ucl = d/3/rhow/visc*taucr;
Rel = ucl.*d/visc;
U = zeros(size(d));
if wcflag==0 %current only case (sub-functions work with vectors)
    U = ucl.*(Rel<=2000);
    [ucs,Cds] = Ucs(d,rhow,taucr);
    [ucr,Cdr] = Ucr(d,d50,rhow,taucr);
    U = U + min(ucs,ucr).*(Rel>2000);
    if ucr<ucs, Cd = Cdr; else, Cd = Cds; end
else  %wave-current case (sub-functions do not work with vectors)
    % obtain wave height
    for idx = 1:length(d)
        if d(idx) > Hthreshold && Hs(idx)>0
            Hrms = Hs(idx)/sqrt(2);   %Root mean square wave height
            kp   = 2*pi/celerity(Tp(idx),d(idx))/Tp(idx);
            uwv  = pi*Hrms/Tp(idx)/sinh(kp*d(idx));
            Aw   = uwv*Tp(idx)/2/pi;
            Rew  = uwv*Aw/visc;
            Rec  = 2000 + (5.92*10^5*Rew)^0.35;
            if Rel(idx)<=Rec && Rew<=1.5*10^5
                U(idx) =ucl;
                Cd = 0.0001615*exp(6*Rel^-0.08);
            else
                [uwcs,~,Cdms] = Uwcs(d(idx),visc,rhow,taucr,uwv,Rew);
                [uwcr,~,Cdmr] = Uwcr(d(idx),d50,rhow,taucr,uwv,Aw);
                if d50==0
                    U(idx) = uwcs;
                    Cd = Cdms;
                else
                    U(idx) = min(uwcs,uwcr);
                    if ucr<ucs, Cd = Cdmr; else, Cd = Cdms; end
                end
            end
        else
            [ucs,Cds]  = Ucs(d(idx),rhow,taucr);
            [ucr,Cdr]  = Ucr(d(idx),d50,rhow,taucr);
            if d50==0
                U(idx) = ucs;
                Cd = Cds;
            else
                U(idx) = min(ucs,ucr);
                if ucr<ucs, Cd = Cdr; else, Cd = Cds; end
            end
        end
    end
end
%
%--------------------------------------------------------------------------
% Current velocity for the smooth turbulent current only case
function [ucs,cds] = Ucs(d,rhow,taucr)
    visc = 1.34e-6;
    a    = 0.0001615; b=6; c=-0.08;
    fact = taucr/rhow/a.*(d/visc).^2;
    A    = b*c/2*(fact).^(c/2);
    LW   = lambertw(A);
%     ucs  = sqrt(((visc./d).^2).*exp(log(fact)-2*LW/c)).*(d>0);
%     for idx=find(isnan(ucs)) ucs(idx)=0; end
    % the above has a divide by zero and a loop, better (and quicker) to use the following:
    ucs = zeros(size(d)); cds = ucs;
    ind = d>0;
    ucs(ind)  = sqrt(((visc./d(ind)).^2).*exp(log(fact(ind))-2*LW(ind)/c));
    cds(ind)  = a*exp(b*(ucs(ind).*d(ind)/visc).^c);
%
%--------------------------------------------------------------------------
% Current velocity for the rough turbulent current only case
function [ucr,cdr] = Ucr(d,d50,rhow,taucr)
    zo  = d50/12;
%     ucr = sqrt(taucr/rhow./Cd).*(d>0);
%     for idx=find(isnan(ucr)) ucr(idx)=0; end   
    % again, the above has a divide by zero and a loop, better (and quicker) to use the following:
    cdr  = zeros(size(d));
    ucr = zeros(size(d));
    ind = d>0;
    cdr(ind)  = (0.4./(log(d(ind)/zo)-1)).^2;
    ucr(ind) = sqrt(taucr/rhow./cdr(ind));
%
%--------------------------------------------------------------------------
% Wave-current velocity for the smooth turbulent case
function [uwcs,Cdws,Cdm] = Uwcs(d,visc,rhow,taucr,uwv,Rew)
   as   = 0.24;
   fws  = 0.0521*Rew.^-0.187;
   tauw = 0.5*rhow*fws*uwv^2;
   %check whether the wave threshold already exceeds taucr
   if tauw > sqrt(2)*taucr
       warning('Wave shear stress exceeds threshold')
       uwcs = 0; Cdws = 0; Cdm = 0;
       return
   else
       taurm = sqrt(taucr^2-0.5*tauw^2);
   end
   uwcs  = 1; %initialisation value
   delU = 1;
   while delU > 0.001
       Uc   = uwcs;
       Rec  = Uc*d/visc;
       Cdws  = 0.0001615.*exp(6*Rec.^-0.08);
       T1   = 9*as*Rew*sqrt(fws/2)*(Cdws^2*(Uc/uwv)^4+(fws/2)^2)^0.25;
       T2   = Rec/Rew*uwv/Uc/as*sqrt(2/fws);
       T3   = (Cdws^2+(fws/2)^2*(uwv/Uc)^4)^0.25;
       A1   = T3*(log(T2)-1)/2/log(T1);
       A2   = 0.4*T3/log(T1);
       Cdm  = (sqrt(A1^2+A2)-A1)^2;
       uwcs  = sqrt(taurm/rhow/Cdm);
       delU = abs(uwcs-Uc);
   end
%
%--------------------------------------------------------------------------
% Wave-current velocity for the rough turbulent case
function [uwcr,Cdwr,Cdm] = Uwcr(d,d50,rhow,taucr,uwv,Aw)
   ar   = 0.24;
   zo   = d50/12;
   fwr  = 1.39*(Aw/zo)^-0.52;
   tauw = 0.5*rhow*fwr*uwv^2;
   %check whether the wave threshold already exceeds taucr
   if tauw > sqrt(2)*taucr
       warning('Wave shear stress exceeds threshold')
       uwcr = 0;  Cdwr = 0; Cdm = 0;
       return
   else
       taurm = sqrt(taucr^2-0.5*tauw^2);
   end
   Cdwr  = (0.4/(log(d/zo)-1))^2;
   T1   = max(ar*sqrt(fwr/2)*Aw/zo,12);
   T2   = d/T1/zo;  
   uwcr  = 1; %initialisation value
   delU = 1;
   while delU > 0.001
       Uc   = uwcr;
       T3   = (Cdwr^2+(fwr/2)^2*(uwv/Uc)^4)^0.25;
       A1   = T3*(log(T2)-1)/2/log(T1);
       A2   = 0.4*T3/log(T1);
       Cdm  = (sqrt(A1^2+A2)-A1)^2;
       uwcr  = sqrt(taurm/rhow/Cdm);
       delU = abs(uwcr-Uc);
   end