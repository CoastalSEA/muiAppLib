function Q = littoraldrift(Hsi,Tp,Diri,depi,theta,mb,d50,Kc,g,rhos,rhow,visc)
%
%-------function help------------------------------------------------------
% NAME
%   drift.m
% PURPOSE
%   sediment transport drift rates for sand and shingle
% USAGE
%   Q = drift(Hsi,Tp,Diri,depi,theta,mb,d50,Kc,dt)
% INPUTS
%   Hsi  - inshore significant wave height (m)
%   Tp   - peak wave period (s)
%   Diri - wave direction (degrees TN)
%   depi - inshore water depth (m)
%   theta- angle of shoreline fron north (degrees TN)
%   mb   - beach slope, CoatalTools uses beach slope=1:mb
%   d50  - grain size d50 (m)
%   Kc   - drift coefficient (-) Kc=0.0006~1/(rhos-rhow), Inman using
%   Hrms, Kc=0.77/(rhois-rhow) and SPM using Hs, Kc=0.39/(rhos-rhow)=0.0002
%   g    - acceleration due to gravity (m/s2)
%   rhos - density of sediment (kg/m3)
%   rhow - density of water (kg/m3)
%   visc - viscosity of water (m2/s)
% RESULTS
%   Q - 4 column matrix with 4 estimates of drift in m3/s
%          (i) original CERC formula (SPM, 1984)
%         (ii) SANDS formulae - Dynamics of Marine Sands, Soulsby
%        (iii) Kamphuis formula
%         (iv) Damgaard & Soulsby (shingle)
% SEE ALSO
% getalp.m, refraction.m beachtransportratio.m, xshore_bailard.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%--------------------------------------------------------------------------
% 

% original code replaced with code now used in refaction calling getalp
% rads = pi()/180;
% idoff = Hsi==0 & isnan(Diri);  %index for wave directions that are offshore
%                                %refraction assigns Hsi=NaN to gaps and 0 to
%                                %offshore. Use index to set qs=0 for waves 
%                                %in an offshore direction   
% beta = (Diri-theta)*rads;      %angle between wave direction and contour
% idx = find(beta<0);            %index of negative values
% beta(idx) = -(beta(idx)+pi()); %possible angles when pi<theta<2.pi
% bsgn = beta-pi()/2;  %negative when beta<90 == waves right to left when looking at shore from sea
% alpi = bsgn.*sign(bsgn);       %angle between wave crest and bed contour (rads)

% alpi - angle between wave crest and contour (rads)
% quad - quadrant of the wave direction relative to the contour
[alpi,quad] = getalp(Diri,theta);
bsgn = zeros(size(quad));
bsgn(quad==1 | quad==4) = -1;       %waves right to left when looking at shore from sea
bsgn(quad==2 | quad==3) = +1;       %waves left to right when looking at shore from sea
idoff = quad==3 | quad==4;          %offshore waves
alpi(idoff) = NaN;

ci = celerity(Tp,depi);             %wave celerity
pidLi = 4*pi()*depi./(ci.*Tp);      %factor 4pi.d/(c.Tp)
cgi = ci/2.*(1+pidLi./sinh(pidLi)); %group celerity

% Original CERC formula (SPM, 1994)
Q(:,1) = (Kc*sign(bsgn).*(Hsi.^2.*cgi).*sin(2*alpi));

% SANDS formulae - Dynamics of Marine Sands, Soulsby
sm1 = rhos/rhow-1;   %relative density minus 1
% CERC with adjustments for wave breaking and using cg=sqrt(gh)
% for fine sediments with d50<0.6mm
Q(:,2) = (0.023*sign(bsgn)*sqrt(g).*(Hsi.^(5/2)).*sin(2*alpi)/sm1);

% Kamphuis formula
Q(:,3) = (7.3*sign(bsgn).*(Hsi.^2).*(Tp.^1.5).*(sin(2*alpi).^(0.6)).*((1./mb).^0.75)*d50^-0.25);

% Damgaard & Soulsby (shingle)
Q(:,4) = shingle_transport(Hsi,Tp,alpi,sm1,visc,d50,1./mb,g).*sign(bsgn);

Q(idoff,:) = 0;
%
% Transport functions-----------------------------------------------------
%
    function Qls = shingle_transport(Hsi,Tp,alp1,sm1,visc,d50,mb,g)
        % shingle littoral drift based on formulation of 
        % Damgaard & Soulsby,1997: see p196-8 of Sand Manual
        % test input variables based on example in manual
        % d50=0.01; Hs=1; Tp=6; alp1=10*rads; mb=0.1;
        % default assignment overwritten based on th* conditions below
        Qls1 = zeros(size(Hsi)); %#ok<PREALL> 
        Qls2 = zeros(size(Hsi)); %#ok<PREALL>
        
        Hrms  = Hsi./sqrt(2);  %RMS wave heights
        Ds    = d50*(g*(sm1)/visc^2)^(1/3); %dimensionless grain size
        thcr  = 0.3/(1+1.2*Ds)+0.055*(1-exp(-0.02*Ds)); %threshold Shields parameter
        thcrf = 16.7*thcr*sm1*d50./Hrms./sin(2*alp1)./tan(mb); %threshold factor
        dirf  = (0.95-0.19*cos(2*alp1)).*sin(2*alp1);
        thwr  = 0.15*Hrms.^0.75/sm1/g^0.25./sqrt(Tp*d50);
        thwsf = 0.004*Hrms.^1.2/sm1^1.4/g^0.2./Tp.^0.4/d50;
        thw   = max(thwr,thwsf);
        thm   = 0.1*Hrms.*sin(2*alp1).*tan(mb)/sm1/d50;
        thmx  = sqrt((thm+thw.*sin(alp1)).^2+(thw.*cos(alp1)).^2);
        %conditional relationships for longshore transport
        Qls1=0.19*sqrt(g*tan(mb)).*(sin(2*alp1)).^1.5.*Hrms.^2.5.*(1-thcrf)/12/sm1.*(thcrf<1);
        Qls2a=0.24*g^(3/8)*d50^0.25*dirf.*Hrms.^(19/8)/12/sm1./Tp.^0.25.*(thwr>=thwsf).*(thmx>thcr);
        Qls2b=0.046*g^0.4*dirf.*Hrms.^(13/5)/12/sm1^1.2./(pi()*Tp).^0.2.*(thwr<thwsf).*(thmx>thcr);
        Qls2 = Qls2a+Qls2b;
        Qls = max(real(Qls1),real(Qls2));
    end
end

