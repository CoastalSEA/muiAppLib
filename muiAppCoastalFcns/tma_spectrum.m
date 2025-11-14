function [Hs,Tp,T2] = tma_spectrum(Uw,zw,Fch,df,ds)
%
%-------function help------------------------------------------------------
% NAME
%   tma_spectrum.m
% PURPOSE
%   Calculate the TMA spectrum for waves that are depth limited
% USAGE
%   [Hs Tp T2] = tma_spectrum(Uw,zw,Fch,df,ds)
% INPUTS
%   Uw - wind speed (m/s) <can be vector>
%   zw - elevation of wind speed measurement (m)
%   Fch- dominant fetch length (m) <can be vector>
%   df - average water depth over fetch (m)
%   ds - depths at site (m) > can be a vector <can be vector>
%           NB vector data must be the same length
% RESULTS
%   Hs - significant wave height (m)
%   Tp - peak wave period (s)
%   T2 - spectral zero upcross period or mean sero crossing wave period (s)
% NOTES
%   CERC-84-7 (p10-12) 
%   Bouws et al, 1985,
%   Donelan, Hamilton and Hui, R.Soc, 1985, eqn 6.1)
%   Hunt,ASCE,WW4,1974,p457-459
% SEE ALSO
%   uses celerity.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    %supress integral warnings
    wid = 'MATLAB:integral:NonFiniteValue';
    warning('off',wid)
    
    g = 9.81;
    % Adjust wind speed to 10m using power law profile
    U = Uw*(10/zw)^(1/7);
    
    U(U<0.1) = NaN;       %wind speed too small to calculate wave height
    %sort out which input variables are vector
    intflag = check_vector_lengths(U,Fch,ds);
    if isempty(intflag)
        warndlg('Vector inputs are a different length');
        Hs=[]; Tp=[]; T2=[];
        return;
    end
    
    % Calculate the peak period using the method of Donelan et al,1985
    Tp = 0.54*g^-0.77*U.^0.54.*Fch.^0.23;
    fp = 1./Tp;
    Lp = celerity(Tp,df).*Tp; %use Hunt equation to obtain celerity
    kp = 2*pi*U.^2./g./Lp;
    F0 = @(f) kit_limit(f,ds).*jonswap(f,fp,kp);
    m0 = integral(F0,min(fp/10),max(10*fp),'ArrayValued',intflag);
    F2 = @(f) (f.^2).*kit_limit(f,ds).*jonswap(f,fp,kp);
    m2 = integral(F2,min(fp/10),max(10*fp),'ArrayValued',intflag);
    Hs = 4*sqrt(m0);
    T2 = sqrt(m0./m2);
    
    if length(Hs)~=length(Tp)  
        %if depth varies Hs can be vector when Tp is scalar
        Tp = ones(size(Hs))*Tp;
    end
    
    %restore warnings
    warning('on',wid)
end
% 
%--------------------------------------------------------------------------
%
function E = jonswap(f,fp,kp)
% USAGE
% E = jonswap(f,fp,kp)
% 
% PURPOSE
% Calculate the JONSWAP spectrum energy density at frequency f
% 
% INPUTS
% f  - wave frequency (1/s)
% fp - peak frequency (1/s)
% kp - dimensionless wave number corresponding to the peak frequency of the
% TMA spectrum and the depth (-)
%
% RESULTS
% E - energy density in the Jonswap spectrum at frequency f
%
    g = 9.81;
    alp = 0.0078*kp.^0.49;
    gam = 2.47*kp.^0.39;    
    sig = 0.07.*(f<=fp);        % if f<=fp (vectorised)
    sig = sig + 0.09.*(f>fp);   % else
    E = zeros(length(fp),length(f));
    for i=1:length(f)
        q = exp((-(f(i)-fp).^2)./(2*sig(:,i).^2.*fp.^2));
        E(:,i) = g^2*(2*pi)^-4*alp.*(f(i).^-5).*(exp(-1.25*(f(i)./fp).^-4)).*(gam.^q);
    end
end
% 
%--------------------------------------------------------------------------
% 
function phi=kit_limit(f,d)
% USAGE
% phi=kit_limit(f,d)
% 
% PURPOSE
% Calculate the Kitaigorodskii limit to the spectrum at frquency f
% 
% INPUTS
% f - wave frequency (1/s)
% d - water depth (m)
%
% RESULTS
% phi - frquency dependent Kitaigorodskii adjustment to be applied to the 
% JONSWAP spectrum to take accound of depth limiting effects
%
    g = 9.81;
    omega = 2*pi*f.*sqrt(d/g);
    phi = 0.5*omega.^2.*(omega<=1);   % if omega<=1 (vectorised)
    phi = phi + 1.*(omega>2);         % else
    phi = phi + (1 - 0.5*(2-omega).^2).*(omega>1 & omega<=2);
end
