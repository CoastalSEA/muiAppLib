function gamma = estimate_jonswap_gamma(inp,f,gamma_0,istma)
%
%-------function help------------------------------------------------------
% NAME
%   estimate_jonswap_gamma.m
% PURPOSE
%   estimate the value of the JONSWAP peakiness parameter from the wave 
%   parameters provided by model hindcasts such as the Copernicus reanalysis
% USAGE
%   gamma = estimate_jonswap_gamma(inp,f,istma);
% INPUTS
%   inp -a struct that includes:
%       Hs - significant wave height [4sqrt(mo)] (m) <scalar or vector>
%       Tp - peak wave period (s)  <can be vector - same length as Hs>
%       T1 - mean wave period, m0/m1 (s) 
%       T2 - second moment wave period, sqrt(mo/m2) (s)
%       T10- inverse wave wave period, m-1/m0 (s)
%   f - spectrum frequency vector
%   gamma_0 - if gamma_0>0, return gamma = gamma_0, otherwise
%             if gamma_0<=0 use this as the default estimate of gamma
%             which is used if solution not found.
%   istma - logical true if form being used is TMA rather than Jonswap
% OUTPUT
%   gamma - JONSWAP peakiness parameter
% NOTES
%   This code provides a number of options for estimating gamma but is hard
%   coded to use use spectral moments and Hs, Tp, T1, T2, T-10 to define an
%   objective function and obtain the minimum within bounds. Other options
%   were developed to explore different ways of solving the problem and
%   some perform better than others depending on the form of the spectrum.
% SEE ALSO
%   used in ctWaveSpectrum to provide gamma value for use in wave_spectrum
%   
% Author: Ian Townend
% CoastalSEA (c) March 2026
%--------------------------------------------------------------------------
% 
    if gamma_0>0, gamma = gamma_0; return; end

    bnds = [0,7];         %hard coded**    
    lambda = 0.1;         %prior weighting - hard coded******
    gamma_option = 1;     %options to estimate gamma - hard coded******
    sigmas = [0.07,0.09]; %default values of sigma
    gamma0 = 3.3;         %default gamma (overwritten by negative input value)
    if gamma_0<0, gamma0 = -gamma_0; end
    
    fp = 1/inp.Tp;
    %estimate the best fit for gamma. Does not always fall
    %within realistic region and gamma0 is used  
    switch gamma_option
        case 1
            %use spectral moments and Hs, Tp, T1, T2, T-10
            gamma = fit_gamma_multimoment(inp,f,bnds,sigmas,gamma0,lambda);
        case 2
            %estimate gamma using Hs, T1 and steepness
            [gamma,~] = gamma_estimation(inp,f,fp,sigmas,bnds,istma);
        case 3
            %estimate gamma using Hs, Tp and the sea state component value of T1
            [gamma, ~] = jonswap_gamma_estimation(inp,bnds);
        case 4
            %estimate gamma using the mean of the regressions for T1, T2, and T10
            [gamma, ~] = carter_gamma_from_periods(inp,bnds);
            if isnan(gamma), gamma = gamma0; end
        case 5
            %estimate gamma using the integral of the spectrum either side of T1
            gamma = gamma_integral_estimation(inp,f,fp,sigmas,bnds,istma);
        case 6
            %simply use the relationship for T1
            I0I1 = inp.T1/inp.Tp;                   %Carter eq(23)
            if I0I1>0.7 && I0I1<0.9
                gamma = 45.3*(I0I1)^14.6;           %fit to MIAS No.4, Table 1
            else
                gamma = gamma0;
            end                            
    end

    %use gamma0 if estimate is on a boundary. This can make the
    %overall fit of a timeseries poorer but it removes the
    %dependency on the choice of bounds.
    if gamma<=bnds(1)+0.01 || gamma>=bnds(2)-0.01,  gamma = gamma0; end  %use default
end

%%
function gamma_est = fit_gamma_multimoment(inp,f,bnds,sigmas,g0,lambda)%(Hs, Tp, T1_obs, T2_obs, T10_obs)
    %estimate gamma using spectral moments and Hs, Tp, T1, T2, T-10

    % Inputs
    if any(ismatch(inp.Properties.VariableNames,'T1'))  %COP reanalysis includes all parameters
        Hs = inp.Hs; Tp = inp.Tp; T1_obs = inp.T1; T2_obs = inp.T2; T10_obs = inp.T10;
    else                 %Basic wave buoy data only has Tp and Tz (==T2).
        Hs = inp.Hs; Tp = inp.Tp; T1_obs = 0; T2_obs = inp.Tz; T10_obs = 0;
    end
    % flim = [0.025,0.58]; %frequency limits
    % f = [flim(1):0.005:0.1,0.11:0.01:flim(2)];  %observed frequency intervals (spt format)
    % Weights and prior
    w1  = 1; w2 = 1; w10 = 1;
    usePrior = true;

    % Objective in gamma
    obj = @(g) objective_gamma(g, Hs, Tp, f, ...
                               T1_obs, T2_obs, T10_obs, ...
                               w1, w2, w10, usePrior, sigmas, g0, lambda);
    % Bounded search
    gamma_est = fminbnd(obj, bnds(1),bnds(2), optimset('TolX',1e-6,'Display','off'));
end

%%
function J = objective_gamma(gamma, Hs, Tp, f, ...
                             T1_obs, T2_obs, T10_obs, ...
                             w1, w2, w10, usePrior, sigmas, gamma0, lambda)
    %objective function for fit_gamma_multimoment
    fp = 1./Tp;
    sigma = sigmas(1)*(f<=fp) + sigmas(2)*(f>fp);

    % JONSWAP shape (alpha chosen to match Hs via m0)
    % Start with unit alpha, then rescale
    Sshape = jonswap_shape_unit(f, fp, gamma, sigma);
    m0_unit = trapz(f, Sshape);
    alpha = (Hs/4)^2 / m0_unit;   % enforce m0 -> Hs

    S = alpha * Sshape;

    % Moments
    m0  = trapz(f, S);
    m1  = trapz(f, f    .* S);
    m2  = trapz(f, f.^2 .* S);
    m_1 = trapz(f, f.^(-1) .* S);

    T1_mod  = m0 / m1;
    T2_mod  = sqrt(m0 / m2);
    T10_mod = m_1 / m0;

    e1 = 0; e2 = 0; e10 = 0;
    if T1_obs>0, e1  = (T1_mod  - T1_obs ) / T1_obs; end    
    if T2_obs>0, e2  = (T2_mod  - T2_obs ) / T2_obs; end
    if T10_obs>0, e10 = (T10_mod - T10_obs) / T10_obs; end

    J = w1*e1.^2 + w2*e2.^2 + w10*e10.^2;

    if usePrior
        J = J + lambda * ((gamma - gamma0)/gamma0).^2;
    end
end

%%
function Sshape = jonswap_shape_unit(f, fp, gamma, sigma)
    % Unit-alpha JONSWAP shape (m^2/Hz per unit alpha)
    r  = (f./fp);
    Sshape = r.^(-5) .* exp(-1.25 * r.^(-4)) .* ...
             gamma.^(exp(- ( (r - 1).^2 ) ./ (2*sigma.^2)));
end

%%
function gamma = gamma_estimation(inp,f,fp,sigmas,bnds,istma)
    %estimate gamma using Hs, T1 and the wave steepness
    % Inputs
    Hs = inp.Hs; T1 = inp.T1; ds = inp.ds; 
    lb = bnds(1)-1e-3; ub = bnds(2)+1e-3; 

    % fun_I0 = @(gamma,sigma) spectral_moment(0,gamma,sigma);                         %Carter eq(18)
    % fun_alp = @(g,s) (2*pi).^4*Hs.^2.*fp.^4./(16*g.^2.*fun_I0(g,s)); %Carter eq(20)  
    obj1 = @(gamma) abs(computeT1(fp,gamma,sigmas) - T1);
    obj2 = @(gamma) abs(computeHs(f,fp,gamma,sigmas,Hs,ds,istma) - Hs);
    obj = @(gamma) abs(obj2(gamma)/obj1(gamma)^2 -Hs/T1^2);      %steepness
    % search gamma in [lb,ub]    
    gamma = fminbnd(obj, lb, ub, optimset('TolX',1e-6,'Display','off'));
end

%%
function T1 = computeT1(fp,gamma,sigmas)
    %estimate T1 for given spectrum parameters 
    I0 = spectral_moment(0,gamma,sigmas);                %Carter eq(18)
    I1 = spectral_moment(1,gamma,sigmas);
    T1 = I0/I1/fp;
end

%%
function Hs = computeHs(f,fp,gamma,sigmas,Hs,ds,istma)
    %estimate T1 for given spectrum parameters 
    g= 9.81;
    % const = g^2*(2*pi)^-4;
    I0 = spectral_moment(0,gamma,sigmas);           %Carter eq(18)
    alpha = (2*pi).^4*Hs.^2.*fp.^4./(16*g.^2.*I0); %Carter eq(20) 
    % Hs = sqrt(16*alpha*const*(fp^-4)*I0);
    S = getSpectrum(f,fp,alpha,gamma,sigmas,ds,istma);
    Hs = sqrt(16*trapz(f,S));
end

%%
function [gamma, gamma_sensitivity] = jonswap_gamma_estimation(inp,bnds,sigmas)
    %estimate gamma using Hs, Tp and the sea state component value of T1
    Hs = inp.Hs; Tp = inp.Tp; T1 = inp.T1;
    flim = [0.025,0.58]; %frequency limits
    f = [flim(1):0.005:0.1,0.11:0.01:flim(2)];  %obbserved frequency intervals (spt format)

    %f = linspace(1/(5*Tp), 2/Tp, 2000);   % frequency grid (adjust as needed)
    ncomp = numel(Hs);
    gamma = nan(size(Hs));
    gamma_sensitivity = nan(size(Hs,2),3); % store gamma for alt sigma sets

    for i=1:ncomp
        Hs = Hs(i);
        T1_obs = T1(i);

        % objective: for given gamma compute model T1 and return error
        obj = @(gamma) abs( compute_T1_from_jonswap(Hs, Tp, gamma, f, sigmas) - T1_obs );

        % search gamma in [lb,ub]
        lb = bnds(1)-1e-3; ub = bnds(2)+1e-3; 
        gamma_est = fminbnd(obj,lb,ub,optimset('TolX',1e-6,'Display','off'));
        gamma(i) = gamma_est;

        % sensitivity: try three sigma variants
        gamma_sensitivity(i,1) = fminbnd(@(g) abs(compute_T1_from_jonswap(Hs,Tp,g,f,0.06,0.08)-T1_obs),lb,ub);
        gamma_sensitivity(i,2) = fminbnd(@(g) abs(compute_T1_from_jonswap(Hs,Tp,g,f,0.07,0.09)-T1_obs),lb,ub);
        gamma_sensitivity(i,3) = fminbnd(@(g) abs(compute_T1_from_jonswap(Hs,Tp,g,f,0.08,0.10)-T1_obs),lb,ub);
    end
end

%%
function T1 = compute_T1_from_jonswap(Hs, Tp, gamma, f, sigmas)
    fp = 1/Tp;
    sigma = sigmas(1)*(f<=fp) + sigmas(2)*(f>fp);
    % JONSWAP shape (Goda / common parameterisation)
    A = 0.0624 ./ (0.230 + 0.0336*gamma - 0.185./(1.9+gamma)); % empirical alpha prefactor
    S = A .* Hs.^2 .* Tp.^(-4) .* f.^(-5) .* exp(-1.25*(Tp.*f).^(-4)) .* ...
        gamma.^(exp(-((Tp.*f - 1).^2)./(2*sigma.^2)));

    m0 = trapz(f, S);
    m1 = trapz(f, f .* S);
    T1 = m0 ./ m1;
end

%%
function [gamma_hat, gamma_parts] = carter_gamma_from_periods(inp,bnds)
    %estimate gamma using the mean of the regressions for T1, T2, and T10
    R1  = inp.T1 ./inp.Tp;
    R2  = inp.T2 ./inp.Tp;
    R10 = inp.T10./inp.Tp;

    % Carter-style regressions (from your excerpt)
    gamma1  = 45.3 .* (R1 ).^14.6;
    gamma2  = 69.7 .* (R2 ).^12.2;
    gamma10 = 34.4 .* (R10).^23.0;

    gamma_parts = [gamma1(:), gamma2(:), gamma10(:)];

    % Basic validity mask (1–7 range, say)
    valid = gamma_parts >= bnds(1) & gamma_parts <= bnds(2);

    % Weighted mean: equal weights for valid entries
    w = double(valid);
    w(w==0) = NaN;                 % ignore invalid
    gamma_hat = mean(gamma_parts.*w, 2,'omitnan');
end

%%
function gamma = gamma_integral_estimation(inp,fp,sigmas,bnds,istma)
    %estimate gamma using the integral of the spectrum either side of T1
    obj = @(gamma) compute_objective_function(inp,fp,gamma,sigmas,istma);

    gamma = fminbnd(obj,bnds(1),bnds(2),optimset('TolX',1e-6,'Display','off'));
end

%%
function obj = compute_objective_function(inp,fp,gamma,sigmas,istma)
    %objective function for gamma_integral_estimation
    f1 = 1/inp.T1;
    ds = inp.ds;
    Hs = inp.Hs;
    g = 9.81;
    I0 = @(gam) spectral_moment(0,gam,sigmas);                  %Carter eq(18)
    alpha = @(gam) (2*pi).^4*Hs.^2.*fp.^4./(16*g.^2.*I0(gam)); %Carter eq(20) 
    Func = @(f) getSpectrum(f,fp,alpha(gamma),gamma,sigmas,ds,istma);
    In1 = integral(Func,0,f1,'RelTol',1e-3,'AbsTol',1e-3);
    In2 = integral(Func,f1,Inf,'RelTol',1e-3,'AbsTol',1e-3);
    obj = abs(In1-In2);
end

%%
function In = spectral_moment(nm,gamma,sigmas)
    %find the nm spectral moment of the Jonswap spectrum for given gamma
    sigma = @(f) sigmas(1).*(f<=1) + sigmas(2).*(f>1);   
    q = @(f) exp((-(f-1).^2)./(2.*sigma(f).^2));
    Func = @(f) (f.^nm).*(f.^-5).*exp(-1.25.*f.^-4).*gamma.^q(f);
    In = integral(Func,0,Inf,'RelTol',1e-3,'AbsTol',1e-3);
end

%%
function S = getSpectrum(f,fp,alpha,gamma,sigmas,ds,istma)
    %construct the Jonswap or TMA spectrum
    g = 9.81;
    sigma = @(f) sigmas(1).*(f<=fp) + sigmas(2).*(f>fp);   %Hughes eq(5 & 25), or
    q = @(f) exp((-(f-fp).^2)./(2.*sigma(f).^2.*fp.^2)); %Carter eq(16)
    const = g^2*(2*pi)^-4;
    Jonswap = @(f) const*alpha.*(f.^-5).*(exp(-1.25*(f./fp).^-4)).*(gamma.^q(f));

    S = zeros(length(fp),length(f));
    if istma                                          %TMA spectrum
        if isempty(ds) || ds<=0
            ds = (g./fp./2./pi)./fp/2;                %use deep water 
        end 
        Func = @(f) kit_limit(f,ds).*Jonswap(f);      %Hughes eq(15)
        for i=1:length(f)
            S(:,i) = Func(f(i));                          
        end

    else                                              %Jonswap spectrum
        for i=1:length(f)
            S(:,i) = Jonswap(f(i));
        end
    end
end

%%
function Lp = wave_length(inp)
    %find peak period based on wind or wave input
    g = 9.81;
    Tp = inp.Tp;
    %get the wave length at peak frequency
    Lp = (g*Tp./2./pi).*Tp;                   %use deep water celerity
    if isfield(inp,'df') && ~isempty(inp.df)
        if isscalar(inp.df) && inp.df>0       %scalar depth over fetch
            Lp = celerity(Tp,inp.df).*Tp;     %use Hunt eq.for celerity
        elseif length(inp.df)==length(Tp)     %vector depths over fetch
            for i=1:length(inp.df)
                if inp.df(i)>0
                    Lp(i) = celerity(Tp(i),inp.df(i)).*Tp(i);  %use Hunt eq.for celerity
                end
            end
        elseif length(inp.df)~=length(Tp)
            errordlg('Input vectors to wave_spectrum must be the same length')
        else
            %defaults to deepwater for all other combinations
        end
    end 
end