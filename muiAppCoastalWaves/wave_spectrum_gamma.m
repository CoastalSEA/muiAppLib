function gamma = wave_spectrum_gamma(SG,f,bnd,dir,islog)
%
%-------function help------------------------------------------------------
% NAME
%   wave_spectrum_gamma.m
% PURPOSE
%   Estimate the JONSWAP gamma from a frequency spectrum, S(f), or a
%   frequency–direction spectrum S(dir,f)
% USAGE
%   gamma = wave_spectrum_gamma(SG,f,bnd,dir,islog);  %for 2D spectrum
%   gamma = wave_spectrum_gamma(SG,f,bnd,[],islog);   %for 1D spectrum
%   gamma = wave_spectrum_gamma(SG,f,bnd);  %for 1D spectrum linear option
% INPUTS
%   SG - Spectrum, 1D frequency, SG(f) or 2D frequency–direction spectrum
%        SG(dir,f)
%   f - frequency vector (Hz)
%   bnd - vector [1x2] for lower and upper bounds of gamma search
%   dir - direction vector (degrees) - only required if 2D spectrum used
%   islog - logical true to use log-space misfit, false for linear option
% OUTPUT
%    gamma - best-fit JONSWAP peak enhancement factor
% SEE ALSO
%   wave_spectrum.m. Called from ctWaveSpectraPlots.estimateSpectrumGamma
%   
% CoastalSEA (c) March 2026
%--------------------------------------------------------------------------
%  
    if nargin<5, islog = false; end

    %collapse spectrum to 1D if 2D input
    if ~isvector(SG)
        if size(SG,1) ~= numel(dir), SG = SG.'; end
        theta = deg2rad(mod(dir(:),360));  % [I x 1]
        wp = trapz_weights_periodic(theta);
        Sf = sum(SG.*wp);
    else
        Sf = SG(:);
    end

    % --- 2. Identify peak frequency ---
    [~, idxp] = max(Sf);
    fp = f(idxp);
    %Tp = 1/fp;

    % --- 3. Fit gamma by minimising spectral misfit ---
    if islog
        obj = @(gamma) jonswap_misfit_log(gamma, f, Sf, fp);
    else
        obj = @(gamma) jonswap_misfit(gamma, f, Sf, fp);
    end

    % Search gamma in a realistic range
    gamma = fminbnd(obj,bnd(1),bnd(2),optimset('TolX',1e-6,'Display','off'));
end

%%
function J = jonswap_misfit(gamma, f, Sf_obs, fp)
    % Standard sigma values
    sigma1 = 0.07; 
    sigma2 = 0.09;

    % Build JONSWAP shape with unit alpha
    sigma = sigma1*(f <= fp) + sigma2*(f > fp);
    r = f./fp;

    Sshape = r.^(-5) .* exp(-1.25 * r.^(-4)) .* ...
             gamma.^(exp(-((r - 1).^2) ./ (2*sigma.^2)));

    % Scale alpha so that total energy matches observed m0
    m0_obs = trapz(f, Sf_obs);
    m0_shape = trapz(f, Sshape);
    alpha = m0_obs / m0_shape;

    Sf_mod = alpha * Sshape;

    % Misfit: L2 norm in spectral space
    J = trapz(f, (Sf_mod - Sf_obs).^2);
end

%%
function J = jonswap_misfit_log(gamma, f, Sf_obs, fp)
    % Standard sigma values
    sigma1 = 0.07;
    sigma2 = 0.09;

    % Build JONSWAP shape with unit alpha
    sigma = sigma1*(f <= fp) + sigma2*(f > fp);
    r = f./fp;

    Sshape = r.^(-5) .* exp(-1.25 * r.^(-4)) .* ...
             gamma.^(exp(-((r - 1).^2) ./ (2*sigma.^2)));

    % Scale alpha so that total energy matches observed m0
    m0_obs   = trapz(f, Sf_obs);
    m0_shape = trapz(f, Sshape);
    alpha    = m0_obs / m0_shape;

    Sf_mod = alpha * Sshape;

    % --- Log-space misfit ---
    % Avoid log(0): clip very small values
    eps_val  = 1e-12;
    Sf_obs_c = max(Sf_obs, eps_val);
    Sf_mod_c = max(Sf_mod, eps_val);

    log_diff = log(Sf_mod_c) - log(Sf_obs_c);

    % Optionally restrict to a band around the peak (often helpful)
    mask = (f > 0.5*fp) & (f < 2.5*fp);
    log_diff = log_diff(mask);
    f_int    = f(mask);
    J = trapz(f_int, log_diff.^2);

    %J = trapz(f, log_diff.^2);
end