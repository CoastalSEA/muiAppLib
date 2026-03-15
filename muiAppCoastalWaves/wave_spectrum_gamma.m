function gamma = wave_spectrum_gamma(SG,f,dir)
%
%-------function help------------------------------------------------------
% NAME
%   wave_spectrum_gamma
% PURPOSE
%   Estiamte the JONSWAP gamma from a frequency spectrum, S(f), or a
%   frequency–direction spectrum S(dir,f)
% USAGE
%   gamma = wave_spectrum_gamma(SG,f,dir);
% INPUTS
%   SG - Spectrum, 1D frequency, SG(f) or 2D frequency–direction spectrum
%        SG(dir,f)
%   f - frequency vector (Hz)
%   dir - direction vector (degrees) - only required if 2D spectrum used
% OUTPUT
%    gamma - best-fit JONSWAP peak enhancement factor
% SEE ALSO
%   
%   
% CoastalSEA (c) March 2026
%--------------------------------------------------------------------------
%  
    lb = 0.1; ub = 7.0;    %bounds of search
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
    obj = @(gamma) jonswap_misfit(gamma, f, Sf, fp);

    % Search gamma in a realistic range
    gamma = fminbnd(obj, lb, ub, optimset('TolX',1e-6,'Display','off'));
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
