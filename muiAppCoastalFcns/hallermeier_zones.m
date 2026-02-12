function [hc,hb,smp] = hallermeier_zones(H,T,stdH,d50,rhow,rhos)
%-------header-------------------------------------------------------------
% NAME
%   bmvprofile.m
% PURPOSE
%   To compute the limits of the surf and shoal zones and returning the
%   surf zone limit, hb and the closure depth or shoal limit, hc.
% USAGE
%   [hc,hb] = hallermeier_zones(H,T,stdH,d50,rhow,rhos)
% INPUTS (except for plot)
%   Hs   - mean significant wave height (m)
%   Ts   - mean significant wave period (s)
%   stdH - standard deviation of wave height (m)
%   d50  - median sediment grain size (m)
%   rhow - density of sea water (kg/m3)
%   rhos - density of sediment (kg/m3)
% OUTPUT
%   hc - seaward closure depth or shoal zone limit
%   hb - landward surf zone limit
%   smp - struct with simplified estimates for
%         hc1 = (H - 0.3*stdH)*T/sqrt(g/5000/d50) Eq.7
%         hc2 = 0.018*H*T*sqrt(g/(d50*gamma))     Hallermeier, '83
%         hb = 2*H+11*stdH                        Eq.6
% NOTES
%   see Hallermeier, 1981, A profile zonation for seasonal sand beaches from 
%   wave climate. Coastal Engineering, 4 (3), 253-277.
%   "It is recommended that tidal effects be considered, as in 
%   Hallermeier (1978), by taking the calculated water depths as being with 
%   respect to mean low water level at a locality."
% SEE ALSO
%   test_hallermeier.m script for comparison results from Table II of paper
%   used in Sim_BMVmodel.m
%
% AUTHOR
%   Ian Townend
% COPYRIGHT
%    CoastalSEA (c) April 2020
%--------------------------------------------------------------------------
%
    g = 9.81;
    gamma = (rhos-rhow)/rhow;
    gTpi = g.*T.^2/(4*pi()^2);

    %seaward shoal zone limit - closure depth, Table 1 eqns for di
    %limit of significant on/offshore transport
    denom = 8*gamma*g*d50.*T.^2; 
    etac = asinh(sqrt(pi()^2*(H-0.3*stdH).^2./denom));
    hc = etac.*tanh(etac).*gTpi;

    % % Landward surf‑zone limit (Table 1 equations for dl)
    etab = NaN(size(H));    
    for i = 1:numel(H)

        % Skip invalid inputs early
        if isnan(H(i)) || isnan(stdH(i)) || isnan(gTpi(i)) || isnan(etac(i))
            continue
        end

        % Compute target value
        htprop = (H(i) + 5.6*stdH(i)).^2 ./ (0.03 * gamma * 4 * gTpi(i).^2);

        % Define function of etb only
        etafun = @(etb) etb .* sinh(etb).^2 .* tanh(etb) - htprop;

        % Bracket: use [0, etac(i)] instead of a single guess
        a = 0;
        b = max(etac(i), 1e-6);   % ensure positive upper bound

        % Evaluate endpoints safely
        fa = etafun(a);
        fb = etafun(b);

        if ~isfinite(fa) || ~isfinite(fb)
            continue
        end

        % Require a sign change for robustness
        if sign(fa) == sign(fb)
            continue
        end

        % Solve
        [root,~,exitflag] = fzero(etafun, [a b]);

        if exitflag > 0 && isfinite(root)
            etab(i) = root;
        end
    end

    hb = etab.*tanh(etab).*gTpi;  

    %simplified equations
    smp.hc1 = (H - 0.3*stdH).*T./sqrt(g/5000./d50); %Eq.7
    smp.hc2 = 0.018.*H.*T.*sqrt(g./(d50*gamma));    %Hallermeier, 1983
    smp.hb = 2*H+11.*stdH;                          %Eq.6