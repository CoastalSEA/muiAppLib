function [hc,hb] = hallermeier_zones(H,T,stdH,d50,rhow,rhos)
%-------header-------------------------------------------------------------
% NAME
%   bmvprofile.m
% PURPOSE
%   To compute the limits of the surf and shoal zones and returning the
%   surf zone limit, hb and the closure depth or shoal limit, hc.
% USAGE
%   [hc,hb] = hallermeierZones(H,T,stdH,d50,rhow,rhos)
% INPUTS (except for plot)
%   Hs   - mean wave height (m)
%   Tp   - mean peak wave period (s)
%   stdH - standard deviation of wave height (m)
%   d50  - median sediment grain size (m)
%   rhow - density of sea water (kg/m3)
%   rhos - density of sediment (kg/m3)
% OUTPUT
%   hc - seaward closure depth or shoal zone limit
%   hb - landward surf zone limit
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

    %seaward shoal zone limiti - closure depth 
    %limit of significant on/offshore transport
    denom = 8*gamma*g*d50.*T.^2;
    etac = asinh(sqrt(pi()^2*(H-0.3*stdH).^2./denom));
    hc = etac.*tanh(etac).*gTpi;

    %landward surf zone limit
    etab = zeros(size(H));
    for i=1:length(H)
        htprop = (H(i)+5.6*stdH(i)).^2./(0.03*gamma*4*gTpi(i).^2);
        etafun = @(etb) etb.*sinh(etb).^2.*tanh(etb)-htprop;
        etab(i) = fzero(etafun,etac(i)/4);
    end
    hb = etab.*tanh(etab).*gTpi;