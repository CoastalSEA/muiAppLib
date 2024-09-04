function [eqdepth,sedload] = interpdepthload(sm,cn,aws,c0,dslr)
%
%-------function help------------------------------------------------------
% NAME
%   interpdepthload.m
% PURPOSE
%  Iterate to find the equilibrium depth and sediment load,qm(s^-1)
% USAGE
%   [depth,sedload] = interpdepthload(sm,cn,aws,qm0,dslr)
% INPUTS
%   sm - Saltmarsh instance or struct of Saltmarsh properties, with:
%        NumSpecies - number of species
%        MinSpDepth - minimum depth for each species (m)
%        MaxSpDepth - maximum depth for each species (m)
%        MaxBiomass - maximum biomass for each species (kg/m2)
%        SpeciesProduct - species productivity (m2/kg/yr)
%        SettlingAlpha - coefficient for biomass dependent enhanced settling rate (m/s)
%        SettlingBeta - exponent for biomass dependent enhanced settling rate  
%   cn - struct of abbreviated Constants values (eg cn.g for Gravity)
%   aws(1) = vertical exchange of element with no biology (m/s)
%   aws(2) = value for depths greater than dmx (optional, default=aws(1))  
%   c0 - background concentration adjacent to marsh (-)   
%   dslr - rate of sea level rise (m/s)
% OUTPUTS
%   eqdeqth - equilibrium depth (m)
%   sedload - sediment load, qm (s^-1)
% SEE ALSO
%   see Saltmarsh class in Asmita
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
%

    %concentration over marsh as a function of depth   
    cdep = sm.MarshDepthConc.Depth;
    conc = sm.MarshDepthConc.Concentration;
    mdep(1) = cdep(find(conc>0,1,'first'));  %minimum depth in concentration array
    mdep(2) = max(cdep);                     %maximum depth in concentration array
    fundep = @(x) abs(getNewDepth(sm,cn,aws,cdep,conc,mdep,dslr,c0,x)-x);
    % options = optimset('FunValCheck','on','TolX',1e-6,'Display','iter','PlotFcns',@optimplotfval);
    options = optimset('FunValCheck','on','TolX',1e-4);
    thedepth = fminbnd(fundep,mdep(1),mdep(2),options);
    [eqdepth,sedload] = getNewDepth(sm,cn,aws,cdep,conc,mdep,dslr,c0,thedepth);
end
%%
function [newdepth,sedload] = getNewDepth(sm,cn,aws,cdep,conc,mdep,dslr,c0,depth)
    %find the marsh equilbrium depth and sediment load for given initial
    %depth concentration over the marsh and rate of sea level rise.
    if depth>=mdep(1) && depth<=mdep(2) 
        cem = interp1q(cdep,conc,depth);
        wsm = bioenhancedsettling(sm,depth,aws);
        sedload = wsm*cem/depth;
    elseif depth<mdep(1)
        sedload = 0;
    else %depth>maxdepth
        sedload = aws(2)*c0/mdep(2);
    end
    newdepth = morris_eqdepth(sm,cn,sedload,dslr);
    if isempty(newdepth), newdepth = 0; end
end

