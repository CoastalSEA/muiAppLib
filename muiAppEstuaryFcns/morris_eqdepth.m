function deq = morris_eqdepth(sm,cn,qm,dslr)
%
%-------function help------------------------------------------------------
% NAME
%   morris_eqdepth.m
% PURPOSE
%   Solve the Morris equation for equilibrium depth
% USAGE
%   deq = morris_eqdepth(sm,cn,qm,dslr)
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
%   qm - sediment loading (s^-1)
%   dslr - rate of sea level rise (m/s)
% OUTPUTS
%   deq - equilibrium depth of marsh (returns real root or empty if no root found)
% NOTES
%   As qm->0, solution only found if kbm is sufficient to keep pace with
%   dslr. similarly if kbm->0 then qm must be large enough to offset dslr.
%   Hence. no solution implies that the marsh has drowned.
% SEE ALSO
%   see Saltmarsh class in Asmita
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
%
    kbm = sm.SpeciesProduct;  
    kbm = kbm/cn.y2s;           %rate of biomass production (m^2/kg/s)
    %kbm = kbm/acb;             %correct for bed density
    Bc = morris_biocoeffs(sm);  %biomass coefficients
    %
    if isrow(kbm), kbm = kbm'; end  %force a column vector
    KBM = repmat(kbm,1,3);      %replicate kbm for each coefficient
    kbc = Bc.*KBM;              %multiply each coefficient by kbm
    cfs = sum(kbc,1);           %sum each coefficient for all species
    cfs(3)=cfs(3)+qm;           %add in the sediment loading
    coefs = [cfs -dslr];        %construct polynomial
    coefs = coefs(:,[2 1 3 4]); %swap coefs a and b to give
                                %f(D)=kbD^3+kaD^2+(qm+kc)D-slr 
                                %where kb=kbm*b summed for all species, etc.                              
    deq = roots(coefs);         %find roots
    deq = deq(deq==real(deq));  %find real values including when imaginary part is zero
    deq = min(deq(deq>0));      %select minimum root>0 
    %alternative solution and check plot. For this problem, when roots 
    %returns imaginary solution the imaginary part is usually small and 
    %the real part a valid root. 
    % polyfunc = @(x) coefs(1)*x.^3+coefs(2)*x.^2+coefs(3)*x+coefs(4);
    % deqch = fzero(polyfunc, mean(sm.MaxSpDepth)/2);
    % checkPlot(sm,coefs,deqch)
end
%
function checkPlot(sm,coefs,deq) %#ok<DEFNU> 
    depths = min(sm.MinSpDepth)-0.1:0.001:max(sm.MaxSpDepth)+0.1;
    polyfunc = @(x) coefs(1)*x.^3+coefs(2)*x.^2+coefs(3)*x+coefs(4);
    hf = figure('Tag','PlotFig');  
    ax =axes(hf);    
    plot(ax,depths,polyfunc(depths))
    hold on
    plot(ax,deq,polyfunc(deq),'xr')
    hold off
    ax.XAxisLocation = 'origin';
end
