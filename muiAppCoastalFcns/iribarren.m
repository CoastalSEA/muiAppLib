function [Iri,Itype] = iribarren(Hs,Tp,bs,g)
%
%-------function help------------------------------------------------------
% NAME
%   iribarren.m
% PURPOSE
%   function to calculate the Iribarren number which characterise the 
%   wave breaker type 
% USAGE
%   [Ira,Itype] - iribarren(Hs,Tp,bs,g)
% INPUTS
%   Hs  - incident significant wave height (m)
%   Tp  - peak wave period (s)
%   bs  - beach slope(1:bs)
%   g   - acceleration due to gravity (m/s2)
% RESULTS
%   Ira - Iribarren number as numeric array
%   Itype - Type of breaker as cell array
%          - 3=surging, 2=plunging, 1=spilling  (order reversed Apr 2020)
% NOTES
%   Irabarren number: tan(alpha)/sqrt(Hb/Lo)
% SEE ALSO
%   Dimensionless fall velocity as used in beach_type.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
Len = g.*Tp.^2/2/pi();               %deep water wave length

Steep = Hs./Len;                     %wave steepness
check = Steep>0.14 | Steep <0.0001;  %check that wave is not too steep (1/7)
Steep(check) = NaN;                  %or too flat (eg Hs<<0.1)
Iri = (1./bs)./sqrt(Steep);          %Iribarren number
Iri(Iri==Inf) = NaN;                 %remove Infinite values
nrec = length(Iri);
Itype = zeros(nrec,1);   err = 1e-2;
%breaker types based on breaking wave heights nearshore. Deepwater values
%are different <0.5 =spilling; >0.5 plunging <3.3; >3.3 surging 
for ij=1:nrec
    if Iri(ij)>2
        Itype(ij) = 3-err;   %surging
    elseif Iri(ij)<=2 && Iri(ij)>0.4
        Itype(ij) = 2-err;   %plunging
    else
        Itype(ij) = 1-err;   %spilling
    end
end


