function Hs1 = shoaling(Hs0,Tp,dep0,dep1)
%
%-------function help------------------------------------------------------% NAME
%   shoaling.m
% PURPOSE
%   Plane bed wave shoaling using linear wave theory
% USAGE
%   Hs1 = shoaling(Hs0,Tp,dep0,dep1)
% INPUTS
%   Hs0  - input significant wave height (m)
%   Tp   - peak wave period (s)
%   dep0 - input water depth (m)
%   dep1 - output water depth (m)
% RESULTS
%   Hs1 - output wave height
% NOTES
%   USACE Shore Protection Manual
% SEE ALSO
%   calls celerity.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2018
%--------------------------------------------------------------------------
%

% find wave celerity at input and output depths
c0 = celerity(Tp,dep0);
ci = celerity(Tp,dep1);

% group celerity and shoaling coefficient (linear wave theory)
pidL0 = 4*pi()*dep0./(c0.*Tp);
cg0 = c0/2.*(1+pidL0./sinh(pidL0));
pidLi = 4*pi()*dep1./(ci.*Tp);
cgi = ci/2.*(1+pidLi./sinh(pidLi));
Ks  = sqrt(cg0./cgi);
%output wave as multiple of input wave height and shoaling coefficient
Hs1 = Hs0.*Ks;