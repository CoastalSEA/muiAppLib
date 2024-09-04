function [Hsb1,Hsb2] = hs_break(Hs,Tp,bs,dep,beta,g,hbflag)
%
%-------function help------------------------------------------------------% NAME
%   hs_break.m
% PURPOSE
%   Significant wave height after breaking for given water depth d
% USAGE
%   [Hsb1,Hsb2] = hs_break(swl,Hs,Tp,bs,bl,beta,g,flag)
% INPUTS
%   Hs   - significant wave height before breaking (m)
%   Tp   - peak period (s)
%   bs   - bed slope (1:bs)
%   dep  - bed level(m above datum)
%   beta - breaking coefficient (default=1)
%   g    - acceleration due to gravity(m/s2)
%   hbflag - choice of algorithm:, 0,1,2,3 used in Hb_break
%           0: simple ratio of 0.78 (maximum for a solitary wave)
%           1: SPM breaking on a slope
%           2: SPM breaking on a slope some distance in front of toe
%           3: SPM modified for overtopping some distance in front of toe
% RESULTS
%   Hsb1  - a lower bound estimate of the breaking value of Hs
%   Hsb2  - an upper bound estimate of the breaking value of Hs
% NOTES
%   Uses equations proposed by Holmes for a lower bound estimate and Tucker,
%   Carr and Pitt for an upper bound estimate
%   For energy equivalence should use Tp and Hrms (see Soulsby 1997, p69)
%   so use Tp to get Hb. The output values of Hsb based on Hb can then be
%   adjusted to an RMS value if requried (Hs/sqrt(2))
% SEE ALSO
%   Function hb_break.m which is used by this function
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%--------------------------------------------------------------------------
%

% find the breaker height for a single wave
Hb = hb_break(dep,bs,Tp,g,hbflag);

%distinguish between 0 and NaN values of Hs
idx = Hs==0;
% now find the upper and lower bound estimates of Hs
HbHs2 = 2*(Hb./Hs).^2;
fnHbHs = exp(-HbHs2);
% redistribution over the whole of the pdf below Hb: Holmes - a lower bound
Hsb1 = sqrt(beta^2*Hs.^2.*(1 - fnHbHs.*(1+HbHs2))./(1-fnHbHs));
Hsb1(idx) = 0;
% redistribution at Hb: Tucker, Carr & Pitt - an upper bound
Hsb2 = sqrt(beta^2*Hs.^2.*(1-fnHbHs));
Hsb2(idx) = 0;