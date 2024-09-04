function [alp,quad] = getalp(phi,theta)
%
%-------function help------------------------------------------------------
% NAME
%   getalp.m
% PURPOSE
%   Find the angle between the wave crest and the bed contour
% USAGE
%   [alp,quad] = getalp(phi,theta)
% INPUT
%   phi   - wave direction (degTN)
%   theta - angle of the bed contour (degTN) - looking seaward
% OUTPUT
%   alp  - angle between wave crest and contour (rads)
%   quad - quadrant of the wave direction relative to the contour
% NOTES
%   The contour is measured from TN to the contour looking seaward
%   quads 1 and 2 are onshore, quads 3 and 4 are offshore
%   quads 1 and 4 are waves from right to left when looking from seaward
%   quads 2 and 3 are waves from left to right
% SEE ALSO
%   used in refraction.m, beachtransportratio.m, littoraldrift.m, 
%   xshore_bailard.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%----------------------------------------------------------------------
% 
    beta = mod(phi-theta,360);
    limit = 0:90:360; %quadrants relative to contour
    quad = zeros(size(phi)); alp = quad;
    for i=1:length(phi)
        if beta(i)>=limit(1) && beta(i)<limit(2)
            quad(i) = 1;     %from theta to theta+90
            alp(i) = 90-beta(i);
        elseif beta(i)>=limit(2) && beta(i)<limit(3)
            quad(i) = 2;     %from theta+90 to theta+180
            alp(i) = beta(i)-90;
        elseif beta(i)>=limit(3) && beta(i)<limit(4)
            quad(i) = 3;     %from theta+180 to theta+270
            alp(i) = 270-beta(i);
        elseif beta(i)>=limit(4) && beta(i)<=limit(5)
            quad(i) = 4;     %from theta+270 to theta
            alp(i) = beta(i)-270;
        else
            alp(i) = NaN; quad(i) = NaN;
        end
    end
    alp = alp*pi()/180;
end