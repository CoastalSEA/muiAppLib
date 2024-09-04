function mp = profileslope(dep,swl,z1km,ubs)
%
%-------function help------------------------------------------------------
% NAME
%   profileslope.m
% PURPOSE
%   Calculate the bed slope at some depth within the surf zone
% USAGE
%   mp = profileslope(dep,swl,z1km,ubs)
% INPUTS
%   dep - depth of water at point (m)
%   swl - still water level (mOD)
%   z1km - bed level 1km from shore (mOD) [OR y,z coordinates]
%   ubs  - upper beach slope (1:bs)
% OUTPUTS
%   mp - slope of profile at 'dep'
% NOTES
%   Uses a Dean profile and invokes the coefficient 'A' from the depth
%   at a given distance offshore. 
%
% Author: Ian Townend
% CoastalSEA (c)April 2020
%----------------------------------------------------------------------
%
    if length(z1km)>1
        y1km = z1km(1);       %y and z have been specified
        z1km = z1km(2);
    else
        y1km = 1000;          %default distance
    end
    
    A = (0-z1km)/y1km^(2/3);
    z = swl-dep;
    mp = (-z).^0.5./(2/3*A^1.5);
    mp(z>=0) = ubs;    %restore correct upper beach slope
    mp(mp<ubs) = ubs;  %set ubs as the maximum slope (ie cannot be steeper)
end