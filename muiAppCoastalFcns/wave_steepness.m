function [t,S] = wave_steepness(H,T,dep,t)
%
%-------function help------------------------------------------------------
% NAME
%   wave_steepness.m
% PURPOSE
%   compute the wave steepness for given wave height, wave period and water
%   depth
% USAGE
%   S = wave_steepness(Hs,Tp,dep)
% INPUTS: 
%   H - wave height (m)
%   T - wave period (s)
%   dep - water depth (m)
%   t - input time (datetime)
% OUTPUTS 
%   t - input time is returned (e.g. for use in muiUserModel)
%   S - wave steepness based on linear wave theory with depth dependent
%   wave length.
% REQUIREMENTS
%   uses celerity.m
%
% Author: Ian Townend
% CoastalSEA (c) March 2023
%----------------------------------------------------------------------
%
    if length(H)~=length(T)
        warndlg('Wave height and period are different lengths')
        S = []; return;
    end

    c = celerity(T,dep);
    L = c.*T;
    S = H./L;
end