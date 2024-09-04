function c=celerity(T,d)
%
%-------function help------------------------------------------------------
% NAME
%   celerity.m 
% PURPOSE
%   Calculate the wave celerity using Hunt's equation
% USAGE
%   c=celerity(T,d)
% INPUTS
%   T - wave period (s)
%   d - water depth (m
% RESULTS
%   c - wave celerity (m/s)
% NOTES
%   For method see Hunt,ASCE,WW4,1974,p457-459
%
% Author: Ian Townend
% CoastalSEA (c)June 2016
%--------------------------------------------------------------------------
%
    g = 9.81;
    %trap negative depths
    d(d<0) = NaN;
    %
    f = 1./T;
    y = (2*pi*f).^2.*d/g;
    A = (y+(1+0.6522*y+0.4622*y.^2+0.0864*y.^4+0.0675*y.^5).^-1).^-1;
    c = sqrt(g*d.*A);