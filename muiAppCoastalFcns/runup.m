function [R2,setup,swash,gR2] = runup(bs,H,T,depi,wv)
%
%-------function help------------------------------------------------------
% NAME
%   runup.m
% PURPOSE
%   wave runup magnitude
% USAGE
%   [R2,setup,swash,gR2] = runup(beta,H,T);      %assumes input is offshore 
%                                                %Hs and Tp
%   [R2,setup,swash,gR2] = runup(beta,H,T,depi); %applies inverse shoaling to  
%                                                %nearshore input wave conditions H
%   [R2,setup,swash,gR2] = runup(beta,H,T,1,wv); %ditto & uses T10 instead of 
%                                                %Tp if included in wv
%   Note: 1 is dummy depth if wv includes nearshore water depth wv.depi, 
%   if not, include depths as depi to use nearshore wave input option
% INPUTS
%   beta- beach foreshore slope (1:bs)
%   H   - offshore significant wave height (m)
%   T   - offshore wave period (s) - paper is not clear whether this is peak or mean
%   depi - inshore water depth - empty if input is offshore wave data
%   wv - dstable of detailed wave parameters from wave model including T10
%        (inverse wave period Tm-1,0)
% RESULTS
%   R2 - 2% runup (m)
%   setup - setup component of runup (m)
%   swash - swash component of runup (m)
%   gR2 - 2% runup using the eqn for gravel beaches by Poate etal (m)
% NOTES
%   based on equation proposed by Stockdon etal, Coastal Eng. 2006.
%   gravel beach equation using Poate etal, Coastal Eng.2016
%
% Author: Ian Townend
% CoastalSEA (c)June 2018
%--------------------------------------------------------------------------
%
    g = 9.81;       %acceleration due to gravity
    beta = 1./bs;   %beach slope as tan(alpha) or 1/bs
    
    if nargin<4, depi = []; end
    
    if nargin==5                               %additional wave parameters
        vars = wv.VariableNames;
        if any(strcmp(vars,'T10')), T = wv.T10; end
        Cp = 0.39;                             %Poate coefficient for T10
        %depths supplied in wv struct rather than depi
        if ~isempty(depi) && any(strcmp(vars,'depi')), depi = wv.depi; end        
    else
        Cp = 0.33;                             %Poate coefficient for Tp
    end
    
    if ~isempty(depi)
        [~,Ks] = shoaling(1,T,100,depi);
        H = H./Ks;                              %apply inverse shoaling
    end    
    
    %Compute setup and swash using Stockdon etal
    HoLo = H.*(g*T.^2/2/pi());
    setup = 0.35*beta.*sqrt(HoLo);
    swash = sqrt(HoLo.*(0.56*beta.^2+0.004))/2;
    R2 = 1.1*(setup+swash);                     %R2 Stockdon et al
    
    %Compute gravel runup using Poate etal
    gR2 = Cp*sqrt(beta).*T.*H;                  %R2 Poate et al (gravel+Tp or T10)
end