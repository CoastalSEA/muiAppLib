function [acem,asmt] = marsh_conc(c0,zm,wsm,wl,options)
%
%-------function help------------------------------------------------------
% NAME
%   marsh_conc.m
% PURPOSE
%   Calculate the concentration over a marsh flat using method proposed 
%   by Krone(1987).
% USAGE
%   [acem,asmt] = marsh_conc(c0,zm,wsm,wl,options)
% INPUTS
%   c0 - concentration in adjacent element (uses cE for flat elements)
%   zm - elevation of marsh (assumes symmetric tide)
%   wsm - settling rate over marsh 
%   wl - struct used in simple_tide (set unwanted constituents to zero)
%        MSL0                      mean tidel level to ordnance datum (mOD)
%        TidalAmp                  tidal elevation amplitude (m)
%        ElevPhase                 phase of elevation (ie k.x) (rads)
%        VelocityAmp               tidal velocity amplitude (m/s)
%        VelocityPhase             phase of velocity (ie k.x+phi) (rads)
%        M2amplitude               M2 tidal amplitufe (m)
%        S2amplitude               S2 tidal amplitufe (m)
%        O1amplitude               O1 tidal amplitufe (m)
%   options - struct with the following settings (default values shown
%        tsn  - duration of spring-neap cycle (days) [default = 14.77 days]
%        delt - time step (secs)  [default = 1 sec]
%        dmin - minimum depth used in calculations (m) [default = 0.01m];
% OUTPUTS
%   acem - average concentration over spring-neap cycle
%   asmt - proportion of time submerged over sp-np cycle
% NOTE
%   NB: result sensitive to time step and dmin
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
%    
    if nargin<5
        tsn  = 14.77;   delt = 1;   dmin = 0.01;  %default settings
    else
        tsn = options.tsn; delt = options.delt;  dmin = options.dmin;
    end
    %Initialise counters and sums
    scsn = 0; %initialise sum of marsh concentration over sp-np cycle
    sumT = 0; %intialise sum of time submerged over sp=np cycle
    dcflg= 0; %flag denotes downcrossing of water level relative to marsh
              %level. dcflg=0 water below, dcflg=1 water above marsh surface
    conc = 0; %time varying concentration over marsh
    sumc = 0;   %sum of concentration over individual up-crossings
    tide = simple_tide(wl,tsn,delt/60); %uses tsn in days and tint in minutes
    % ptime = tide.t/3600/24;    
    % figure;  plot(ptime,tide.z);
    nrec = length(tide.t);
    dcc = zeros(nrec,1);
    for it = 1:nrec
        if tide.z(it) > zm
            if tide.dz(it)>0
                cs = c0;          %on rising tide use conc in adjacent element
            else
                cs = conc;        %on falling tide use conc over marsh from previous time step
            end
            %
            flodep = tide.z(it)-zm;
            if flodep<dmin         %flow depth too small or flat is dry
                conc = 0;
            else
                dc = ((-wsm-tide.dz(it))*conc+cs*tide.dz(it))*delt/flodep;
                if conc+dc<0       %trap negative concentration
                    conc = 0;
                elseif dc>0        %increment concentration by change
                    conc = conc+dc;
                end
            end
            dcc(it,1) = conc;       %only used for test plot
            sumc = sumc+conc*delt;  %sum over individual up-crossings
            scsn = scsn+conc*delt;  %total concentration over sp-np cycle
            sumT = sumT+delt;       %total submerged time over sp-np cycle
            dcflg= 1;
        elseif tide.z(it)<zm && dcflg==1
            sumc = 0; conc = 0; dcflg = 0;             
        end
    end

    % figure; plot(ptime,dcc);
    Tsn = tsn*24*3600;
    acem = scsn/Tsn; %av. concentration over spring-neap cycle
    asmt = sumT/Tsn; %proportion of time submerged over sp-np cycle
end