function tide = simple_tide(inp,rundur,tint)
%
%-------function help------------------------------------------------------
% NAME
%   simple_tide.m
% PURPOSE
%   Function to compute a tidal water level time series using main
%   constituents scaled to required tidal amplitude  
% USAGE
%   tide = simple_tide(inp,rundur,tint)
% INPUTS
%   inp - structure as defined below
%         MSL0                      mean tidel level relative to datum (mAD)
%         TidalAmp                  tidal elevation amplitude (m)
%         ElevPhase                 phase of elevation (ie k.x) (degs)
%         VelocityAmp               tidal velocity amplitude (m/s)
%         VelocityPhase             phase of velocity (ie k.x+phi) (degs)
%         M2amplitude               M2 tidal amplitude (m)
%         S2amplitude               S2 tidal amplitude (m)
%         O1amplitude               O1 tidal amplitude (m)
%   rundur - run duration (days)
%   tint - time interval for simulation (mins)
% OUTPUTS
%   tide - structure for time series containing
%        > time at Tintervals for Duration (s)
%        > tidal elevation relative to 0 (mAD)
%        > vertical velocity - rate of change of elevation dz/dt (m/s)
%        > horizointal velocity (m/s)
% NOTES
%   Velocity phase is kx+phi. and phi is the 
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    dur  = rundur*24*3600;           %run duration (seconds)
    tint = tint*60;                  %time interval for simulation (secs)
    z0   = inp.MSL0;                 %mean tidel level to ordnance datum (mOD)
    amp  = inp.TidalAmp;             %tidal elevation amplitude (m)
    pha1 = inp.ElevPhase*pi/180;     %phase of elevation (ie k.x) (rads)
    U    = inp.VelocityAmp;          %tidal velocity amplitude (m/s)
    pha2 = inp.VelocityPhase*pi/180; %phase of velocity (ie k.x+phi) (rads)
    aM2  = inp.M2amplitude;          %M2 tidal amplitude (m)
    aS2  = inp.S2amplitude;          %S2 tidal amplitude (m)
    aO1  = inp.O1amplitude;          %O1 tidal amplitude (m)

    fact = amp/(aM2+aS2+aO1);
    aM2s = fact*aM2;       %scaled amplitude of M2 constituent
    wM2  = 1.405e-4;       %angular speed of M2 constituent (rads/s)
    aS2s = fact*aS2;       %scaled amplitude of S2 constituent
    wS2  = 1.455e-4;       %angular speed of S2 constituent (rads/s)
    aO1s = fact*aO1;       %scaled amplitude of O1 constituent
    wO1  = wM2/2;          %angular speed of O1 constituent (rads/s)

    tide.t = (0:tint:dur)';
    tide.z = aM2s*cos(wM2*tide.t-pha1)+...
                     aS2s*cos(wS2*tide.t-pha1)+...
                     aO1s*cos(wO1*tide.t-pha1)+z0;

    tide.dz = -wM2*aM2s*sin(wM2*tide.t-pha1)...
                    -wS2*aS2s*sin(wS2*tide.t-pha1)...
                    -wO1*aO1s*sin(wO1*tide.t-pha1);

    tide.u = -(aM2s*sin(wM2*tide.t-pha2)+...
                     aS2s*sin(wS2*tide.t-pha2)+...
                     aO1s*sin(wO1*tide.t-pha2))*U/amp;

    % output_plot(tide);
end
%%
function output_plot(v) %#ok<DEFNU> 
    figure('Tag','PlotFig');
    tdays = days(seconds(v.t));
    subplot(3,1,1)
    plot(tdays,v.z);
    ylabel('Elevation (mAD)')
    subplot(3,1,2)
    plot(tdays,v.u)
    ylabel('Horizontal Velocity (m/s)')
    subplot(3,1,3)
    plot(tdays,v.dz)
    ylabel('Vertical velocity (m/s)')
    xlabel('Time (d)')
end




