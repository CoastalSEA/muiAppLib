function tide = simple_tide(inp)
%
%-------function help------------------------------------------------------
% NAME
%   simple_tide.m
% PURPOSE
%   Function to compute a tidal water level time series using main
%   constituents scaled to required tidal amplitude  
% USAGE
%   tide = simple_tide(inp)
% INPUTS
%   inp - structure as defined below
%         Duration                  duration (days)
%         Tinterval                 time interval for simulation (mins)
%         z0                        mean tidel level to ordnance datum (mOD)
%         TidalAmp                  tidal elevation amplitude (m)
%         ElevPhase                 phase of elevation (ie k.x) (deg)
%         VelocityAmp               tidal velocity amplitude (m/s)
%         VelocityPhase             phase of velocity (ie k.x+phi) (deg)
%         aM2                       M2 tidal amplitude (m)
%         aS2                       S2 tidal amplitude (m)
%         aO1                       O1 tidal amplitude (m)
% OUTPUTS
%   tide - structure for time series containing
%        > time at Tintervals for Duration
%        > tidal elevation relative to 0 mOD
%        > vertical velocity - rate of change of elevation dz/dt (m/s)
%        > horizointal velocity (m/s)
%
% Author: Ian Townend
% CoastalSEA (c)June 2020, corrected Nov 2023
%--------------------------------------------------------------------------
%
    Dur  = inp.Duration*24*3600;       %duration (seconds)
    tint = inp.Tinterval*60;           %time interval for simulation (secs)
    z0   = inp.z0;
    amp  = inp.TidalAmp;
    pha1 = inp.ElevPhase*pi/180;
    U    = inp.VelocityAmp;
    pha2 = inp.VelocityPhase*pi/180;
    aM2  = inp.aM2;        %M2 tidal amplitude (m)
    aS2  = inp.aS2;        %S2 tidal amplitude (m)
    aO1  = inp.aO1;        %O1 tidal amplitude (m)
    
    fact = amp/(aM2+aS2+aO1);
    aM2s = fact*aM2;       %scaled amplitude of M2 constituent
    wM2  = 1.405e-4;       %angular speed of M2 constituent (rads/s)
    aS2s = fact*aS2;       %scaled amplitude of S2 constituent
    wS2  = 1.455e-4;       %angular speed of S2 constituent (rads/s)
    aO1s = fact*aO1;       %scaled amplitude of O1 constituent
    wO1  = wM2/2;          %angular speed of O1 constituent (rads/s)
    
    tide.t = (0:tint:Dur)';
    tide.z = aM2s*cos(wM2*tide.t-pha1)+...
                     aS2s*cos(wS2*tide.t-pha1)+...
                     aO1s*cos(wO1*tide.t-pha1)+z0;
    
    tide.dz = -wM2*aM2s*sin(wM2*tide.t-pha1)...
                    -wS2*aS2s*sin(wS2*tide.t-pha1)...
                    -wO1*aO1s*sin(wO1*tide.t-pha1);
    
    tide.u = -(aM2s*sin(wM2*tide.t-pha2)+...
                     aS2s*sin(wS2*tide.t-pha2)+...
                     aO1s*sin(wO1*tide.t-pha2))*U/amp;
end





