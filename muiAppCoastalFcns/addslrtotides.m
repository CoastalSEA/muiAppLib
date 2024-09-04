function [t,var] = addslrtotides(wl,t,delta,exprate,pivotyear)
%
%-------function help------------------------------------------------------
% NAME
%   addslrtotides.m
% PURPOSE
%   Add sea level rise to water levels to get the combined signal
% USAGE
%   var = addslrtotides(wl,t,delta,exprate,pivotyear)
% INPUTS
%   wl - water level data
%   t - time,
%   delta - scale parameter
%   exprate - rate for exponential
%   pivotyear - (optional) is used to set slr to zero in a defined year and 
%               must be within the period defined by t.
% OUTPUT
%   t - input time is returned (e.g. for use in muiUserModel)
%   var - modified wl with slr added
% NOTES
%   Uses: slr(t)=delta.exp(exprate(t - 1900))where t is in calendar years. 
%   using delta = 0.001m and exprate = 0.011
%   provides a rate of 1mm/year prior to 1900, increasing to 
%   3mm/year by 2000, 5.2mm/year by 2050 and 9mm/year by 2100.
%   On this basis the total change from the year 2000 is 0.07m by 2020, 
%   0.2m by 2050 and 0.55m by 2100.
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%
startyear = datetime(1900,1,1,0,0,0);
tyears = years(t-startyear);
% integrate the rate equation from 1900 to get slr at each time interval
slr = delta/exprate*(exp(exprate*(tyears))-1);
if nargin>4
    %if a pivotyear is defined make slr vary relative to pivotyear
    pivotyear = datetime(pivotyear,1,1,0,0,0);
    pivotime = years(pivotyear-startyear);
    if pivotime<tyears(1) || pivotime>tyears(end)
        var = {'Pivot year must be within selected time period'};
        return;
    end
    idpivot = tyears==pivotime;   %assumes that pivotime is a value in tyears
    slrpivot = slr(idpivot);
else
    slrpivot = 0;
end
wlslr = wl+slr-slrpivot;
% if length(wlslr)==length(wl)
    %record is same length as wl so return variable
var = wlslr;
% else %within modelUI time is handled as a dimension so this is not needed?
%     %record is a different length so return time and variable
%     var = [t,wlslr]; 
% end
 
           