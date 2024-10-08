function [slr,dslr] = sealevelrise(yrs,pivotyear,dslrvec,option)
%
%-------function help------------------------------------------------------
% NAME
%   sealevelrise.m
% PURPOSE
%   Function to compute sea level rise 
% USAGE
%   slr = sealevelrise(time,datetime,dslr,option)
% INPUTS
%   yrs - vector of time in years
%   pivotyear - pivot year for zero sea level rise (ie negative for years 
%               before pivot year) - note: this is not the same as yr0
%   dslrvec - vector containing rate parameters
%              linear - rate of sea level rise (m/yr)
%              exponential - [dslr, yr0, dslr0], where:
%                  dslr - exponential rate increase (m/yr), 
%                  yr0 - start of exponential (year), 
%                  dslr0 - linear rate pre start of exponential
%              Holocene+Modern double exponential with parameters
%              [dslr, yr0, dslr0, hscale, hshape, hoffset], where:
%                  first 3 as defined above and
%                  hscale, hshape, hoffset - define an exponential going
%                  back to circa 6000 BC. Values of 1.32,2400,0.6 reproduce
%                  the Humber Holocene record to 1900 with reasonable
%                  agreement.
%                  Typical values are: [0.011 1900 0.0005 1.32 2400 0.6]
%              if dslrvec is empty default values used [0.011,1900,0.001] 
%   option - determine how slr is calculated
%              1 - linear slr
%              2 - exponential from yr0, linear before
%              3 - Holocene exponential to yr0 and then option 2
%              4, etc - user can add additional options as required
% OUTPUTS
%   slr = magnitude of sea level rise from start time to time now (m)
%   dslr - rate of sea level rise at time t (m/yr)
%          only changes if rate is not linear, otherwise = dslrvec
%          exponential rate is the gradient +/-1 year about time t
% EXAMPLE
%   yrs = 1800:1:2100;
%   [slr,dslr] = sealevelrise(yrs,2000,[],2); %exponential slr from 1900,
%                                             %with zero in the year 2000 
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    switch option
        case 1    %linear rate of slr
            elapsedyears = yrs-pivotyear;
            slr = dslrvec*elapsedyears;
            dslr = dslrvec;
            
        case {2,3}    %exponential rate of slr
            if isempty(dslrvec) || length(dslrvec)<3
                %exponential increase post 1900 that approximates the Defra guidance
                %slr = exp(etr*T) post 1900, where T is in years 
                %For a value of dslr0 = 1mm/yr: 
                %  slr is 1mm/yr in 1900; 3mm/yr in 2000; 5mm/yr in 2050; 
                %  9mm/yr in 2100 
                expslr = 0.011;      %exponential rate of increase
                yr0 = 1900;          %year in which exponential increase starts (yr)
                dslr0 = 0.001;       %linear rate of slr pre 1900
            else
                expslr = dslrvec(1); %exponential rate of increase
                yr0 = dslrvec(2);    %year in which exponential increase starts (yr)
                dslr0 = dslrvec(3);  %linear rate of slr pre 1900
            end
            
            if option==2
                slr = combinedLinExp(yrs,yr0,dslr0,expslr);
                slrpivot   = combinedLinExp(pivotyear,yr0,dslr0,expslr);                
            else
                slr = combinedHoloceneModern(yrs,yr0,dslr0,expslr,dslrvec);                                                        
                slrpivot = combinedHoloceneModern(pivotyear,yr0,dslr0,...
                                                        expslr,dslrvec);                                    
            end
            %offset the curve so slr=0 in the pivotyear
            slr = slr-slrpivot;  
            
            %compute the rate of slr at time t
            if isscalar(yrs) %single year requested
                gradyr = [yrs-1,yrs,yrs+1];
                gradslr = combinedLinExp(gradyr,yr0,dslr0,expslr);
                dslr = mean(diff(gradslr));
            else             %vector of time values
                dslr = diff(slr)./diff(yrs);
                dslr = [dslr(1),dslr];
            end
            
        case 4    %user option
            %insert user code
            %no slr for 100 years and then 5mm/yr from year 100
            pivotyear = 100;
            dslr = 0;
            newrate = 0.005;
            idx = yrs>pivotyear;
            slr = yrs*dslr;
            if any(idx)
                slr = (yrs-pivotyear)*newrate*idx;
                dslr = newrate;
            end

        otherwise
            slr = []; dslr = [];
            warndlg('Undefined option in sealevelrise')            
    end
end
%%
function slr = combinedLinExp(yrs,yr0,dslr0,expslr)
    %compute slr for the combined linear and exponential variation

    %compute any linear component of slr
    slr = zeros(size(yrs));            
    idx = yrs>yr0;   %index of records after switch to exponential
    slr0 = dslr0*(yr0-yrs(1));   %starting value based on offset from yr0 
    slr(~idx) = dslr0*(yrs(~idx)-yrs(1))-slr0; %linear slr

    %compute the exponential
    slr(idx) = dslr0/expslr*(exp(expslr*(yrs(idx)-yr0))-1);
end
%%
function slr = combinedHoloceneModern(yrs,yr0,dslr0,expslr,dslrvec)
    %compute slr for the combined double exponential variation
    % dslrvec = [expslr,yr0,dslr0,1.32,2400,0.6] for the Humber Holocene
    holoscale = dslrvec(4); holoshape = dslrvec(5); holoffset = dslrvec(6);
    %compute any Holocene component of slr
    slr = zeros(size(yrs));            
    idx = yrs>yr0;   %index of records after switch to exponential
    slr(~idx) = -(holoscale*exp(-yrs(~idx)/holoshape)-holoffset); %inverse exponential

    %compute the exponential
    if any(idx)
        slr(idx) = dslr0/expslr*(exp(expslr*(yrs(idx)-yr0))-1);
    end
end
            