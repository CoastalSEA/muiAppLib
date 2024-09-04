function theta = shore_orientation(sE,sN,bE,bN)
%
%-------function help------------------------------------------------------
% NAME
%   shore_orientation.m
% PURPOSE
%   Find the orientation for a series of coordinates that define a line
% USAGE
%   theta = shore_orientation(~,sE,sN,bE,bN)
% INPUTS
%   sE,sN - coordinates of the shoreline
%   bE,bN - coordinates of the baseline, landwards of shorelin 
% OUTPUT
%   theta - average of the current and next direction (ie the average 
%           either side of each profile point)
% NOTES  
%   Uses baseline to determine direction of shore and reverses line if the
%   approximately normal profile lines are not one quadrant in advance of
%   the shoreline quadrant. Return    
% SEE ALSO
%   sortENdata2line.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
   ang = NaN(size(sE));
   %find orientation based on relative position of baseline
   sgnE = sign(sE-bE);
   sgnN = sign(sN-bN);  

   qn1 = sgnE>=0 & sgnN>=0;  %0-90degrees     no adjustment needed
   qn2 = sgnE>=0 & sgnN<0;   %90-180 degrees  pi-ang
   qn3 = sgnE<0 & sgnN<0;    %180-270 degrees pi+ang
   qn4 = sgnE<0 & sgnN>=0;   %279-360 degrees 2pi-ang   

   qnorm = (mode(qn1*1+qn2*2+qn3*3+qn4*4));

   delE = diff(sE);          %easting differences along line
   delE = [delE(1),delE];    %pad the start of the line
   delN = diff(sN);          %northing differences along line
   delN = [delN(1),delN];    %pad the start of the line

   ql1 = delE>=0 & delN>=0;  %0-90 degTN     no adjustment needed
   ql2 = delE>=0 & delN<0;   %90-180 degTN   pi-ang
   ql3 = delE<0 & delN<0;    %180-270 degTN  pi+ang
   ql4 = delE<0 & delN>=0;   %270-360 degTN  2pi-ang

   qline = mod((mode(ql1*1+ql2*2+ql3*3+ql4*4)+1),5);
   if qline==0, qline=1; end

   if qline~=qnorm
       revE = fliplr(sE);
       revN = fliplr(sN);

       delE = diff(revE);        %easting differences along line
       delE = [delE(1),delE];    %pad the start of the line
       delN = diff(revN);        %northing differences along line
       delN = [delN(1),delN];    %pad the start of the line

       ql1 = delE>=0 & delN>=0;  %0-90 degTN     no adjustment needed
       ql2 = delE>=0 & delN<0;   %90-180 degTN   pi-ang
       ql3 = delE<0 & delN<0;    %180-270 degTN  pi+ang
       ql4 = delE<0 & delN>=0;   %270-360 degTN  2pi-ang
   end

   en = delE./delN;              %ratio of change in eastings/northings   
   ang(ql1) = rad2deg(atan(en(ql1)));
   ang(ql2) = rad2deg(pi()+atan(en(ql2)));
   ang(ql3) = rad2deg(pi()+atan(en(ql3)));
   ang(ql4) = rad2deg(2*pi()+atan(en(ql4))); 

   if any(ql1) && any(ql4)        %check whether there are points +/-360
       q360 = mode(ql1*1+ql4*4);  %find the dominant quadrant
       if q360==1                 %add or subract 360
           ang(ql4) = ang(ql4)-360;
       else
           ang(ql1) = ang(ql1)+360;
       end
   end
   theta = movmean(ang,[0,1]);    %moving average of current and next angle
                                  %this is a mean of the angles either side
                                  %of each point
   theta(theta<0) = theta(theta<0)+360;     %restore the quadrants if
   theta(theta>360) = theta(theta>360)-360; %changed by averaging    
end 