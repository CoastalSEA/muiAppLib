function [Hsi,Diri]= refraction(Hs0,Tp,Dir0,depths,theta,Kf,isshore)
%
%-------function help------------------------------------------------------
% NAME
%   refraction.m
% PURPOSE
%   plane bed wave refraction and shoaling using linear wave theory
% USAGE
%   [Hsi,Diri] = refraction(Hs0,Tp,Dir0,[dep1,dep2],theta,Kf,isshore)
% INPUTS
%   Hs0  - offshore significant wave height (m)
%   Tp   - peak wave period (s)
%   Dir0 - wave direction (degrees TN)
%   depths -  water depths (m): minimum is [dep0,depi]  - can be matrix dep(i,n)
%   theta- angle of contours from north (degrees TN) - can be array theta(n)
%   Kf   - friction coefficient (default=1)
%   isshore - flag to remove offshore propagating waves
% RESULTS
%   Hsi - inshore wave height
%   Diri - inshore wave angle
% NOTES
%   USACE Shore Protection Manual
%   NB changed Feb 2019 to pass in a matrix of depths and bed contours 
%   (rows are cases (eg time steps) and columns are the number of spatial
%   intervals)
%   Minimum input required is dep = [d0,d1] and theta = theta1.
% SEE ALSO
%   calls celerity.m, getalp.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%---------------------------------------------------------------------- 
% Test cases
% Case 1: Hs0 = 1; Tp=8.5; Dir0=210; dep0=10; depi=3; theta=60;
% results: alpi=30.5;Kr=0.76;Ks=1.22;Hsi=0.93;Diri = 180.5 if Kf=1
% [Hsi,Diri] = refraction(1,8.5,210,[10,3],60,1,false);
% Case 2: theta = 320: 
% [Hsi,Diri] = refraction(1,8.5,130,[10,3],320,1,false);
% results: Hsi = 0.56; Diri = 85.3
%--------------------------------------------------------------------------
    
    Hsi = []; Diri = [];
    %handle various input styles
    ndepth = size(depths,2);  %number of depth intervals per time step
    if ndepth<2
        warndlg('Error in refraction. Minimum of two depths required');
        return;
    end
    ntheta = size(theta,2); %number of angles (must be same as ndepth but constant over time)
    if ntheta<2
        theta(2) = theta;   %only inshore value supplied (set offshore value to be same)
    end
    
    %handle multi-contour specification using n-depths and n-thetas
    %or set up intervals between start and end values
    if ndepth>2  
        if ndepth~=ntheta
            txt = sprintf('Error in refraction\nMult-contour specification requires the same number of depths and angles');
            warndlg(txt);
        else
            %nint = ndepth - the number of user specified intervals
        end
    elseif theta(1)==theta(2)
        %nint = 2 and the start and end values are the same
        if any(depths<0)              %fzero in hs_surf can call function 
            Hsi = 0; Diri = 0;        %with negative depths
            return;
        end
    else     
        %interpolate specified start and end values to nint increments
        nint = 5; %number of intervals used to iterate refraction solution
        depths = getdepthinc(depths,nint);   
        theta = theta(1):(theta(end)-theta(1))/nint:theta(end);   
    end

    % inshore wave parameters
    [Kr,Ks,Diri,quad] = getwavecoeffs(Tp,depths,Dir0,theta);
    Hsi = Hs0.*Kr.*Ks.*Kf;

    %if wave is propagating onshore to a shoreline, remove offshore going waves
    %littoraldriftstats uses qs=0 to identify offshore going waves
    if isshore
        %remove waves in quadrants 3 and 4
        idx = quad==3 | quad==4;
        Diri(idx) = NaN;
        Hsi(idx) = 0;       
    end
end
%%
function phi = getdirection(alp,theta,quad)
    %find the wave direction, phi, for given wave crest to contour angle
    %alp is the angle between the wave crest and the contour
    %theta is the angle of the bed contour (degTN) - looking seaward
    %quad is the quadrant of the input wave direction relative to contour
    alp = alp*180/pi();
    phi = NaN(size(alp));
    for i=1:length(alp)
        switch quad(i)
            case 1
                phi(i) = theta+90-alp(i);
            case 2
                phi(i) = theta+90+alp(i);
            case 3
                phi(i) = theta+270-alp(i);
            case 4
                phi(i) = theta+270+alp(i);
            otherwise
                %phi is NaN
        end
    end
    idn = ~isnan(phi);
    phi(idn) = mod(phi(idn),360);   %arguments must be real
end
%%
function depths = getdepthinc(depth2,nint)
    %set up array of intervals for a vector of start and end depths (end=2)
    %depth2 - nx2 array of start and end depths for n time iintervals
    %nint - number of intervals between start and end
    m = size(depth2,1);
    depths = NaN(m,nint+1);
    %ensure that depths are not negative and not the same
    depth2(depth2<0.2) = 0.2;      
    dep1 = depth2(:,1);
    dep2 = depth2(:,2);   
    %find records that are not nan
    [idr,~] = find(~isnan(dep1) & ~isnan(dep2)); 
    for i=idr'
        %define intervals for each depth range in time series
        if dep1(i)==dep2(i), dep2(i) = dep1(i)+0.5; end %force a difference
        dinc = (dep2(i)-dep1(i))/nint;
        depths(i,:) = dep1(i):dinc:dep2(i); 
    end 
end
%%
function [Kr,Ks,dirs,quad,alpi] = getwavecoeffs(Tp,depths,Dir0,theta)
    %compute the refraction coefficient for given wave period and direction
    %as the wave travels from depth(i)->depth(i+1) and theta(i)->theta(i+1)
    %Tp - array pf wave periods (s) for n time intervals
    %depths - nxm array of two or more depths for each time interval
    %Dir0 - array pf wave direction (degTN) for n time intervals
    %theta - angle of bed contour to match number of depths
    dirs = Dir0;
    [alp0,quad] = getalp(Dir0,theta(1));    
    nint = size(depths,2);    
    for i=2:nint
        [alpo,quad] = getalp(dirs,theta(i));  
        c0 = celerity(Tp,depths(:,i-1));
        ci = celerity(Tp,depths(:,i));    
        alpi = getalpin(c0,ci,alpo);
        dirs = getdirection(alpi,theta(i),quad);
    end
    Kr = getrefcoeff(alp0,alpi);
    Ks = getshoalcoeff(Tp,c0,ci,depths);
end
%%
function Kr = getrefcoeff(alp0,alpi)
    %refraction coefficient, Kr, using Snells law to obtain 
    %alp0 - start angle between wave crest and contour (rads)
    %alpi - end angle between wave crest and contour (rads)
    Kr = sqrt(cos(alp0)./cos(alpi));
    %set Kr=1 for waves with wave crest close to normal to contour
    %1.4-1.7 rads (~80-100deg) and 4.5-4.9 (~260-280deg)
    %note that getalpin uses the combination of ci/c0 and sin(alpo) to
    %identify conditions that become complex, and sets alpi=pi/2 or 3pi/2,
    %hence these all fall within this range
    ida = (alpi>1.4 & alpi<1.75) | (alpi>4.5 & alpi<4.9);
    Kr(ida) = 1;
end
%%
function Ks = getshoalcoeff(Tp,c0,ci,depth)
    %shoaling coefficient, Ks (linear wave theory)
    %Tp - vector of wave period at n time intervals 
    %c0,ci - vectors of celerity for start and end depths at n time intervals 
    %depth - nx2 array of start and end depths
    pidL0 = 4*pi()*depth(:,1)./(c0.*Tp);
    cg0 = c0/2.*(1+pidL0./sinh(pidL0));
    pidLi = 4*pi()*depth(:,end)./(ci.*Tp);
    cgi = ci/2.*(1+pidLi./sinh(pidLi));
    Ks  = sqrt(cg0./cgi);
end
%%
function alpi = getalpin(c0,ci,alpo)
    %compute the end angle, alpi, between wave crest and contour (rads)
    %c0,ci - vectors of celerity for start and end depths at n time intervals 
    %alpo - start angle between wave crest and contour (rads) 
    %quad - quadrant of the wave direction relative to the contour    
    var = ci./c0.*sin(alpo);
    idx = var<=1;  %asin is complex for values >1, oblique wave angles
    alpi = NaN(size(var));
    alpi(idx) = asin(var(idx)); 
end