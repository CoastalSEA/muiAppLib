function Fch = eff_fetch(F_len, F_dir, Dir, dirn, eflg)
% 
%-------function help------------------------------------------------------
% NAME
%   eff_fetch.m
% PURPOSE
%	Compute the effective fetch for each mean direction based on the 
%   directional distribution function 
% USAGE
%   Fch = eff_fetch(F_len, F_dir, Dir, dirn, eflg)
% INPUTS
%   F_len - Fetch length array (m)
%   F_dir - Fetch direction array (deg) NB: code assumes EQUAL INTERVALS
%   Dir   - Mean wind direction (deg)
%   dirn  - Directional index
%   eflg  - eflg=0 use Donelan; eflg=1 use SPM; eflg=-1 use fetch along 
%           wind direction
% RESULTS
%   Fch = Effective fetch length
% NOTES
%   Donelan, Hamilton and Hui, R.Soc, 1985, A315, 509-562
%   US Army Corps, Shore Protection Manual, 1984
% SEE ALSO
%   used in tma_calc.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2016
%---------------------------------------------------------------------------
%
if eflg<0
    Fch = interp1(F_dir, F_len, Dir); %interpolate fetch for wind direction
    return;
end

beta = 2.28;                %frequency dependent coefficient (from Donelan et al)
gam  = 1;                   %SPM scaling coefficient
%
Dir0 = Dir-90;              %start angle; offset from mean direction by 90 degrees
Dir0 = Dir0 + 360*(Dir0<0); %correct if less than 0 degrees 
% 
dphi = F_dir(2)-F_dir(1);   %interval for fetch directions
%
% convert to radians
d2r = pi/180;
Dir=Dir*d2r; Dir0=Dir0*d2r; dphi=dphi*d2r; F_dir=F_dir*d2r;
%correct for Phi greater than 360. 
Phi  = mod(Dir0:dphi:Dir0+pi,2*pi);

% handle angles greater than provided in F_dir
if any(Phi>max(F_dir))
    n_dirs = ceil((max(Phi)- max(F_dir))/dphi);
    F_len = cat(1,F_len,zeros(n_dirs,1));
    F_add = zeros(n_dirs,1);
    for i=1:n_dirs
        F_add(i)=max(F_dir)+i*dphi;
    end
    F_dir = cat(1,F_dir,F_add);
end
%NB interp1 requires monotonic input
Flen = interp1(F_dir, F_len, Phi);
ndir = length(Dir);
Fch = zeros(ndir,1);
for i=1:ndir
    ang =Phi-Dir(i);
    %
    if eflg==0              %use hyperbolic secant as per Donelan
        Fch(i) = sum(Flen.*(sech(beta.*(ang))).^dirn);
        Fch(i) = Fch(i)./sum((sech(beta.*(ang))).^dirn);
    else                    %use cosine distribution as per SPM
        Fch(i) = sum(Flen.*(cos(gam.*(ang))).^dirn);
        Fch(i) = Fch(i)./sum((cos(beta.*(ang))).^dirn);
    end
end


