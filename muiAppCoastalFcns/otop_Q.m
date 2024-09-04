function Q = otop_Q(swl,Hs,Tz,alp,bs,g,structure,isquiet)
%
%-------function help------------------------------------------------------
% NAME
%   otop_Q.m
% PURPOSE
%   Function to calculate the overtopping discharge for a simple sloping
%   structure based on the Method proposed by Owen
% USAGE
%   Q = otop_Q(swl,Hs,Tz,alp,bs,g,structure,isquiet)
% INPUTS
%   hydraulic conditions
%   swl - still water level (m above datum)
%   Hs  - incident significant wave height (m)
%   Tz  - mean zero crossing wave period (s)
%   alp - wave angle in degrees from normal to the wall
%   bs  - nearshore bed slope (1:m)
%   g   - acceleration due to gravity (m/s^2)
%   structure - struct that defines form of sea wall
%   isquiet - flag to show waitbar - show if false (optional)
% RESULTS
%   Q - overtopping discharge in m3/s/m-run of seawall  
% NOTES
%   Ref: Owen M W, 1980, Design of seawalls allowing for wave overtopping, 
%   Report No: EX 924, H R Wallingford.
%   TR W178 provides method to include crest width and wave return walls -
%   **not yet implemented*****************************************
% SEE ALSO
%   otop_C to get crest level for a given overtopping discharge
%   coeff_AB for model coefficients
%
% Author: Ian Townend
% CoastalSEA (c)June 2016
%--------------------------------------------------------------------------
%
if nargin<8
    isquiet = false;
end
% read in definition of structure from struct
cl=structure.cl;       %crest level (m above datum)
cw=structure.cw;       %crest width (m)
uws=structure.uws;     %upper wall slope (1:m)
bl=structure.bl;       %berm level (m above datum)
bw=structure.bw;       %berm width (m)
lws=structure.lws;     %lower wall slope (1:m)
tl=structure.tl;       %toe level (m above datum)
r=structure.r;         %wall roughness

% crude check for depth limiting conditions
if ~isquiet, hw = waitbar(0,'Processing'); end
dep = swl-tl;
hbflag = 3; %specific case in hb_break for overtopping, which uses Tz
Hb = hb_break(dep,bs,Tz,g,hbflag);
Hsb = min(Hb,Hs,'includenan');

if ~isquiet, waitbar(0.2); end
nrec = length(Hs);
Q = zeros(nrec,1); A = Q; B = Q;
for i=1:nrec 
    % define effective wall slope (ws)
    bd = swl(i)-bl;
    cws = ((bl-tl)*lws+(cl-bl)*uws+bw)/(cl-tl);
    if bd<0 && bw>0           %water below berm level and berm exists
        %use lower wall slope if berm exists and effective slope>5
        %should probably check effective freeboard for each water level
        if cws>5, ws = lws; else, ws = cws; end 
    elseif bw==0 && uws~=lws  %no berm width use av.upper + lower wall slope
        ws = cws;
    elseif bd>4               %berm too deep to be effective
        ws = cws;  bw = 0; bd = 0; 
    else                      %all other cases use lower wall slope
        ws = lws;
    end
    %
    if ws>5 || bw>80
        warndlg('Wall slope or berm width out of range');
        Q = [];
        close(hw)
        return;
    end
    % get wall coefficients A and B
    [A(i),B(i)] = coeff_AB(ws,bw,bd);  
    if ~isquiet, waitbar(0.2+i/nrec*0.7); end
end

% correct coefficients for wave angle
alp = abs(alp);
A = A.*(1+alp.*(0.0781-alp.*(0.0031-2.7995*10^-5*alp)));
B = B.*(1-alp.*(0.0001+alp.*(1.9349*10^-5-2.1016*10^-6*alp)));
if ~isquiet, waitbar(1); end
% calculate overtopping rate
Rc = cl-swl;    %structure freeboard
Rs = Rc./(Tz.*sqrt(g*Hsb));
Qs = A.*exp(-(B./r).*Rs);
Qo  = g*Qs.*Tz.*Hsb;    %overtopping discharge (m3/s/m-run of seawall)
Q = Qo.*(Qo>1e-6);
if ~isquiet, close(hw); end
