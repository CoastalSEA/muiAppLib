function ratio = beachtransportratio(Diri,theta,isvector)
%
%-------function help------------------------------------------------------
% NAME
%   beachtransportratio.m
% PURPOSE
%   Compute the ratio of onshore to alongshore wave action given by tan(alp)
% USAGE
%   ratio = beachtransportratio(Diri,theta,is)
%INPUT
%   Diri - inshore wave directions (degTN)
%   theta - bed contour (degTN)
%   isvector - true if alongshore direction is to be retained in output (optional)
%OUTPUT
%   ratio - tan(alpha) where alpha is angle between wave crest and bed contour 
%           note: ratio ignores direction along shore (i.e. always +ve)
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%----------------------------------------------------------------------
%
    if nargin<3, isvector = false; end
    %static variables 
    rads = pi()/180;
    
    prompt = {'Max angle from normal (deg)- eg 85:'};
    dlgtitle = 'Max wave angle';
    numlines = 1;
    defaultvalues{1} = num2str(85);
    useInp=inputdlg(prompt,dlgtitle,numlines,defaultvalues);
    if isempty(useInp), ratio = 'User cancelled'; return; end %user cancelled
    maxangle = str2double(useInp{1})*rads;
    
    %find the angle between the wave crest and the bed contour (always positive
    %so direction along shore is lost).
    [alpi,quads] = getalp(Diri,theta(1));
    alpi(alpi>maxangle) = NaN;
    
    ratio = nan(size(alpi));  %to assign to timeseries, return a column vector
    if isvector
        %to retain direction with left to right positive when looking from offshore
        idon_neg = quads==1;
        ratio(idon_neg) = tan(-alpi(idon_neg));
        idon_pos = quads==2;
        ratio(idon_pos) = tan(alpi(idon_pos));
    else    
        idon = quads==1 | quads==2;
        ratio(idon) = tan(alpi(idon)); 
    end
end
