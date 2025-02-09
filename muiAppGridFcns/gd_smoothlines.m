function smoothlines = gd_smoothlines(lines,method,window,degree,npnts)
%
%-------function help------------------------------------------------------
% NAME
%   gd_smoothlines.m
% PURPOSE
%   smooth one or more line segments using either moving average or the
%   Savitzky-Golay method
% USAGE
%   smoothlines = gd_smoothlines(lines,method,win,deg,npnts)
% INPUTS
%   lines - struct of x,y vectors to be smoothed. NaN used as line separator
%           x and y need to be row vectors of points [2xN]
%   method - 'movmean' for moving average; or 'sgolay' for Savitzky-Golay smoothing;
%   window - window size to use for moving average (1 or 2 elements)
%   degree - Savitzky-Golay degree (<window)
%   npnts - minimum number of points required to apply smoothing
% OUTPUTS
%   smoothlines - input lines with the smoothing applied
% NOTES
%    uses smoothdata function. default is to 'omitnan'
% SEE ALSO
%   called in gd_boundary
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
% 
    idN = [0,find(isnan(lines.x))];
    smlines = [];
    %find each line and smooth
    for i=1:length(idN)-1
        aline = lines.y(idN(i)+1:idN(i+1)); %include trailing NaN
        if length(aline)>npnts
            if strcmp(method,'sgolay')
                aline = smoothdata(aline,method,window,'Degree',degree);
            else %use movmean
                aline = smoothdata(aline,method,window);
            end
        end
        smlines = [smlines,aline(1:end-1),NaN]; %#ok<AGROW>
    end
    smoothlines.x = lines.x;
    smoothlines.y = smlines;
end
