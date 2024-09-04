function [xdist,slope,xs,zs] = slope_points(x,z,zlevel,delz,xmin)
%
%-------function help------------------------------------------------------
% NAME
%   slope_points.m
% PURPOSE
%   Find point and slope on line (eg a shore profile) nearest to zlevel
% USAGE
%   [xs,zs,xdist,slope] = slope_points(x,z,zlevel,delz,xmin)
% INPUTS
%   x - vector or matrix of horizontal distance, each row is a section-line
%   z - vector or matrix vertical data (y-axis), each row is a section-line
%   zlevel - vertical level to be found on section-line
%   delz - +/- offset distance relative to zlevel to uuse to compute slope
%   xmin - minimum x value to start search (use only x>xmin)
% OUTPUT
%   xdist - distance to zlevel on section-line
%   slope - slope at zlevel over points that are +/-delz
%   xs - array of distances used for each section line (zlevel and +/-zdel)
%   zs - array of levels used for each section line (zlevel and +/-zdel)
%   3 rows for each section line ie xs(3,nrec) and zs(3,nrec)  
% SEE ALSO
%   used in CT_BeachAnalysis.m function getShorelinePositions
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%        
    nstep = size(x,1);            
    xdist = NaN(nstep,1); slope = xdist; xs = NaN(3,nstep); zs = xs;  
    for i=1:nstep
        xi = x(i,:);
        zi = z(i,:);
        [~,Ix] = min(abs(xi-xmin));
        [~,Ix0] = max(zi(Ix:end)); %crest seawards of ymin
        Ix0 = Ix+Ix0-1;
        [~,Iz0] = min(abs(zi(Ix0:end)-zlevel));       %index for zlevel          
        [~,Imn] = min(abs(zi(Ix0:end)-zlevel-delz));  %index for -delz
        [~,Imx] = min(abs(zi(Ix0:end)-zlevel+delz));  %index for +delz
        idx = [Imn+Ix0-1,Iz0+Ix0-1,Imx+Ix0-1];
        xs(:,i) = xi(idx);   %distances of selected points (-delz,z0,+delz)
        zs(:,i) = zi(idx);   %levels of selected points (-delz,z0,+delz)
        idd = Imn+Ix0-1:1:Imx+Ix0-1;
        zdata = zi(idd);
        xdata = xi(idd);
        %check for duplicates
        [~,ia,ib]=unique(zdata);
        icz=setdiff(ib,ia);
        while ~isempty(icz)
            tol = 0.0001;
            zdata(icz) = zdata(icz)+tol;
            [~,ia,ib]=uniquetol(zdata,tol);
            icz=setdiff(ib,ia);
        end
        %interpolate to find point
        if length(zdata)<2
            xdist(i) = NaN;  slope(i) = NaN;
        else
            try
                xdist(i) = interp1(zdata,xdata,zlevel,'linear');
            catch
                wtxt = sprintf('Unable to interpolate profile %d in "slope_points"',i);
                warndlg(wtxt);
                return;
            end
            slope(i) = -1/mean(gradient(zs(:,i),xs(:,i)));
        end
    end
end