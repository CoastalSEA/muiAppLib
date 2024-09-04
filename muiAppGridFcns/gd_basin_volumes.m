function hyps = gd_basin_volumes(hdst,isplot)
%
%-------function help------------------------------------------------------
% NAME
%   gd_basin_volumes.m
% PURPOSE
%   compute area and volume totals over x-z from hypsometry of gridded elevation data
% USAGE
%   hyps = gd_basin_volumes(hdst,isplot)
% INPUTS
%   hdst - dstable of surface area datagenerated by gd_basin_volumes.m
%   isplot - plot outputs, plot if true (optional, default is no plot)
% OUTPUT
%   hyps - zcentre, xdist, zsurf, zvol, xzsurf, xzvol as a struct, where
%          zcentre: vertical elevation of interval mid-point
%          xdist: distance along x-axis
%          zsurf: total surface plan area for each vertical interval, z
%          zvol: total volumes below each vertical interval, z
%          xzsurf: cumulative area from head and bed, as function of x and z
%          xzvol: cumulative area from head and bed, as function of x and z
%   isplot - true if plots are to be produced
% NOTES
%   uses the dstable of plan areas as a function of x and z, created using
%   gd_basin_hypsometry. Function used to explore methods of estimating the
%   prism in ms_userfunction. 
% SEE ALSO
%   gd_basin_hypsometry which outputs area and volume for whole channel
%   and uses gd_basin_properties which uses the same methods as this
%   function.
%
% Author: Ian Townend
% CoastalSEA (c) July 2022
%--------------------------------------------------------------------------
%
    if nargin<2, isplot = false; end

    hyps.zcentre = hdst.Dimensions.Z;
    hyps.xdist = hdst.Dimensions.X;
    zhist = squeeze(hdst.SArea);
    histint = hdst.UserData.histint;
    %delx = abs(hyps.xdist(2)-hyps.xdist(1));
    
    zarea = sum(zhist,1); %sum surface area over x to give total area at each z
    hyps.zsurf = cumsum(zarea); %cumulative surface area at or below each elevation z
    %total volume below each elevation z, trapezoidal summation handles ends
    hyps.zvol = cumtrapz(histint,hyps.zsurf);
    
    %cumulative surface area at or below each elevation z for each x inteval
    %such that hyps.zsurf = sum(hyps.xzsurf,1);
    hyps.xzsurf = cumsum(zhist,2);
    %cumulative volume at or below each elevation z for each x inteval
    %such that hyps.zvol = sum(hyps.xzvol,1);
    hyps.xzvol = cumtrapz(histint,hyps.xzsurf,2);

    if isplot
        plot_area_volume(hdst,hyps) %check plot
    end
end

%%
function plot_area_volume(hdst,hyps)
    %check plots of area and volume outputs
    hf1 = figure('Name','Hypsometry','Tag','PlotFig');
    ax1 = axes(hf1);
    SA = squeeze(hdst.SArea)'; 
    SA(SA==0) = NaN;
    surf(ax1,hdst.Dimensions.X,hdst.Dimensions.Z,SA);
    hc = colorbar;
    xlabel('Distance (m)');
    ylabel('Elevation (mAD)')
    zlabel('Surface Area (m^2)')
    hc.Label.String = 'Surface Area (m^2)';
    title('Surface area distribution')
    
    %conventional full basin plots of hypsometry
    hf2 = figure('Name','Hypsometry','Tag','PlotFig');
    ax2 = axes(hf2);
    plot(ax2,hyps.zsurf,hyps.zcentre,'DisplayName','Surface area');
    hold on
    plot(ax2,hyps.zvol,hyps.zcentre,'DisplayName','Volume');
    hold off
    xlabel('Volume and Surface Area');
    ylabel('Elevation (mAD)');
    legend
    title('Basin hypsometry')

    %plots of area and volume as function of x and z, cumulative from head
    figure('Name','Hypsometry','Tag','PlotFig');
    sa1 = subplot(1,2,1);
    contourPlot(sa1,hdst,hyps.xzsurf,'Surface Area (m^2)')

    sa2 = subplot(1,2,2);
    contourPlot(sa2,hdst,hyps.xzvol,'Volume (m^3)')
    sgtitle('Cumulative variation from head and bed')
end
%%
function contourPlot(ax,hdst,var,plotvar)
    if isgraphics(ax,'figure')
        ax = axes(ax);
    end
     var(var==0) = NaN;
    contourf(ax,hdst.Dimensions.X,hdst.Dimensions.Z,var');
    hc = colorbar(ax);
    xlabel('Distance (m)');
    ylabel('Elevation (mAD)')
    zlabel(plotvar)
    hc.Label.String = plotvar;
end