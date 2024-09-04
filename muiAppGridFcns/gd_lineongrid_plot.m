function hfig = gd_lineongrid_plot(grid,cline,plotxt)
%
%-------function help------------------------------------------------------
% NAME
%   gd_lineongrid_plot.m
% PURPOSE
%   plots a defined line onto a countour or surface plot of a grid (e.g a
%   channel centre-line). If z values are included the plot is a 3D plot.
% USAGE
%   gd_lineongrid_plot(grid,cline,plotxt)
% INPUTS
%   grid - struct of the source grid provided by getGrid in GDinterface 
%   cline - x,y,z struct of centre-line path to use for transformation
%   plotxt - test to be used for the plot title
% OUTPUT
%   plot of grid with the selected line (eg centreline in in 2D or thalweg 
%   in 3D) superimposed on a contour plot or mesh of the grid).
%   hfig - handge to figure
% NOTES
%   if cline includes z values a 3D plot is generated
% SEE ALSO
%   used in gd_xy2sn.m and FGDinterface.addValleyBase in ChannelForm
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2022
%--------------------------------------------------------------------------
%
    hfig = figure('Name','line_on_grid','Units','normalized','Tag','PlotFig');
    ax = axes(hfig);
    if isfield(cline,'z') && ~isempty(cline.z)     %z values defined - 3D       
        channel_z = interp2(grid.x,grid.y,grid.z',cline.x,cline.y,'linear');
        plot3(ax,cline.x,cline.y,channel_z,'-.b','LineWidth',2);
        hold on
        plot3(ax,cline.x,cline.y,cline.z,'-.r','LineWidth',2);
        contour(ax,grid.x,grid.y,grid.z');
        hold off
    else                                           %no z values - 2D
        contourf(ax,grid.x,grid.y,grid.z');
        hold on
        plot(ax,cline.x,cline.y,'-.r','LineWidth',2);
        hold off
    end
    gd_ax_dir(ax,grid.x,grid.y);
    title(plotxt)
    xlabel('x-axis (m)')
    ylabel('y-axis (m)')
    legend('grid','centre line')
end