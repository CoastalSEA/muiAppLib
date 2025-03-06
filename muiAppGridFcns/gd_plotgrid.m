function ax = gd_plotgrid(hfig,grid)
%
%-------function help------------------------------------------------------
% NAME
%   gd_plotgrid.m
% PURPOSE
%   create plot of gridded surface
% USAGE
%   ax = gd_plotgrid(hfig,grid)
% INPUTS
%   hfig - figure handle
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
% OUTPUTS
%   ax - handle to axes object
% SEE ALSO
%   called in GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    hp = findobj(hfig,'Type','axes');
    delete(hp);
    ax = axes(hfig);
    hpc = pcolor(ax,grid.x,grid.y,grid.z');
    hpc.Tag = 'PlotGrid';
    ax = gd_ax_dir(ax,grid.x,grid.y);
    shading interp
    axis equal tight
    xlabel('X-dimension (m)'); 
    ylabel('Y-dimension (m)') 
    
    mnmx = minmax(grid.z);
    if abs(diff(mnmx))>0 
        %only call if there is some variation in the grid
        gd_colormap([mnmx(1),mnmx(2)]);
    end
    cb = colorbar;
    cb.Label.String = 'Elevation (mAD)';  
end 