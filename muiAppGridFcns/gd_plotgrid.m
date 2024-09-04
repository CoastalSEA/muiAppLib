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
    pcolor(ax,grid.x,grid.y,grid.z')
    ax = gd_ax_dir(ax,grid.x,grid.y);
    shading interp
    colormap('parula')
    colorbar
    xlabel('Length (m)'); 
    ylabel('Width (m)')
end 