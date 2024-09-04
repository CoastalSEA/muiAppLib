function sngrid = gd_sn2xy(grid,cline,isplot)
%
%-------function help------------------------------------------------------
% NAME
%   gd_sn2xy.m
% PURPOSE
%   map grid from curvilinear to cartesian coordinates to return 
%   elevations on the source cartesian grid
% USAGE
%   sngrid = gd_sn2xy(grid,cline,isplot)
% INPUTS
%   grid - struct of the source grid provided by getGrid in GDinterface 
%   cline - x,y struct of centre-line path to use for transformation
%   isplot - logical flag, true to plot output
% OUTPUT
%   sngrid - cartesian grid with elevations mapped from a curvilinear grid
%   defined by the centre-line path. 
% NOTES
%   If a straight channel has been mapped to a meander form using gd_xy2sn, 
%   this function reverses the process to return the original straight channel
%   Uses xy2sn and sn2xy by Bart Vermeulen (2022). Cartesian to Curvilinear 
%   coordinate forward and backward transformation. 
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation 
% SEE ALSO
%   used in GDinterface and the ChannelForm App
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2022
%--------------------------------------------------------------------------
%
    sngrid = grid;
    [X,Y] = meshgrid(grid.x,grid.y);

    N = numel(X);
    xi = reshape(X,[N,1]); yi = reshape(Y,[N,1]);     
    zi= grid.z';            %invert to match x-y dimensions  
    zi = reshape(zi,[N,1]);   
    %xy2sn output S,N relative to the origin (first point) of the channel 
    %centre line. Therefore assume that sn2xy requires the Sx,Ny input to
    %be realtive to the same [0,0] origin. Transform xi,yi using the
    %centrline origin 0xSN0,ySN0] before calling sn2xy.
    xSN0 = cline.x(1);   %origin of centreline used in xy2sn
    ySN0 = cline.y(1);
    xsgn = sign(grid.x(2)-grid.x(1));  %direction of axes (-ve descending)
    ysgn = sign(grid.y(2)-grid.y(1));    
    Sx = xSN0+xsgn*xi;
    Ny = ySN0-ysgn*yi;  
    [xc, yc] = sn2xy(cline.x,cline.y,Sx,Ny);
    %map elevations from the curvilinear grid back to the cartesian grid
    F = scatteredInterpolant(xc,yc,zi,'linear','nearest');
    newz = F(X,Y);
    sngrid.z = newz';  %restore to original dimensions
    
    if isplot
        plotgrids(xc,yc,zi,grid,newz,cline)
    end
end
%%
function plotgrids(xc,yc,zi,grid,newz,cline)
    %generate check plots using (i) curvilinear points  and (ii) cartesian grid
    %with transformed surface
    hf = figure('Name','sn2xy','Units','normalized','Tag','StatFig');

    ax = gd_plotgrid(hf,grid);
    subplot(3,1,1,ax);
    title(sprintf('Source grid: %s',grid.desc))
    
    ax2 = subplot(3,1,2);
    scatter(ax2,xc,yc,[],zi,'f');
    gd_ax_dir(ax2,grid.x,grid.y);
    title('Curvilinear surface converted to Cartesian coordinates')
    xlabel('x-axis (m)')
    ylabel('y-axis (m)')  
    
    ax3 = subplot(3,1,3);
    contourf(ax3,grid.x,grid.y,newz);
    hold on
    plot(ax3,cline.x,cline.y,'-.r');  %**was -cline.y *** needs checking with form model
    hold off
    gd_ax_dir(ax3,grid.x,grid.y);
    title('Cartesian grid of interpolated surface')
    xlabel('x-axis (m)')
    ylabel('y-axis (m)')
end