function [sngrid,cline] = gd_xy2sn(grid,cline,isxy,isplot)
%
%-------function help------------------------------------------------------
% NAME
%   gd_xy2sn.m
% PURPOSE
%   map grid from cartesian to curvilinear coordinates with option to return 
%   elevations on the source cartesian grid, or as a curvilinear grid
% USAGE
%   sngrid = gd_xy2sn(grid,cline,isxy,isplot)
% INPUTS
%   grid - struct of the source grid provided by getGrid in GDinterface 
%   cline - x,y struct of centre-line path to use for transformation
%   isxy - logical flag, true if the curvilinear z values are to be mapped 
%          onto the source cartesian grid, false to return the curvilinear
%          grid coordinates
%   isplot - logical flag, true to plot output
% OUTPUT
%   sngrid - cartesian grid with elevation mapped from curvilinear grid, or
%            curvilinear grid
%   tol - struct containing 'maxiter' for the maximum number of iterations
%         and 'tolerance' for the tolerance used in the call to xy2sn
% NOTES
%   uses xy2sn and sn2xy by Bart Vermeulen (2022). Cartesian to Curvilinear 
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

    Num = numel(X);
    xi = reshape(X,[Num,1]); yi = reshape(Y,[Num,1]);  
    zi= grid.z';            %invert to match x-y dimensions  
    zi = reshape(zi,[Num,1]);   
    
    clx = cline.x; cly = cline.y;
    if isempty(cline.maxiter) || isempty(cline.tolerance)
        %plot selected line on grid and give user option to quit
        hf = gd_lineongrid_plot(grid,cline,'Selected centre line');
        answer = questdlg('Accept centreline and continue?','Curvilinear','Accept','Quit','Accept');
        if strcmp(answer,'Quit')
            sngrid = []; delete(hf); return; 
        end

        promptxt = {'Maximum number of iterations','Tolerance (m)'};
        tol = num2str(1e-3*nanmean(sqrt(sum(diff(clx).^2))));
        vals = inputdlg(promptxt,'xy2sn',1,{'500',tol});
        if isempty(vals), sngrid = []; return; end
        cline.maxiter = str2double(vals{1});
        cline.tolerance = str2double(vals{2});
    end
    
    [S, N] = xy2sn(clx,cly,xi,yi,'MaximumIterations',cline.maxiter,...
              'Tolerance',cline.tolerance);
    %S and N are returned from origin of [0,0]. Adjust to grid origin
    xSN0 = cline.x(1);   %origin of centreline used in xy2sn
    ySN0 = cline.y(1);
    xsgn = sign(grid.x(2)-grid.x(1));  %direction of axes (-ve descending)
    ysgn = sign(grid.y(2)-grid.y(1));    
    S = xSN0+xsgn*S;
    N = ySN0+ysgn*N;   
    if isxy
        %map elevations from the curvilinear grid back to the cartesian grid
        F = scatteredInterpolant(S,N,zi,'linear','nearest');
        [X,Y] = meshgrid(grid.x,grid.y);   
        newz = F(X,Y);
        sngrid.z = newz';  %restore to original dimensions
    else
        %return as curvilinear grid tuples
        sngrid.x = S;
        sngrid.y = N;
        sngrid.z = zi;
    end
    
    if isplot
        plotgrids(S,N,zi,grid,newz,cline)
    end
end
%%
function plotgrids(S,N,zi,grid,newz,cline)
    %generate check plots using (i) curvilinear points  and (ii) cartesian grid
    %with transformed surface
    hf = figure('Name','xy2sn','Units','normalized','Tag','StatFig');
    
    ax = gd_plotgrid(hf,grid);
    subplot(3,1,1,ax);
    title(sprintf('Source grid: %s',grid.desc))
    
    ax2 = subplot(3,1,2);
    scatter(ax2,S,N,[],zi,'f');
    gd_ax_dir(ax2,grid.x,grid.y);
    title('Cartesian surface converted to Curvlinear coordinates')
    xlabel('S-axis (m)')
    ylabel('N-axis (m)')  
    
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