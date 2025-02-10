function clines = gd_getcontour(grid,zlevel,isplt)
%
%-------function help------------------------------------------------------
% NAME
%   gd_getcontour.m
% PURPOSE
%   extract a contour at a defined level
% USAGE
%   clines = gd_getcontour(grid,zlevel)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   zlevel - level of contour to be extracted
%   isplt - logical flag true to create plot - optional, default is false
% OUTPUTS
%   clines - struct of x,y vectors defining the contour
% SEE ALSO
%   called in gd_boundary
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
    if nargin<3
        isplt = false;
    end
    [X,Y] = meshgrid(grid.x,grid.y);
    
    hf = figure('Name','Search','Tag','PlotFig','Visible','off');
    
    if isnan(zlevel)
        Z = isnan(grid.z');
        if all(~Z,'all')
            warndlg('No NaNs in grid to use for NaN mask')
            clines = []; return;
        end
        zlevel = 0.5;   %not sure why this value works?
    else
        Z = grid.z';
    end

    C = contour(X, Y, Z, [zlevel,zlevel], 'LineColor', 'k', 'LineWidth', 2);
    if isempty(C), clines = []; return; end

    clines = cleanContours(C,zlevel);
    
     delete(hf)  %delete hidden figure
    
    if isplt
        getPlot(grid,clines);    
    end

end
%%
function ax = getPlot(grid,lines)
    hf = figure('Name','Search','Tag','PlotFig','Visible','on');
    ax = gd_plotgrid(hf,grid);
    hold on 
    plot(ax,lines.x,lines.y,'r','LineWidth', 2);
    hold off
end
%%
function lines = cleanContours(C,level)
    %extract line segments separate with NaN values. return as points struct
    idx = find(C(1,2:end)==level);      %should always return first value
    lines.x = C(1,2:end);
    lines.x(idx) = NaN;
    lines.y = C(2,2:end);
    lines.x(idx) = NaN;

    lines.x(end+1) = NaN;               %add trailing NaN
    lines.y(end+1) = NaN;
    lines = strucfun(@transpose,lines); %output as column vectors
end




























