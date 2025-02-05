function clines = gd_getcontour(grid,zlevel,isplt)
%
%-------function help------------------------------------------------------
% NAME
%   gd_getcontour.m
% PURPOSE
%   
% USAGE
%   points = gd_getcontour(grid,zlevel)
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
    idx = find(C(1,2:end)==level);   %should always return first value
    lines.x = C(1,2:end);
    lines.x(idx) = NaN;
    lines.y = C(2,2:end);
    lines.x(idx) = NaN;

    nrec = length(lines.x);
    idN = find(isnan(lines.x));
    idN = [1,idN,nrec]; 
    xlines = [];  ylines = [];
    %find each line and smooth
    for i=1:length(idN)-1
        xline = lines.x(idN(i):idN(i+1));
        yline = lines.y(idN(i):idN(i+1));
        %first line has no leading nan
        if i>1, xline = xline(2:end); yline = yline(2:end); end
        %last line has not trailing nan
        if idN(i)==nrec, lend = nrec; else, lend = length(xline)-1; end

        if all(xline==0) || all(yline==0)
            %nothing to add
        else    
            xlines = [xlines,xline(1:lend),NaN]; %#ok<AGROW>
            ylines = [ylines,yline(1:lend),NaN]; %#ok<AGROW>
        end
    end
    lines.x = xlines(1:end-1);   %remove trainling NaN
    lines.y = ylines(1:end-1);
end




























