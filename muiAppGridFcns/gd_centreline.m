function cline= gd_centreline(grid,mobj,props,clines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_centreline.m
% PURPOSE
%   create a centreline of a channel using function a_star to trace the
%   shortest path between start and end points whilst finding the deepest
%   points (ie a thalweg).
% USAGE
%   cline = gd_centreline(grid,mobj,props,clines);
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   mobj - mui model instance
%   props - struct for maximum water level and depth exponent (maxwl and dexp)
%   clines - any previously defined lines (optional)
% OUTPUTS
%   cline - struct of x,y vectors defining the centre-line
% NOTEs
% Finds indices of the grid used to find deepest points in channel and
% resolution depends on the xy spacing of the grid used.
% SEE ALSO
%   called in GD_sections; uses gd_selectpoints and a_star
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
    if nargin<4
        clines = [];
    end    

    [X,Y] = meshgrid(grid.x,grid.y);
    N = numel(X);
    xy = [reshape(X,[N,1]),reshape(Y,[N,1])];
    Z = grid.z';

    %accessible map (water) and use -Z as the cost map
    water = true(size(Z));
    water(isnan(Z) | Z>props.maxwl) = false;
    paneltxt = 'Define end points for a channel centre-line. Close window to Quit';
    promptxt3 = {'Select start of path','Select end of path'};
    gridmasked = grid;        gridmasked.z(~water') = NaN;
    points = gd_selectpoints(gridmasked,paneltxt,promptxt3,clines,2,0,true); %2 points, output as struct array, delete figure
    if isempty(points), cline = []; return; end
    
    %index of nearest grid point to selected start end end points    
    start = dsearchn(xy,[points(1).x,points(1).y]); 
    goal = dsearchn(xy,[points(2).x,points(2).y]);
    
    hwb = progressbar(mobj,[],'Computing centre-line');
    %find the shortest path taking account of the cost (depths)
    %Z(Z>maxwl) = 0;
    costfnc = @(z) -(min(z,[],'all')-z).^props.dexp; %weighted inverse depth to favour staying deep
    thalweg = a_star(water, costfnc(Z), start, goal);
    [idy,idx] = ind2sub(size(Z),thalweg); %convert indices to x,y subscripts
    progressbar(mobj,hwb);

    cline.x = [flipud(grid.x(idx));NaN]; %return points in order from start point
    cline.y = [flipud(grid.y(idy));NaN]; %as column vectors terminated with NaNs
end    