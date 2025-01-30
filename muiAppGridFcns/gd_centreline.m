function [cline,hf] = gd_centreline(grid,mobj,props,clines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_centreline.m
% PURPOSE
%   create a centreline of a channel using function a_star to trace the
%   shortest path between start and end points whilst finding the deepest
%   points (ie a thalweg).
% USAGE
%   cline = gd_centreline(grid,mobj,clines);
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   mobj - mui model instance
%   props - struct for maximum water level and depth exponent (maxwl and dexp)
%   clines - any previously defined lines (optional)
% OUTPUTS
%   cline - struct of x,y vectors defining the centre-line
% SEE ALSO
%   called in ??? uses gd_selectpoints and a_star
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
if nargin<4
    clines = [];
end
    %     gridclasses = {'EDBimport'}; %add other classes if needed
    %     promptxt1 = {'Select Case to plot (Cancel to quit):','Select timestep:'};
    %     [obj,~,irec] = selectCaseDatasetRow(mobj.Cases,[],...
    %                                                  gridclasses,promptxt1,1);
    %     if isempty(obj) || isempty(irec), return; end
    % 
        % desc = sprintf('%s at %s',obj.Data.Grid.Description,char(obj.Data.Grid.RowNames(irec)));
    %     grid = getGrid(obj,irec);   %grid for selected year
    [X,Y] = meshgrid(grid.x,grid.y);
    N = numel(X);
    xy = [reshape(X,[N,1]),reshape(Y,[N,1])];
    Z = grid.z';    

%     %get maximum water level to define 
%     promptxt2 = {'Maximum accessible water level?','Depth exponent'};
%     defaults = {num2str(max(Z,[],'all')),'5'};
%     answer = inputdlg(promptxt2,'Water level',1,defaults);
%     if isempty(answer), return; end %user cancelled
%     maxwl = str2double(answer{1});
%     dexp = str2double(answer{2});

    %accessible map (water) and use -Z as the cost map
    water = true(size(Z));
    water(isnan(Z) | Z>props.maxwl) = false;
    paneltxt = 'Define end points for a channel centre-line. Close window to Quit';
    promptxt3 = {'Select start of path','Select end of path'};
    gridmasked = grid;        gridmasked.z(~water') = NaN;
    points = gd_selectpoints(gridmasked,paneltxt,promptxt3,clines,2,0,true); %2 points, struct array, delete figure
    if isempty(points), cline = []; hf = []; return; end
    
    %index of nearest grid point to selected start end end points    
    start = dsearchn(xy,[points(1).x,points(1).y]); 
    goal = dsearchn(xy,[points(2).x,points(2).y]);
    
    hwb = progressbar(mobj,[],'Computing centre-line');
    %find the shortest path taking account of the cost (depths)
    %Z(Z>maxwl) = 0;
    costfnc = @(z) -(min(z,[],'all')-z).^props.dexp; %weighted inverse depth to favour staying deep
    thalweg = a_star(water, costfnc(Z), start, goal);
    [idy,idx] = ind2sub(size(Z),thalweg);
    progressbar(mobj,hwb);

    cline.x = grid.x(idx);
    cline.y = grid.y(idy);
    %plot base map of initial grid selection and defined mask
    hf = figure('Name','Thalwegs','Units','normalized','Tag','PlotFig');  
    hf.Position = [0,0,1,1];
    ax = gd_plotgrid(hf,grid);
    axis equal  %assume geographical projection or grid of similar dimensions
    axis tight
    colormap(ax,'gray');
    lines = {'-','--',':','-.'};
    
    hs = findobj(ax.Children,'Type','surface');
    hs.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on
    hp = plot(ax,[points(:).x],[points(:).y],'ok','MarkerSize',8,'MarkerFaceColor','w','Tag','mypoints');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    plot(ax,cline.x,cline.y,'r','LineStyle',lines{1},'LineWidth',1,...
                                                      'DisplayName',grid.desc);
    hold off
    title('Centre-line between defined start and end points')
    legend
end    