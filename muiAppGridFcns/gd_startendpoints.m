function points = gd_startendpoints(grid,isdel)
%npts,
%-------function help------------------------------------------------------
% NAME
%   gd_startendpoints.m
% PURPOSE
%   Accept figure to interactively select start and end points on a grid
% USAGE
%   points = gd_startendpoints(grid,isdel)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%  points - struct with x and y fields defining selected start and end points
% SEE ALSO
%   REPLACED BY gd_selectpoints.m which is more generic and has improved
%   interactive editing options.            
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    if nargin<3, isdel = false; end
    
    figtitle = sprintf('Select points');
    promptxt = 'Select path end points?';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Accept','Edit','Quit'};
    position = [0.2,0.4,0.4,0.4];
    [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames,position);
    ax = gd_plotgrid(h_plt,grid);
    hold on
    points.st = [grid.x(2), grid.y(end)/2];
    points.nd = [grid.x(end-1), grid.y(end)/2];
    hpst = plot(ax,points.st(1),points.st(2),'ok','MarkerSize',8,'MarkerFaceColor','w');
    hpnd = plot(ax,points.nd(1),points.nd(2),'xk','MarkerSize',10);
    hold off
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            points = []; return;             
        elseif strcmp(h_but.Tag,'Edit') 
            points = getCoords(points);
            h_but.Tag = 'reset';
            hold on
            hpst.XData = points.st(1); hpst.YData = points.st(2); 
            hpnd.XData = points.nd(1); hpnd.YData = points.nd(2);
            hold off
        elseif strcmp(h_but.Tag,'Quit') 
            points= [];
            ok = 1; delete(h_but);   %keep figure but delete buttons
        else   %user accepted
            ok = 1; delete(h_but);   %keep figure but delete buttons
        end        
    end

    %delete figure if isdel has been set by call.
    if isdel
        delete(h_plt.Parent)
    end
end
%%
function points = getCoords(points)
    %Prompt user to define co-ordinates for start-end points
    promptxt = {'Start co-ordinates (x,y)','End co-ordinates (x,y)'};
    cst = num2str(points.st);
    cnd = num2str(points.nd);
    defaults = {cst,cnd};
    title = 'Define co-ordinates';
    answer = inputdlg(promptxt,title,1,defaults);
    if ~isempty(answer)
        points.st = str2num(answer{1}); %#ok<ST2NM> %vector coordinates
        points.nd = str2num(answer{2}); %#ok<ST2NM>
    end
end