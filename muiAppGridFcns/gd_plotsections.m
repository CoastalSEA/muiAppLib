function gd_plotsections(grid)
%
%-------function help------------------------------------------------------
% NAME
%   gd_plotsections.m
% PURPOSE
%   display grid and allow user to interactively define start and
%   end points of a section line to be plotted in a figure.   
% USAGE
%   gd_plotsections(grid)
% INPUTS
%   grid - struct of x,y,z values that define grid 
% OUTPUT
%   one or more figures of the user defined sections
% NOTES

% SEE ALSO
%   used in GDinterface.plotSections.
%   
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%
    figtitle = sprintf('Plot Sections');
    promptxt = 'Define section, Use Edit to modify and New/Add to plot';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Define','Edit','New','Add','Clear','Quit'};
    position = [0.2,0.4,0.4,0.4];
    [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames,position);
    h_plt.UserData = false;
    ax = gd_plotgrid(h_plt,grid);
    msg = 'Define end points of section first';
    points = [];
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            return; 
        elseif  strcmp(h_but.Tag,'Define')
            points = setLine(ax);
        elseif strcmp(h_but.Tag,'Edit') 
            if isempty(points)
                warndlg(msg); continue
            end
            points = editCoords(ax,points);
        elseif strcmp(h_but.Tag,'New')
            if isempty(points)
                warndlg(msg); continue
            end
            hf_section = getfigure(grid,points);
            hf_section.UserData = h_plt;
            figure(h_plt.Parent); %make the accept figure the current figure
        elseif strcmp(h_but.Tag,'Add')
            %get handle to all exising section plots
            hfs = findobj('Type','figure','Tag','Sections');
            if isempty(points)
                warndlg(msg); continue
            elseif isempty(hfs)
                %no figure yet created, so create
                getfigure(grid,points)
            else
                %prompt user to select figure to use
                prmptxt = 'Select figure to Add section';
                hd = setdialog(prmptxt);
                 
                waitfor(h_plt,'UserData')
                h_plt.UserData = false;
                delete(hd)
                
                ax_section = gca;                
                hold on
                plot_section(ax_section,points,grid);
                hold off
            end
            figure(h_plt.Parent); %make the accept figure the current figure
        elseif strcmp(h_but.Tag,'Clear')
            hpts = findobj(ax,'Tag','mypoints');
            delete(hpts);
            hlns = findobj(ax,'Tag','mysection');
            delete(hlns);
            htxt = findobj(ax,'Tag','mytext');
            delete(htxt);
        elseif strcmp(h_but.Tag,'Quit') 
            ok = 1; continue; 
        end
        %gd_plotgrid(h_plt,grid);
        h_but.Tag = 'reset';
    end
    delete(h_plt.Parent);  
end
%%
function points = editCoords(ax,points)
    %Prompt user to define co-ordinates for start-end points
    hpts = findobj(ax,'Tag','mypoints');
    delete(hpts);
    promptxt = {'Start co-ordinates (x,y)','End co-ordinates (x,y)'};
    cst = num2str([points(1).x,points(1).y]);
    cnd = num2str([points(2).x,points(2).y]);
    defaults = {cst,cnd};
    title = 'Define co-ordinates';
    answer = inputdlg(promptxt,title,1,defaults);
    if isempty(answer), return; end

    st = str2num(answer{1}); %#ok<ST2NM> %vector coordinates
    nd = str2num(answer{2}); %#ok<ST2NM>
    x = num2cell([st(1),nd(1)]);
    y = num2cell([st(2),nd(2)]);
    [points.x] = x{:};
    [points.y] = y{:};
    hold on
    plot(ax,[points(:).x],[points(:).y],'ok','MarkerSize',4,'MarkerFaceColor','w','Tag','mypoints');
    plot(ax,[points(:).x],[points(:).y],'-r','Tag','mysection')
    hold off
end
%%
function points = setLine(ax)
    %define start and end of section line and plot it
    hpts = findobj(ax,'Tag','mypoints');
    if ~isempty(hpts), delete(hpts(1)); end
    nline = length(findobj(ax,'Tag','mysection'));
    prompt1 = 'Select start point';
    points(1) = gd_setpoint(ax,prompt1,false);
    prompt2 = 'Select end point';
    points(2) = gd_setpoint(ax,prompt2,false);
    hpts = findobj(ax,'Tag','mypoints');
    
    hold on
    plot(ax,[points(:).x],[points(:).y],'-r','Tag','mysection')
    hpts(2).MarkerSize = 7;  %make start point marker bigger
    text(ax,points(1).x,points(1).y,sprintf('%d',nline+1),...
        'HorizontalAlignment','center','FontSize',6,'Tag','mytext');
    hold off
    title(ax,'')
end
%%
function plot_section(ax,points,grid)
    %plot the elevations for the defined section
    %plot(ax,[points(:).x],[points(:).y],'-r','Tag','asection')
    xlen = diff([points(:).x]);
    ylen = diff([points(:).y]);
    slen = hypot(xlen,ylen);                %length of section
    gd = gd_dimensions(grid);
    del = (gd.delx+gd.dely)/2;
    spnts = 1:del:slen;                     %points along section
    xq = points(1).x+spnts/slen*xlen;
    yq = points(1).y+spnts/slen*ylen;  

    %check for NaNs in grid
    grid.z(isnan(grid.z)) = gd.zmax;        %clean grid for interpolation
    [X,Y] = meshgrid(grid.x,grid.y);
    zline = interp2(X,Y,grid.z',xq,yq,'makima');
    
    nline = length(findobj(ax,'Tag','asection'));
    lname = sprintf('Section %d',nline+1);
    lines = {'-','--',':','-.'};
    plot(ax,spnts,zline,'LineStyle',lines{rem(nline,4)+1},...
          'LineWidth',1,'Tag','asection','DisplayName',lname)
    legend
end
%%
function getsectionplot(src,~)
    %WindowButtonDownFcn callback function for section figures. when user 
    % clicks on figure. Resets the UserData value in figure panel (h_plt)
    % which is used by waitfor when adding sections to detect when a 
    % figure selection has sbeen made
    %disp(src.Number);
    src.UserData.UserData = true;
end
%%
function hfig = getfigure(grid,points)
    %create a new section figure and plot points
    hfig = figure('Name','Sections','Units','normalized',...
        'WindowButtonDownFcn',@getsectionplot,'Tag','Sections');
    ax_section = axes(hfig);
    plot_section(ax_section,points,grid);
end