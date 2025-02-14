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
%
% SEE ALSO
%   used in GDinterface.plotSections.
%   
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%
    figtitle = sprintf('Plot Sections');
    promptxt = 'Define section, Use Edit to modify and New/Add to plot. Define and plot one section at a time';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Define','Edit','New','Add','Clear'};
    tooltips = {'Define a new section',...
                'Edit the end points of a section',...
                'Generate a plot using the currently defined sections',...
                'Add a section to an exiting plot',...
                'Clear the current section lines. Close figure window to Quit'};
    position = [0.3,0.4,0.35,0.5];
    [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames,position,0.8,tooltips);
    h_plt.UserData = false;
    ax = gd_plotgrid(h_plt,grid);

    msg = 'Use ''Define'' to initialise a section';
    isxyz = false;
    points = [];
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            return; 

        elseif  strcmp(h_but.Tag,'Define')
            points = setLine(ax);

        elseif strcmp(h_but.Tag,'Edit') 
            if isempty(points), warndlg(msg); continue; end %no sections defined
            promptxt = 'Select point to edit';
            delpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(delpnt)
                promptxt = 'Left click to create points, right click to quit';
                newpnt = gd_setpoint(ax,promptxt,isxyz);
                if ~isempty(newpnt)
                    points = deletepoint(ax,points,delpnt,newpnt);
                end               
            end

        elseif strcmp(h_but.Tag,'New')
            if isempty(points), warndlg(msg); continue; end %no sections defined
            hf_section = getfigure(grid,points);
            hf_section.UserData = h_plt;
            figure(h_plt.Parent); %make the accept figure the current figure

        elseif strcmp(h_but.Tag,'Add')
            %get handle to all exising section plots
            hfs = findobj('Type','figure','Tag','SectionsPlot');
            if isempty(points)
                warndlg(msg); continue                 %no sections defined
            elseif isempty(hfs)                
                getfigure(grid,points)        %no figure created, so create
            else
                %prompt user to select figure to use
                prmptxt = 'Select figure to Add section';
                hd = setdialog(prmptxt);                 
                waitfor(h_plt,'UserData')
                h_plt.UserData = false;
                delete(hd)
                
                ax_section = gca;                
                plot_section(ax_section,points,grid);
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
        h_but = resetbutton(ax,h_but);
        resetpoints(ax);
    end
    delete(h_plt.Parent);  
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
    plot(ax,[points(:).x],[points(:).y],'-r','PickableParts','none','Tag','mysection')
    hpts(2).MarkerSize = 7;  %make start point marker bigger
    text(ax,points(1).x,points(1).y,sprintf('%d',nline+1),'Clipping', 'on',...
        'HorizontalAlignment','center','PickableParts','none','FontSize',6,'Tag','mytext');
    hold off
end

%%
function h_but = resetbutton(ax,h_but)
    %reset button and panel title string to get new input
     if ishandle(h_but)   %check that figure has not been deleted
        ax.Title.String = 'Select a button';
        h_but.Tag = 'reset';        
     end     
end

%%
function resetpoints(ax)
    %gd_getpoint sets point UserData and color when clicked on. Reset in 
    %case user clicked on points without making an action selection
    h_pnts = findobj(ax,'Tag','mypoints');
    if isempty(h_pnts), return; end
    idx = [h_pnts.UserData]>0;
    if any(idx)
        [h_pnts(idx).UserData] = repmat(int32(0),sum(idx),1);
        cellobj = {[h_pnts(idx)]};
        [cellobj{:}.Color] = repmat(zeros(1,3),sum(idx),1);
    end 
end

%%
function points = deletepoint(ax,points,delpoint,newpnt)
    %delete point defined by delpnt if newpnt is empty, otherwise edit
    %point to new value as defined in newpnt
    h_pnts = findobj(ax,'Tag','mypoints');
    idx = [h_pnts(:).XData]==delpoint.x & [h_pnts(:).YData]==delpoint.y; 
    idp = [points(:).x]==delpoint.x & [points(:).y]==delpoint.y; 
    
    if ~any(idp)
        warndlg('Not current section'),          
    elseif isempty(newpnt)          %call to delete point (newpnt is empty)
        answer = questdlg('Confirm deletion','Delete point','Yes','No','Yes');
        if strcmp(answer,'Yes')
            delete([h_pnts(idx)]);  %remove any existing points
            points(idp) = [];
        end
    else                            %call to edit point (newpnt is new point)        
        h_pnts(idx).XData = newpnt.x;
        h_pnts(idx).YData = newpnt.y;        
        points(idp).x = newpnt.x;
        points(idp).y = newpnt.y;        
        updateSection(ax,points,delpoint,newpnt);
        idd = [h_pnts(:).XData]==newpnt.x & [h_pnts(:).YData]==newpnt.y;
        delete([h_pnts(idd & ~idx)]);         %remove new point       
    end 
end
%%
function updateSection(ax,points,delpoint,newpnt)
    %update the linework and text associated with sections on main figure
    %assumes only the last section can be edited. for case of all sections
    %available and editable see gd_sectionlines
    h_lns = findobj(ax,'Tag','mysection');
    nline = length(h_lns);
    idl = [h_lns(:).XData]==delpoint.x & [h_lns(:).YData]==delpoint.y; 
    idl =  any(reshape(idl,[],2)); %rows are lines columns xy points
    delete(h_lns(idl))

    h_txt = findobj(ax,'Tag','mytext'); 
    txtpnts = reshape([h_txt(:).Position],3,[])';
    idt = txtpnts(:,1)==delpoint.x & txtpnts(:,2)==delpoint.y; 
    if any(idt)
        delete(h_txt(idt))
    end
    
    hold on
    plot(ax,[points(:).x],[points(:).y],'-r','PickableParts','none','Tag','mysection')
    if any(idt)
        %h_pnts(2).MarkerSize = 7;  %make start point marker bigger
        text(ax,newpnt.x,newpnt.y,sprintf('%d',nline),'Clipping', 'on',...
            'HorizontalAlignment','center','PickableParts','none','FontSize',6,'Tag','mytext');
    end
    hold off    
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
    hold on
    plot(ax,spnts,zline,'LineStyle',lines{rem(nline,4)+1},...
          'LineWidth',1,'Tag','asection','DisplayName',lname)
    hold off
    legend
end
%%
function getsectionplot(src,~)
    %WindowButtonDownFcn callback function for section figures. when user 
    % clicks on figure. Resets the UserData value in figure panel (h_plt)
    % which is used by waitfor when adding sections to detect when a 
    % figure selection has been made
    % disp(src.Number);
    src.UserData.UserData = true;
end
%%
function hfig = getfigure(grid,points)
    %create a new section figure and plot points
    hfig = figure('Name','Sections','Units','normalized',...
        'WindowButtonDownFcn',@getsectionplot,'Tag','SectionsPlot');
    ax_section = axes(hfig);
    plot_section(ax_section,points,grid);
end