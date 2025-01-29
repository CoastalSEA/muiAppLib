function points = gd_editlines(grid,paneltxt,promptxt,nlines,isdel)
%npts,
%-------function help------------------------------------------------------
% NAME
%   gd_editlines.m
% PURPOSE
%   Accept figure to interactively edit a line
% USAGE
%   points = gd_editlines(grid,paneltxt,promptxt,nlines,isdel);
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   paneltxt- character string used for title of figure
%   promptxt - cell array of prompts to be used for each line
%   nlines - struct or table of x,y vectors to be edited
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - array of structs with x, y and z fields defining selected points,
% NOTES
% 
% SEE ALSO
%   called in
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
    if nargin<5, isdel = false; end
    isxyz = false;                      %assume just lines
    figtitle = sprintf('Edit lines');
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Add','Edit','Insert','Delete','View','Save'};
    tooltips = {'Add point to set',...
                'Edit a point from the set',...
                'Insert one or more point between two existing points of a line',...
                'Delete a point from the set',...
                'Toggle display of connecting lines on and off',...
                'Save digitised points and exit. Close figure window to Quit without saving'};
    position = [0.3,0.4,0.35,0.5];
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position,tooltips);
    ax = gd_plotgrid(h_plt,grid);
    axis equal  %assume geographical projection or grid of similar dimensions
    axis tight
    %get user to define the required points
    [points,outype] = gd_vec2pnt(nlines);
    ax = plotPoints(ax,points,isxyz);                         
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            points = [];
            return;
        elseif strcmp(h_but.Tag,'Add') 
            promptxt = 'Left click to create points, right click to finish';
            newpnt = gd_setpoints(ax,promptxt,isxyz);         
            points = [points,newpnt];                         %#ok<AGROW> 

        elseif strcmp(h_but.Tag,'Edit') 
            promptxt = 'Select point to edit';
            delpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(delpnt)
                promptxt = 'Left click to create points, right click on any point to quit';
                newpnt = gd_setpoint(ax,promptxt,isxyz);     
                points = deletepoint(ax,points,delpnt,newpnt);
            end

        elseif strcmp(h_but.Tag,'Insert')
            pnt = startendpoints(ax,points); %define insertion position
            if ~isempty(pnt)
                promptxt = 'Left click to create points, right click to finish';
                newpnts = gd_setpoints(ax,promptxt,isxyz);    %points to insert
                if ~isempty(newpnts)
                    points = insertpoints(ax,points,pnt,newpnts); %add points
                end
            end

        elseif strcmp(h_but.Tag,'Delete') 
            promptxt = 'Select point to Delete, right click on any point to quit';
            delpnt = gd_getpoint(ax,promptxt);   %get point to delete
            if ~isempty(delpnt)
                points = deletepoint(ax,points,delpnt,[]);  %delete the point                
            end

        elseif strcmp(h_but.Tag,'View') 
            ax = toggle_view(ax,points);

        else   %user accepted
            ok = 1;                                  %accepted so end loop

        end   
        h_but = resetbutton(ax,h_but);
        resetpoints(ax);
    end

    %convert format of output if required
    points = gd_pnt2vec(points,outype);
    
    %delete figure if isdel has been set by call.
    if isdel
        delete(h_plt.Parent)
    end
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
function ax = plotPoints(ax,points,isxyz)
    %plot the imported lines
    hold on
    for i=1:length(points)
        H = plot(ax,points(i).x,points(i).y,'ok','MarkerSize',4,...
                                   'MarkerFaceColor','w','Tag','mypoints');
        H.ButtonDownFcn = {@LineSelected, H};
        H.UserData = int32(0);
        if isxyz
            text(ax,points(i).x,points(i).y,sprintf('  %.1f',point.z),...
                               'Color','white','FontSize',6,'Tag','ztext');
        end
    end
    hold off
    %nested function
        function LineSelected(src, evt, H)
            if evt.Button==1
                H(H==src).Color = 'r';
            elseif evt.Button==3
                H(H==src).Color = 'k';        
            end
            H(H==src).UserData = evt.Button;
        end
end
%%
function pnt = startendpoints(ax,points)
    %find two adjacent points and return in pnt struct: start and adjacent
    if isempty(points), pnt = []; warndlg('No points defined'); return; end
    
    %plot line(s) connecting points
    hold on
    hline = plot(ax,[points(:).x],[points(:).y],'.--k','PickableParts','none',...
                                 'MarkerEdgeColor','r','Tag','insertline');                                                   
    hold off
    
    %get the start and adjacent points for the insertion of new points
    isok = false;
    while ~isok
        promptxt = 'Select first point';
        pnt.start = gd_getpoint(ax,promptxt);   
        if isempty(pnt)                        %user right clicked a point
            pnt = []; delete(hline); 
            return; 
        else
            h_pnts = findobj(ax,'Tag','mypoints'); 
            idx = [h_pnts.UserData]>0;
            h_pnts(idx).UserData = int32(0);   %reset the selected point
        end 
        promptxt = 'Select adjacent point';
        pnt.adjacent = gd_getpoint(ax,promptxt);
        if isempty(pnt)                        %user right clicked a point
            pnt = []; delete(hline); 
            return;  
        elseif ~isequal(pnt.start.x,pnt.adjacent.x) || ~isequal(pnt.start.y,pnt.adjacent.y)
            %the selected points are different
            isok = true;
        end
    end
    delete(hline)
end

%%
function points = insertpoints(ax,points,pnt,newpnts)
    %insert additional points, 'newpnts', between the selected points,
    %'pnt', in the digitised points vector,'points'.
    idpos(1) = find([points(:).x]==pnt.start.x & [points(:).y]==pnt.start.y); 
    idpos(2) = find([points(:).x]==pnt.adjacent.x & [points(:).y]==pnt.adjacent.y); 
    if idpos(2)<idpos(1), idpos = fliplr(idpos); end  %check that in ascending order
    %insert the new points in the correct interval for the line
    %to be in the right order they need to be created in the same direction
    %as the initial line is being created - this is NOT checked.
    points = [points(1:idpos(1)),newpnts,points(idpos(2):end)];

    %restore colour of selected points
    h_pnts = findobj(ax,'Tag','mypoints');
    idpos(1) = find([h_pnts(:).XData]==pnt.start.x & [h_pnts(:).YData]==pnt.start.y,1,'first'); 
    idpos(2) = find([h_pnts(:).XData]==pnt.adjacent.x & [h_pnts(:).YData]==pnt.adjacent.y,1,'first'); 
    h_pnts(idpos(1)).Color = zeros(1,3);
    h_pnts(idpos(2)).Color = zeros(1,3);
end
%%
function points = deletepoint(ax,points,delpoint,newpnt)
    %delete point defined by delpnt if newpnt is empty, otherwise edit
    %point to new value as defined in newpnt
    h_pnts = findobj(ax,'Tag','mypoints');
    idx = [h_pnts(:).XData]==delpoint.x & [h_pnts(:).YData]==delpoint.y; 
    idp = [points(:).x]==delpoint.x & [points(:).y]==delpoint.y; 

    if isempty(newpnt)              %call to delete point (newpnt is empty)
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
    end
end

%%
function ax = toggle_view(ax,points)
    %switch a line of the lines and points defined with the start point
    %of each line emphasised with a red circle marker.
    hline = findobj(ax,'Tag','viewline');
    if isempty(hline)           %toggle line and points on
        hold on
        plot(ax,[points(:).x],[points(:).y],'.--k','PickableParts','none',...
                                   'MarkerEdgeColor','r','Tag','viewline'); 
        idx = [1,find(isnan([points(:).x]))+1];
        plot(ax,[points(idx).x],[points(idx).y],'ob','PickableParts','none',...
                  'MarkerSize',6,'MarkerEdgeColor','r', 'Tag','viewline'); 
        hold off  
    else                        %toggle line and points off
        delete(hline)
    end
end