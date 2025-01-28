function points = gd_digitisepoints(grid,figtxt,outype,isxyz,isdel)
%
%-------function help------------------------------------------------------
% NAMEpnts = 
%   gd_digitisepoints.m
% PURPOSE
%   Accept figure to interactively digitise x,y points on a grid and add
%   z elevations, if required
% USAGE
%   points = gd_digitisepoints(grid,figtxt,outype,isxyz,isdel)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   figtxt - character string used for title of figure
%   outype - format of output - see Outputs for details
%   isxyz - logical flag true to input z values - optional, default is false
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - outype=0: array of structs with x, y and z fields defining selected points,
%            outype=1: Nx2 or Nx3 array.
%            outype=2: struct with x, y (and z) vector fields
%            points = [] if user closes figure, or no points defined
%            outype=1: 
% NOTES
%   Each new line is separated by NaN values in the xyz vectors. When using
%   the View button the start of each line is indicated by a red circle
%   marker.
% SEE ALSO
%   called in GDinterface.getGridLine, similar to gd_selectpoints which
%   only a specified number of x,y points
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    if nargin<5
        isdel = false; 
    elseif nargin<4
        isxyz = false; 
        isdel = false;
    end
 
    figtitle = sprintf('Digitise points');
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'New line','Add','Edit','Insert','Delete','View','Save/Exit'};
    tooltips = {'Start a new line of points',...
                'Add points to the last line created',...
                'Edit a point in any line',...
                'Insert one or more point between two existing points of a line',...
                'Delete a point in any line',...
                'Toggle display of connecting lines on and off',...
                'Save digitised points and exit. Close figure window to Quit without saving'};
    position = [0.2,0.4,0.4,0.4];
    [h_plt,h_but] = acceptfigure(figtitle,figtxt,tag,butnames,position,tooltips);
    ax = gd_plotgrid(h_plt,grid);
    axis equal  %assume geographical projection or grid of similar dimensions
    axis tight
    %get user to define the required points
    points = [];
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but)   %this handles the user deleting figure window 
            points = []; 
            return;    

        elseif strcmp(h_but.Tag,'New line') 
            if ~isempty(points)
                newpnts.x = NaN; newpnts.y = NaN;        %line termination
                if isxyz
                    newpnts.z = NaN;
                end
                points = [points,newpnts];                     %#ok<AGROW> 
            end            
            promptxt = 'Left click to create points, right click to finish';
            newpnts = gd_setpoints(ax,promptxt,isxyz);   %get points to add
            points = [points,newpnts];                         %#ok<AGROW> 
            
        elseif strcmp(h_but.Tag,'Add')             
            promptxt = 'Left click to create points, right click to finish';
            newpnts = gd_setpoints(ax,promptxt,isxyz);   %get points to add
            points = [points,newpnts];                         %#ok<AGROW> 

        elseif strcmp(h_but.Tag,'Edit')
            promptxt = 'Select point to edit';
            delpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(delpnt)
                promptxt = 'Left click to create points, right click to quit';
                newpnt = gd_setpoint(ax,promptxt,isxyz);
                points = deletepoint(ax,points,delpnt,newpnt,isxyz);
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
            promptxt = 'Select point to Delete';
            delpnt = gd_getpoint(ax,promptxt);   %get point to delete
            if ~isempty(delpnt)
                points = deletepoint(ax,points,delpnt,[],isxyz);  %delete the point                
            end

        elseif strcmp(h_but.Tag,'View')
            ax = toggle_view(ax,points);
          
        elseif strcmp(h_but.Tag,'Quit') 
            %no longer used as a button - use close figure
            points = [];
            ok = 1;
            isdel = true;

        else   %user accepted
            h_pnts = findobj(ax,'Tag','mypoints');  %remove construction points
            delete(h_pnts)
            if exist('points','var')                %plot points to be returned
                hold on
                plot(ax,[points(:).x],[points(:).y],'.-k','MarkerEdgeColor','r')                                                   
                hold off
            else                                    %no points to be returned
                points = [];
                isdel = true;
            end
            ok = 1;                                 %accepted so end loop
        end     
         newpnts = resetpoints(ax,points,isxyz);
         h_but = resetbutton(ax,h_but);
    end

    %convert format of output if required
    if outype==1                %return as an array
        points = gd_pnt2vec(points,outype);
    elseif outype==2            %return as a struc of vectors
        points = gd_pnt2vec(points,outype);
    end

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
function newpnts = resetpoints(ax,points,isxyz)
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
    %
    if isxyz
        h_txt = findobj(ax,'Tag','ztext');
        delete(h_txt)
        for i=1:length(points)
            txt= sprintf('  %.1f',points(i).z);
            text(ax,[points(i).x],[points(i).y],txt,...
                'Color','white','FontSize',6,'Tag','ztext');
        end
    end
    newpnts = [];    
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
function points = deletepoint(ax,points,delpoint,newpnt,isxyz)
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
            if isxyz
                h_txt = findobj(ax,'Tag','ztext');
                delete(h_txt(idx))
            end
        end
    else                            %call to edit point (newpnt is new point)
        h_pnts(idx).XData = newpnt.x;
        h_pnts(idx).YData = newpnt.y;        
        points(idp).x = newpnt.x;
        points(idp).y = newpnt.y;
        if isxyz
            points(idp).z = newpnt.z;
            h_txt = findobj(ax,'Tag','ztext');
            delete(h_txt(idp))
        end
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
