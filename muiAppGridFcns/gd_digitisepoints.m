function points = gd_digitisepoints(grid,paneltxt,outype,isxyz,isdel)
%
%-------function help------------------------------------------------------
% NAMEpnts = 
%   gd_digitisepoints.m
% PURPOSE
%   Accept figure to interactively digitise x,y,z points on a grid and edit
%   z elevations, if required
% USAGE
%   points = gd_digitisepoints(grid,paneltxt,outype,isxyz,isdel)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   paneltxt - character string used for title of figure
%   outype - format of output - see Outputs for details
%   isxyz - logical flag true to input z values - optional, default is false
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - outype=0: array of structs with x, y and z fields defining points,
%            outype=1: Nx2 or Nx3 array.
%            outype=2: struct with x, y (and z) vector fields
%            outype=3: table with x, y (and z) vector fields
%            points = [] if user closes figure, or no points defined
% NOTES
%   Each new line is separated by NaN values in the xyz vectors. When using
%   the View button the start of each line is indicated by a red circle
%   marker.
%   Functions is similar to gd_editlines which handles x,y points, whearas
%   gd_digitisepoints handles x,y or x,y,z points
% SEE ALSO
%   called in GDinterface.getGridLine 
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
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
    butnames = {'New','Add','Edit','Insert','Join/Split','Delete','View','Use'};
    tooltips = {'Start a new line of points',...
                'Add points to a selected line',...
                'Edit a point in any line',...
                'Insert one or more point between two existing points of a line',...
                'Join two lines or split existing line (after point selected)',...
                'Delete a point from a line',...
                'Toggle display of connecting lines on and off',...
                'Use digitised points and exit. Close figure window to Quit without saving points/lines'};
    % position = [0.3,0.4,0.35,0.5];
    position = [0,0,1,1];
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position,tooltips);
    ax = gd_plotgrid(h_plt,grid);
    axis equal  %assume geographical projection or grid of similar dimensions
    axis tight
    %get user to define the required points
    points = [];
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            points = [];
            return;

        elseif strcmp(h_but.Tag,'New') 
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
            promptxt = 'Select point at end of line to extend';
            endpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(endpnt)
                promptxt = sprintf('Left click to create points, right click to finish\nFirst new point should be closest to selected end point');
                newpnts = gd_setpoints(ax,promptxt,isxyz);  
                points = addpoints(points,newpnts,endpnt,true);     
            end

        elseif strcmp(h_but.Tag,'Edit') 
            promptxt = 'Select point to edit';
            delpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(delpnt)
                promptxt = 'Left click to create points, right click on any point to quit';
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

        elseif strcmp(h_but.Tag,'Join/Split')
            promptxt = 'Select point within a line to Split line and end point to Join to another line';
            initpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(initpnt)
                location = findPointonLine(points,initpnt);
                if location<5
                    resetpoints(ax,points,isxyz); %reset the selected point
                    promptxt = 'Select point to Join to line';
                    joinpnt =gd_getpoint(ax,promptxt); 
                    if ~isempty(joinpnt)
                        points = joinlines(points,initpnt,joinpnt,isxyz);
                    end
                else
                    points = splitlines(points,initpnt,isxyz);
                end
                delete(findobj(ax,'Tag','viewline'));
                ax = toggle_view(ax,points);
            end

        elseif strcmp(h_but.Tag,'Delete') 
            promptxt = 'Select point to Delete, right click on any point to quit';
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
function newpnts = resetpoints(ax,points,isxyz)
    %gd_getpoint sets point UserData and color when clicked on. Reset in 
    %case user clicked on points without making an action selection
    newpnts = [];    %reset newpnts in case user Quits when adding, etc
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
end

%%
function points = addpoints(points,newpnts,endpnt,isadd)
    %find end to add points to and extend existing line
    idp = find([points(:).x]==endpnt.x & [points(:).y]==endpnt.y); 
    if isadd
        %check direction of new points relative to existing line
        xlen = [endpnt.x-newpnts(1).x,endpnt.x-newpnts(end).x];
        ylen = [endpnt.y-newpnts(1).y,endpnt.y-newpnts(end).y];
        [~,ide] = min(hypot(xlen,ylen));                %index of nearest point
        if ide>1
            newpnts = fliplr(newpnts);
        end
    end

    %find location and add new points
    location = findPointonLine(points,endpnt);
    switch location
        case 1  %first point in vector
            points = [fliplr(newpnts),points];
        case 2  %last point in vector
            points = [points,newpnts];
        case 3  %first point of line (after a NaN)
            points = [points(1:idp-1),fliplr(newpnts),points(idp:end)];
        case 4  %last point of line (before a NaN)
            points = [points(1:idp),newpnts,points(idp+1:end)];
    end
end

%%
function location = findPointonLine(points,endpnt)
    %find the location of a point in a vector of lines
    % location is an index:
    % empty = point not found
    % 1 - first point
    % 2 - last point
    % 3 - first point of line in vector (after a NaN)
    % 4 - last point of line in vector (before a NaN)
    % 5 - point within line 
    nrec = length(points);
    idp = find([points(:).x]==endpnt.x & [points(:).y]==endpnt.y); 

    if sum(idp)==0
        location = []; return;
    elseif idp==1
        location = 1;  %first point in vector
    elseif idp==nrec
        location = 2;  %last point in vector
    elseif isnan(points(idp-1).x)
        location = 3;  %first point of line (after a NaN)
    elseif isnan(points(idp+1).x)
        location = 4;  %last point of line (before a NaN)
    else
        location = 5;  %point found lies on line (ie not an end point of a line segment)
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
function points = joinlines(points,initpnt,joinpnt,isxyz)
    %join two lines based on selected end points
    nrec = length(points);
    msg = 'Cannot join line to itself';
    newpnts.x = NaN; newpnts.y = NaN;        %line termination
    if isxyz, newpnts.z = NaN; end
    idi = find([points(:).x]==initpnt.x & [points(:).y]==initpnt.y); 
    idj = find([points(:).x]==joinpnt.x & [points(:).y]==joinpnt.y); 
    
    %check not trying to join a line to itself
    if idi==idj, warndlg(msg); return; end
    Lineidi = gd_findline(points,initpnt);
    Lineidj = gd_findline(points,joinpnt);
    if Lineidi==Lineidj, warndlg(msg); return; end

    idN = find(isnan([points(:).x]));
    idN = [1,idN,nrec];   
    %extract lines, concatenate and add to end of points array
    [line1,idL1] = getIndex(points,idN,Lineidi);
    [line2,idL2] = getIndex(points,idN,Lineidj);
    idx = [idL1,idL2];
    points(idx) = [];
    if ~isempty(points) && ~isnan(points(end).x)        
        points = [points,newpnts];  %add NaNs to end of existing lines
    end
    %add the two lines. NB this sorts the lines to minimise the distance
    %from the point selected on line1 to the start or end of line 2
    isend = line2(end).x==joinpnt.x & line2(end).y==joinpnt.y;
    if isend, line2 = fliplr(line2); end %line direction based on selected join point       
    newline = addpoints(line1,line2,initpnt,false);

    points = [points,newline];

    %nested function-------------------------------------------------------
    function [aline,idL] = getIndex(points,idN,line)
        %get the indices of the line to be extracted
        % idL - indices of line including trailing NaNs
        idL =idN(line):idN(line+1);
        aline = points(idL);
        if isnan(aline(1).x)
            idL = idN(line)+1:idN(line+1);
            aline(1) = [];
        end
        %
        if ~isnan(aline(end).x)
            aline = [aline,newpnts];
        end
        aline = aline(1:end-1);
    end %------------------------------------------------------------------
end

%%
function points = splitlines(points,initpnt,isxyz)
    %split a line at the defined point
    idi = find([points(:).x]==initpnt.x & [points(:).y]==initpnt.y); 
    newpnts.x = NaN; newpnts.y = NaN;        %line termination
    if isxyz, newpnts.z = NaN; end
    points = [points(1:idi),newpnts,points(idi+1:end)];
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
