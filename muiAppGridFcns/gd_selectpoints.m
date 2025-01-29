function points = gd_selectpoints(grid,paneltxt,promptxt,outype,npts,isdel)
%npts,
%-------function help------------------------------------------------------
% NAME
%   gd_selectpoints.m
% PURPOSE
%   Accept figure to interactively create a specified number of x,y points
%   on a grid
% USAGE
%   points = gd_selectpoints(grid,paneltxt,promptxt,outype,npts,isdel);
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   paneltxt- character string used for title of figure
%   promptxt - cell array of prompts to be used for each point being defined
%              a single char string cell is apended with a count for each point,
%              whereas a cell array should have a length of npts.
%   outype - format of output - see Outputs for details
%   npts - number of points to be selected
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - outype=0: array of structs with x, y and z fields defining selected points,
%            outype=1: Nx2 or Nx3 array.
%            outype=2: struct with x, y (and z) vector fields
%            points = [] if user closes figure, or no points defined
%            outype=1:
% NOTES
%   Captures x,y for the number of points specified in npts. However, the
%   user can add, edit and/or delete the initial input of npts. Edit and
%   delete retain initial order of points, add appends points to struct
% SEE ALSO
%   similar to gd_digitisepoints, which also captures z values at x,y points
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    if nargin<6, isdel = false; end

    figtitle = sprintf('Select points');
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Add','Edit','Delete','Save'};
    tooltips = {'Add point to set',...
                'Edit a point from the set',...
                'Delete a point from the set',...
                'Save digitised points and exit. Close figure window to Quit without saving'};
    position = [0.3,0.4,0.35,0.5];
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position,tooltips);
    ax = gd_plotgrid(h_plt,grid);
    axis equal  %assume geographical projection or grid of similar dimensions
    axis tight
    %get user to define the required points
    points = setPoints(ax,npts,promptxt);    
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            points = [];
            return;
        elseif strcmp(h_but.Tag,'Add') 
            promptxt = 'Left click to create points, right click on any point to quit';
            newpnt = gd_setpoint(ax,promptxt,false);          %isxyz false 
            points = [points,newpnt];                         %#ok<AGROW> 

        elseif strcmp(h_but.Tag,'Edit') 
            promptxt = 'Select point to edit';
            delpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(delpnt)
                promptxt = 'Left click to create points, right click on any point to quit';
                newpnt = gd_setpoint(ax,promptxt,false);      %isxyz false 
                points = deletepoint(ax,points,delpnt,newpnt);
            end

        elseif strcmp(h_but.Tag,'Delete') 
            promptxt = 'Select point to Delete, right click on any point to quit';
            delpnt = gd_getpoint(ax,promptxt);   %get point to delete
            if ~isempty(delpnt)
                points = deletepoint(ax,points,delpnt,[]);  %delete the point                
            end

        elseif strcmp(h_but.Tag,'Quit') 
            %no longer used as a button - use close figure
            points = [];
            ok = 1;

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
function points = setPoints(ax,npts,promptxt)
    %use mouse to select points
    h_pnts = findobj(ax,'Tag','mypoints');
    delete(h_pnts)  %remove any existing points   
    hold on
    for i=1:npts
        if length(promptxt)~=npts  && ~isempty(promptxt)
            prompt = sprintf('%s %d',promptxt{1},i);
        else
            prompt = promptxt{i};
        end
        point = gd_setpoint(ax,prompt,false);
        if isempty(point)
            if isvalid(ax), ax.Title.String = 'Input cancelled'; end
            if ~exist('points','var'), points = []; end
            hold off
            return;   %user cancelled
        else
            points(i) = point; %#ok<AGROW> 
        end
    end
    ax.Title.String = 'Input complete';
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