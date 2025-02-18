function [clines,props] = gd_setcentrelines(grid,mobj,props,inlines,isdel)
%
%-------function help------------------------------------------------------
% NAME
%   gd_setcentrelines.m
% PURPOSE
%   Accept figure to interactively create one or more centreline of a 
%   channel using function a_star to trace the shortest path between 
%   start and end points whilst finding the deepest points (ie a thalweg).
% USAGE
%   cline = gd_setcentrelines(grid,mobj,props,inlines,isdel)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   mobj - mui model instance
%   props - struct for maximum water level and depth exponent (maxwl,dexp,cint)
%   inlines - struct or table of x,y vectors to be edited or the format 
%             of output if no lines are being input (see Outputs for details)
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   clines - struct of x,y vectors defining the centre-lines      
%   props - as input with any updates
% NOTES
%   Finds indices of the grid used to find deepest points in channel and
%   resolution depends on the xy spacing of the grid used. 
% SEE ALSO
%   similar to gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    if nargin<5, isdel = false; end
    isxyz = false;
    figtitle = sprintf('Select points');
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Set line','Delete line','Smooth','Resample'...
        'Add point','Edit point','Insert point','Delete point','View','Use'};
    tooltips = {'Set a new centre-line',...
                'Clip sections to the boundary line',...
                'Delete a centre-line',...
                'Smooth the current contour lines',...
                'Resample the current contour lines at a specified interval',...
                'Add a point to an exsiting line (end point only)',...
                'Edit a point in an existing line',...
                'Insert one or more point between two existing points of a line',...
                'Delete a point in an existing line',...         
                'Toggle display of section lines on and off',...
                'Reset centre-line and sections to intial state',...
                'Use digitised points and exit. Close figure window to Quit without saving'};
    % position = [0.3,0.4,0.35,0.5];
    position = [0,0,1,1];
    paneltxt = 'Set and Edit centre-lines. Close window to Quit';
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position,0.8,tooltips);
    
    %initialise grid
    [ax,grid] = initialisegrid(h_plt,grid,props);
    
    if ~isempty(inlines)
        [c_plines,outype] = gd_lines2points(inlines); 
        clear inlines
        ax = plotCLines(ax,c_plines);
    else
        c_plines = []; outype = 2;
    end
   
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window 
            clines = [];
            return;

        elseif strcmp(h_but.Tag,'Set line') 
            promptxt = {'Select start of path','Select end of path'};
            [points,Hp] = setPoints(ax,2,promptxt);
            if ~isempty(points)                
                c_pline = getCentreLine(mobj,grid,props,points);
                if isempty(c_pline), continue; end
                delete(Hp)                                 %remove set points
                [ax,Hl] = checkPlot(ax,c_pline,c_plines);
                cline = cell2mat(struct2cell(c_pline)');   %struct to [Nx2] array                
                c_plines = setCentreLine(props,cline,c_plines);
                delete(Hl)                                 %remove check line
                ax = plotCLines(ax,c_plines);
            end

        elseif strcmp(h_but.Tag,'Delete line')
            clearLines(ax,{'clines'})
            ax = gd_plotpoints(ax,c_plines,'mylines',2);  %set lines
            promptxt = 'Select line to delete';
            delpline = gd_getpline(ax,promptxt);
            if ~isempty(delpline)
                c_plines =  deletelines(c_plines,delpline);  %delete the line   
            end

        elseif strcmp(h_but.Tag,'Smooth')
            %smooth the shoreline
            c_plines = smoothLines(ax,c_plines);
%             %get the user to define the upper limit to use for the hypsomety
%             promptxt = {'Method (0-moving av, 1-smooth)','Window size',...
%                         'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
%             defaults = {'0','10','4','10'};
%             inp = inputdlg(promptxt,'Boundary',1,defaults);
%             if ~isempty(inp)
%                 idm= logical(str2double(inp{1}));
%                 win = str2num(inp{2}); %#ok<ST2NM> allow input of vector [1x2]
%                 deg = str2double(inp{3});  
%                 npnts = str2double(inp{4});
%                 if idm==0, method = 'movmean'; else, method = 'sgolay'; end
%                 clines = gd_points2lines(c_plines,2); 
%                 clines = gd_smoothlines(clines,method,win,deg,npnts); 
%                 hold on
%                 plot(ax,clines.x,clines.y,'-.g','LineWidth',1,'Tag','slines')
%                 hold off
%                 [c_plines,outype] = gd_lines2points(inlines); 
%                 clear clines
%             end

        elseif strcmp(h_but.Tag,'Resample')
            %resample as selected interval
            props = setInterval(props);
            if ~isempty(props.cint)
                clines = gd_points2lines(c_plines,2);
                clines = resampleLines(clines,props.cint);
                %ax = plotLines(ax,blines);
                c_plines = gd_lines2points(clines); 
                gd_plotpoints(ax,c_plines,'mylines',2); 
            end

        elseif strcmp(h_but.Tag,'Add point') 
            clearLines(ax,{'clines'})
            ax = gd_plotpoints(ax,c_plines,'mypoints',1);  %set points
            promptxt = 'Select point at end of line to extend';
            endpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(endpnt)
                promptxt = sprintf('Left click to create points, right click to finish\nFirst new point should be closest to selected end point');
                newpnts = gd_setpoints(ax,promptxt,isxyz);  
                if ~isempty(newpnts)
                    c_plines = addpoints(c_plines,newpnts,endpnt,true);     
                end
            end

        elseif strcmp(h_but.Tag,'Edit point') 
            clearLines(ax,{'clines'})
            ax = gd_plotpoints(ax,c_plines,'mypoints',1);  %set points
            promptxt = 'Select point to edit';
            delpnt = gd_getpoint(ax,promptxt); 
            if ~isempty(delpnt)
                promptxt = 'Left click to create points, right click on any point to quit';
                newpnt = gd_setpoint(ax,promptxt,isxyz);   
                if ~isempty(newpnt)
                    points = deletepoint(ax,c_plines,delpnt,newpnt);
                end
            end

        elseif strcmp(h_but.Tag,'Insert point')
            clearLines(ax,{'clines'})
            ax = gd_plotpoints(ax,c_plines,'mypoints',1);  %set points
            pnt = startendpoints(ax,c_plines); %define insertion position
            if ~isempty(pnt)
                promptxt = 'Left click to create points, right click to finish';
                newpnts = gd_setpoints(ax,promptxt,isxyz);    %points to insert
                if ~isempty(newpnts)
                    c_plines = insertpoints(ax,c_plines,pnt,newpnts); %add points
                end
            end
        
        elseif strcmp(h_but.Tag,'Delete point') 
            clearLines(ax,{'clines'})
            ax = gd_plotpoints(ax,c_plines,'mypoints',1);  %set points            
            promptxt = 'Select point to Delete, right click on any point to quit';
            delpnt = gd_getpoint(ax,promptxt);   %get point to delete
            if ~isempty(delpnt)
                c_plines = deletepoint(ax,c_plines,delpnt,[]);  %delete the point                
            end

        elseif strcmp(h_but.Tag,'Reorder lines')            
            c_plines = gd_orderlines(ax,c_plines);
            ax = gd_plotpoints(ax,c_plines,'mypoints',1);

        elseif strcmp(h_but.Tag,'View') 
            ax = toggle_view(ax,c_plines);

        else   %user accepteds                          
            ok = 1; 
            delete(h_but);   %keep figure but delete buttons
            title(ax,'')
        end   
        h_but = resetbutton(ax,h_but);
        resetpoints(ax);
    end

    %convert format of output if required
    clines = gd_points2lines(c_plines,outype);  
    
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
        [h_pnts(idx).UserData] = deal(int32(0));
        [h_pnts(idx).Color] = deal([0,0,0]); 
    end 
end

%%
function [points,Hp] = setPoints(ax,npts,promptxt)
    %use mouse to select npts number of points
    h_pnts = findobj(ax,'Tag','mypoints');
    delete(h_pnts)  %remove any existing points   
    hold on
    for i=1:npts
        if length(promptxt)~=npts  && ~isempty(promptxt)
            prompt = sprintf('%s %d',promptxt{1},i);
        else
            prompt = promptxt{i};
        end
        [point,H] = gd_setpoint(ax,prompt,false);
        if isempty(point)
            if isvalid(ax), ax.Title.String = 'Input cancelled'; end
            if ~exist('points','var'), points = []; end
            hold off
            return;   %user cancelled
        else
            points(i) = point; %#ok<AGROW> 
            Hp(i) = H;         %#ok<AGROW> 
        end
    end
    ax.Title.String = 'Input complete';
    hold off
end

%%
function c_pline = getCentreLine(mobj,grid,props,points)
    %create a centreline of a channel using function a_star to trace the
    %shortest path between start and end points whilst finding the deepest
    %points (ie a thalweg)

    %index of nearest grid point to selected start end end points    
    start = dsearchn(grid.xy,[points(1).x,points(1).y]); 
    goal = dsearchn(grid.xy,[points(2).x,points(2).y]);
    
    hwb = progressbar(mobj,[],'Computing centre-line');
    %find the shortest path taking account of the cost (depths)
    %Z(Z>maxwl) = 0;
    costfnc = @(z) -(min(z,[],'all')-z).^props.dexp; %weighted inverse depth to favour staying deep
    thalweg = a_star(grid.water, costfnc(grid.Z), start, goal);
    [idy,idx] = ind2sub(size(grid.Z),thalweg); %convert indices to x,y subscripts
    progressbar(mobj,hwb);

    c_pline.x = [flipud(grid.x(idx));NaN]; %return points in order from start point
    c_pline.y = [flipud(grid.y(idy));NaN]; %as column vectors terminated with NaNs
end

%%
function plines = setCentreLine(props,nline,plines)
    %
    nline = nline(1:end-1,:);                     %remove trailing NaNs
    ansr = questdlg('Accept the centreline?','Centre-line','Yes','No','Yes');
    if strcmp(ansr,'Yes')
        %convert output to specified point pacing
        clength = sum(vecnorm(diff(nline),2,2));  %nline is a column vector [Nx2]
        cpoints = round(clength/props.cint);      %number of points in new line
        newpoints = curvspace(nline,cpoints);      %curvespace uses [Nx2]
        newline = [newpoints;[NaN,NaN]];     %restore trailing NaNs 
        pline = gd_lines2points(newline);
        plines = [plines,pline];
    end
end

%%
function plines = deletelines(plines,delpline)
    %delete line
    cplines = gd_plines2cplines(plines);
    nlines = length(cplines);
    npnts = length(delpline);
    for i=1:nlines
        aline = cplines{1,i};
        if length(aline)==npnts            
            isok = isequal([aline(1:end-1).x],[delpline(1:end-1).x]) && ...
                         isequal([aline(1:end-1).y],[delpline(1:end-1).y]);
            if isok 
                cplines{1,i} = [];
                break
            end
        end
    end
    plines = gd_cplines2plines(cplines);
end

%%
function points = addpoints(points,newpnts,endpnt,isadd)
    %find end to add points to and extend existing line
    %isadd true finds the point that is closest to the endpnt
    idp = find([points(:).x]==endpnt.x & [points(:).y]==endpnt.y); 

    %find the relative position on the line
    isfirst = false; isstart = false;
    if idp==1
        isfirst = true;                   %special case first point in lines
    else
        isstart = isnan(points(idp-1).x);  %end point is at start of a line
    end    
    isend = isnan(points(idp+1).x);        %end point is at end of a line
    
    if isadd
        %check direction of new points relative to existing line
        xlen = [endpnt.x-newpnts(1).x,endpnt.x-newpnts(end).x];
        ylen = [endpnt.y-newpnts(1).y,endpnt.y-newpnts(end).y];
        [~,ide] = min(hypot(xlen,ylen));   %index of nearest point
        if (ide>1 && isend) || (ide==1 && ~isend)
            newpnts = fliplr(newpnts);
        end
    end
    
    %insert the new points (cases ensure NaNs are maintained)
    if isfirst
        points = [newpnts,points];         %in front of all lines        
    elseif isstart
        points = [points(1:idp-1),newpnts,points(idp:end)];
    else
        points = [points(1:idp),newpnts,points(idp+1:end)];
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
    idposs = find([points(:).x]==pnt.start.x & [points(:).y]==pnt.start.y); 
    idposa = find([points(:).x]==pnt.adjacent.x & [points(:).y]==pnt.adjacent.y); 
    %polygons can have the same point twice
    idx = abs(idposs-idposa)==1;  %finds the adjacent point, assuming more than 3 points
    if length(idposs)>1
        idposs = idposs(idx);
    elseif length(idposa)>1
        idposa = idposa(idx);
    end
    idpos = [idposs,idposa];
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
function c_plines = smoothLines(ax,c_plines)
     %smooth the centre-line with option to use or revert to original
    clearLines(ax,{'smoothlines'})   %remove any existing smoothed lines
    %get the user to define the upper limit to use for the hypsomety
    promptxt = {'Method (0-moving av, 1-smooth)','Window size',...
        'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
    defaults = {'0','10','4','10'};
    inp = inputdlg(promptxt,'Boundary',1,defaults);
    if isempty(inp), return; end          %return line unchanged
    idm= logical(str2double(inp{1}));
    win = str2num(inp{2}); %#ok<ST2NM> allow input of vector [1x2]
    deg = str2double(inp{3});  
    npnts = str2double(inp{4});
    if idm==0, method = 'movmean'; else, method = 'sgolay'; end
    c_lines = gd_points2lines(c_plines,2);
    c_lines = gd_smoothlines(c_lines,method,win,deg,npnts);   
    hold on
    plot(ax,c_lines.x,c_lines.y,'-.g','LineWidth',1,'Tag','smoothlines')
    hold off
    answer = questdlg('Accept new line or retain previous version?',...
                      'Sections','Accept','Reject','Reject');

    if strcmp(answer,'Accept')
        c_plines = gd_lines2points(c_lines);
    end
    clearLines(ax,{'smoothlines','clines'}) %remove any existing centre-lines
    plotCLines(ax,c_plines);
end

%%
function props = setInterval(props)
    %prompt user to set point spacing interval for the contour
    promptxt = sprintf('Sampling interval (m)');
    inp = inputdlg({promptxt},'Boundary',1,{num2str(props.cint)});
    if isempty(inp), return; end  %user cancelled
    props.cint = str2double(inp{1});
end

%%
function clines = resampleLines(inlines,cint)
    %resample the contour lines at intervals of cint
    idN = [0;find(isnan(inlines.x))];
    [points,~] = gd_lines2points(inlines);            %convert to points
    nlines = [];
    for i=1:length(idN)-1
        cline = points(idN(i)+1:idN(i+1));            %extract line        
        cline = gd_points2lines(cline(1:end-1),1);    %convert to matrix omit trailing NaN
        if size(cline,1)>1                            %trap single point lines
            clength = sum(vecnorm(diff(cline),2,2));  %cline is a column vector [Nx2]
            cpoints = round(clength/cint);            %number of points in new line
            newcline = curvspace(cline,cpoints);
        else
            newcline = cline;
        end
        nlines = [nlines;newcline;[NaN,NaN]]; %#ok<AGROW>  
    end  
    clines.x = nlines(:,1);    %return struct of column vectors
    clines.y = nlines(:,2);
end

%% ------------------------------------------------------------------------
% Plotting functions
%%-------------------------------------------------------------------------
function [ax,grid] = initialisegrid(h_plt,grid,props)
    %plot grid and add varaiables needed for extracting centre-lines
    [X,Y] = meshgrid(grid.x,grid.y);
    N = numel(X);
    grid.xy = [reshape(X,[N,1]),reshape(Y,[N,1])];
    grid.Z = grid.z';
    %accessible map (water) and use -Z as the cost map
    water = true(size(grid.Z));
    water(isnan(grid.Z) | grid.Z>props.maxwl) = false;
    grid.water = water;
    gridmasked = grid;       
    gridmasked.z(~water') = NaN;
    ax = gd_plotgrid(h_plt,gridmasked);
end

%%
function [ax,H] = checkPlot(ax,cline,nlines)
    %plot base map of initial grid selection and defined mask
    % cline is xy single stuct line and nlines is a set of plines
    points = [cline.x(1),cline.y(1);cline.x(end-1),cline.y(end-1)];
    hold on
    hp = plot(ax,points(:,1),points(:,2),'ok','MarkerSize',8,...
                           'MarkerFaceColor','w','Tag','mypoints');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    hl = plot(ax,cline.x,cline.y,'-r','LineWidth',2,...
                                   'DisplayName','New line');
    if ~isempty(nlines)
        ax = plotCLines(ax,nlines);
    end
    hold off
    H = [hl;hp];
    title('Centre-line between defined start and end points')
end

%%
function ax = plotCLines(ax,points)
    %plot the centreline as non-pickable linework
    clearLines(ax,{'clines'})
    hold on
    plot(ax,[points(:).x],[points(:).y],'+k','MarkerSize',4,...
                              'PickableParts','none','Tag','clines');
    plot(ax,[points(:).x],[points(:).y],'ok','MarkerSize',3,...
                              'PickableParts','none','Tag','clines');
    hold off
end

%%
function ax = toggle_view(ax,cplines)
    %switch a line of the lines and points defined with the start point
    %of each line emphasised with a circle marker.
    hline = findobj(ax,'Tag','mylines');
    if ~iscell(cplines)
        cplines = gd_plines2cplines(cplines);
    end
    if isempty(hline)           %toggle line and points on
        for j=1:length(cplines)                         %call one at a time
            aline = (cplines{1,j});                     %to order numbering
            ax = gd_plotpoints(ax,aline,'mypoints',1);  %set points
            ax = gd_plotpoints(ax,aline,'mylines',2);   %set line 
            ax = gd_plotpoints(ax,aline,num2str(j),3);  %set labels
        end 
    else                        %toggle line and points off
        clearLines(ax,{'mylines','mypoints','mytext'});
    end
end

%%
function clearLines(ax,types)
    % delete one or line types based on Tag names
    for i=1:length(types)
        h_lines = findobj(ax,'Tag',types{i});
        delete(h_lines)
    end
end


