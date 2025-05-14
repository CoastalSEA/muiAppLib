classdef (Abstract = true) PLinterface < handle
%
%-------class help---------------------------------------------------------
% NAME
%   PLinterface.m
% PURPOSE
%   Abstract class to provide a range of methods to manipulate points and
%   lines (usually in conjunction with the grid tools).
% USAGE
%   classes that inherit PLinterface define menus as a subset of the options 
%   defined in the setFigure function, or overload this
%   function if a more customised set of actions is required
%   Any class using the interface must define the Abstract properties for
%   outPoints, outLines and isXYZ.
% NOTES
%   inherits handle
%   the conventions used for points and lines are explained in the
%   documentation (see grid_point_line.m)
% SEE ALSO
%   the following external functions are made use of as part of the tool
%   set used to manipulate points and lines:
%   > point = gd_setpoint(ax,promptxt,isxyz) - create a single point
%   > points = gd_setpoints(ax,promptxt,isxyz) - create a set of points
%   > point = gd_getpoint(ax,promptxt) - select a single point
%   > pline = gd_setpline(ax,promptxt,isxyz) - create a line *********
%   > pline = gd_getpline(ax,promptxt,ispoints) - select a line *******
%
%   functions for converting between points, plines and cplines:
%   > lines = gd_points2lines(points,outype) - points to array of lines
%   > [points,outype] = gd_lines2points(lines) - lines to array of points 
%   > cplines = gd_plines2cplines(plines) - array of plines to a cell array of plines
%   > plines = gd_cplines2plines(cplines) - cell array of plines to a struct array
%
%   other utility functions:
%   > ax = gd_plotpoints(ax,points,tagname,type) - plot as selectable points or lines
%   > lineIndex = gd_findline(plines, queryPoint) - find line that includes queryPoint
%   
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%    
    properties (Hidden,Transient)
        Figure         %handle to class figure
        Axes           %handle to class axes
        Points = []    %holds the current set of points
        pLines = []    %holds the current set of lines
    end

    properties (Dependent)
       NaNpnts         %x,y(z) NaN points used to terminate lines
    end
    
    properties (Abstract)
        outPoints      %saved points for output
        outLines       %saved lines for output
        isXYZ          %logical flag, true if points are lines are to include z values
    end

    methods (Abstract)

    end
%%
    methods
        function obj = setFigure(obj,figtitle,tag,position)
            %initialise a figure and panel for PL classes
            %   figtitle - figure title
            %   tag - figure Tag name
            %   position - figure position; optional, default is [0.37,0.58,0.26,0.34] 
            if isempty(position), position = [0.37,0.58,0.26,0.34]; end

            h_fig = figure('Name',figtitle,'Tag',tag,'NextPlot','add',...
                           'Units','normalized','Position',position,...
                           'Visible','off');
            h_fig.MenuBar = 'none';   %reduce clutter but
            h_fig.ToolBar = 'figure'; %allow access to data tips and save tools
            obj.Figure = h_fig;
        end

%%
        function nanpnts = get.NaNpnts(obj)
            nanpnts = [];
            if isempty(obj.isXYZ)
                warndlg('isXYZ nhas not been set'); return; 
            else
                %initialise the nan points struct
                nanpnts = struct('x',NaN','y',NaN);
                if obj.isXYZ
                    nanpnts.z = NaN;
                end
            end
        end
    end
%%
    methods (Access=protected)
        function obj = setMenu(obj,menu,submenus)            
            %Define the context menu
            % menu - cell array of character vectors defining the main menu
            % submenus - additional submenus to those defined in getCallBacks
            %            defaults are provided for Figure, Point, Line
            if nargin<3, submenus = []; end
            cMenu = uicontextmenu(obj.Figure);

            % Add menu items
            for i=1:length(menu)
                cm(i) = uimenu(cMenu, 'Text', menu(i).label,'Tag','mainPLmenu');                    
            end

            cbs = getCallBacks(obj,submenus);
            % Add submenu items
            for i=1:length(cm)
                cmenu = cm(i);
                label = menu(i).label;
                idmain = find(ismember(cbs.Parent,label));
                for j=1:length(idmain)
                    callback = cbs.Callback{idmain(j)};
                    if isempty(callback)
                        submenu = uimenu(cmenu, 'Text',cbs.Label(idmain(j)),'Tag','subPLmenu');                                  
                        sublabel = cbs.Label(idmain(j));
                        idsub = find(ismember(cbs.Parent,sublabel));
                        for k=1:length(idsub)
                            callback = cbs.Callback{idsub(k)};
                            uimenu(submenu, 'Text',cbs.Label(idsub(k)),...
                                  'Callback',callback,'Tag','subsubPLmenu'); 
                        end
                    else
                        uimenu(cmenu, 'Text',cbs.Label(idmain(j)),...
                                  'Callback',callback,'Tag','subPLmenu'); 
                    end
                end
            end

            % Assign the context menu to the figure
            obj.Figure.UIContextMenu = cMenu;
        end
%%
%--------------------------------------------------------------------------
% Menu callback functions for Points
%--------------------------------------------------------------------------
         function addPoints(obj,~, ~)
            %callback function to add points
            promptxt = sprintf('Add points\nLeft click to create points, right click to finish');
            newpnts = gd_setpoints(obj.Axes,promptxt,'mypoints',obj.isXYZ);%get points to add
            obj.Points = [obj.Points,newpnts];
            obj.Axes.Title.String = 'Select menu option';
         end
         
 %%       
        function editPoint(obj,~, ~)
            %callback function to edit points
            resetMenu(obj)
            promptxt = sprintf('Edit point\nLeft click to select point, right click to quit');
            editpnt = gd_getpoint(obj.Axes,promptxt,'mypoints');           %get points to edit
            while ~isempty(editpnt)
                newpnt = gd_setpoint(obj.Axes,promptxt,'mypoints',obj.isXYZ);   
                if ~isempty(newpnt)
                    obj.Points = editApoint(obj,'Points',editpnt,newpnt,'mypoints');
                end
                resetPoints(obj);
                editpnt = gd_getpoint(obj.Axes,promptxt,'mypoints'); 
            end
            resetMenu(obj,true)
        end

 %%       
        function deletePoint(obj,~, ~)
            %callback function to delete points
            resetMenu(obj)
            promptxt = sprintf('Delete point\nSelect point to Delete, right click to quit');
            delpnt = gd_getpoint(obj.Axes,promptxt,'mypoints');            %get point to delete
            while ~isempty(delpnt)
                obj.Points = deleteApoint(obj,'Points',delpnt);            %delete the point 
                resetPoints(obj);
                delpnt = gd_getpoint(obj.Axes,promptxt,'mypoints');        %get point to delete
            end
            resetMenu(obj,true)
        end
        
%%
%--------------------------------------------------------------------------
% Menu callback functions for Lines
%--------------------------------------------------------------------------      
        function addLine(obj,~, ~)
            %callback function to add lines
            resetMenu(obj)
            ax = obj.Axes;
            promptxt = sprintf('Add lines\nLeft click to create points, right click to finish');
            [newpnts,H] = gd_setpoints(ax,promptxt,'mylines',obj.isXYZ);   %get points to add
            while ~isempty(newpnts)     
                NewLine =  [newpnts,obj.NaNpnts];
                obj.pLines = [obj.pLines,NewLine];
                delete(H);                                                 %delete digitising points
                gd_plotpoints(ax,NewLine,'mylines',2,obj.isXYZ);           %plot as line
                [newpnts,H] = gd_setpoints(ax,promptxt,'mylines',obj.isXYZ);%get points to add
            end
            resetMenu(obj,false)
        end

%%
        function editLine(obj,~, ~)
            %callback function to edit lines
            resetMenu(obj)
            promptxt1 = sprintf('Edit line\nSelect line to Edit, right click on any line to quit');
            promptxt2 = sprintf('Edit point\nLeft click to select point, right click to quit');

            [edline,Hl] = gd_getpline(obj.Axes,promptxt1,'mylines');       %get line to delete
            while ~isempty(edline)
                obj.Axes = gd_plotpoints(obj.Axes,edline,'edpnt',1);       %plot points for line to be edited                
                [editpnt,H] = gd_getpoint(obj.Axes,promptxt2,'edpnt');     %get points to edit
                while ~isempty(editpnt)
                    newpnt = gd_setpoint(obj.Axes,promptxt2,'edpnt',obj.isXYZ); %define new point
                    if ~isempty(newpnt)
                        obj.pLines = editAline(obj,'pLines',editpnt,newpnt,Hl); %replace old with new
                        % labelLines(obj,obj.pLines);                      %update line numbers
                    end
                    delete(H)
                    resetPoints(obj);
                    [editpnt,H] = gd_getpoint(obj.Axes,promptxt2,'edpnt'); %get another point and loop
                end
                clearGraphics(obj,{'edpnt'})                               %clear any edit graphics
                resetLines(obj);                                           %reset line state
                [edline,Hl] = gd_getpline(obj.Axes,promptxt1,'mylines');   %get another line to edit
            end
            resetMenu(obj,false)
        end

%%
        function Extend(obj,~, ~)
            %callback function to add points to the end of a line
            resetMenu(obj)
            promptxt1 = sprintf('Extend line\nSelect point at end of line to extend');   
            promptxt2 = sprintf('Extend line\nLeft click to create points, right click to finish\nFirst new point should be closest to selected end point');
            
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'endpoints',4);            
            endpnt = gd_getpoint(obj.Axes,promptxt1,'endpoints'); 
            while ~isempty(endpnt)
                if ~isempty(endpnt)
                    newpnts = gd_setpoints(obj.Axes,promptxt2,'endpoints',obj.isXYZ);  
                    if ~isempty(newpnts)
                        obj.pLines = extendAline(obj,'pLines',newpnts,endpnt);     
                    end
                end  
                clearGraphics(obj,{'mylines','endpoints'})
                obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2); 
                obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'endpoints',4); 
                endpnt = gd_getpoint(obj.Axes,promptxt1,'endpoints');
            end

            clearGraphics(obj,{'endpoints'})
            resetMenu(obj,false)
        end
        
%%        
        function Insert(obj,~, ~)
            %callback function to insert points into a line
            resetMenu(obj)
            promptxt1 = sprintf('Insert points\nSelect line, right click on any line to quit');
            promptxt2 = sprintf('Insert points\nLeft click to create points, right click to finish');
            prompts = {'Select first point','Select adjacent point'};
        
            insline = gd_getpline(obj.Axes,promptxt1,'mylines');           %get line to use
            while ~isempty(insline)
                [obj.Axes,Hp] = gd_plotpoints(obj.Axes,insline,'insertpnts',1);
                inspnts = getPoints(obj,prompts,'insertpnts',Hp);
                if ~isempty(inspnts)
                    newpnts = gd_setpoints(obj.Axes,promptxt2,'insertpnts',obj.isXYZ); %points to insert
                    if ~isempty(newpnts)
                        obj.pLines = insertPoints(obj,'pLines',inspnts,newpnts); %add points
                    end
                end        
                clearGraphics(obj,{'mylines','insertpnts'})
                obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);
                insline = gd_getpline(obj.Axes,promptxt1,'mylines');       %get line to use
            end
            clearGraphics(obj,{'insertpnts'})
            resetMenu(obj,false)
        end
        
%%        
        function Join(obj,~, ~)
            %callback function to join two lines
            resetMenu(obj)
            prompt1 = sprintf('Join line\nSelect first point of join');   
            prompt2 = sprintf('Join line\nSelect second point of join');
            prompts = {prompt1,prompt2};

            [obj.Axes,Hp] = gd_plotpoints(obj.Axes,obj.pLines,'joinpnts',4); 
            joinpnts = getPoints(obj,prompts,'joinpnts',Hp);
            while ~isempty(joinpnts)                
                obj.pLines = joinLines(obj,'pLines',joinpnts);   
                clearGraphics(obj,{'mylines','joinpnts'}) 
                obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);
                [obj.Axes,Hp] = gd_plotpoints(obj.Axes,obj.pLines,'joinpnts',4); 
                joinpnts = getPoints(obj,prompts,'joinpnts',Hp);               
            end
            clearGraphics(obj,{'joinpnts'})
            resetMenu(obj,false)
        end
 
%%
        function Split(obj,~, ~)
            %callback function to split a line into two lines
            resetMenu(obj)
            promptxt1 = sprintf('Split line\nSelect line to split, right click on any line to quit');
            promptxt2 = sprintf('Split point\nLeft click to select point, right click to quit');

            spline = gd_getpline(obj.Axes,promptxt1,'mylines');            %get line to split
            while ~isempty(spline)
                obj.Axes = gd_plotpoints(obj.Axes,spline,'splitpnt',1);    %plot points for line to be split
                splitpnt = gd_getpoint(obj.Axes,promptxt2,'splitpnt'); %get split point
                if ~isempty(splitpnt)                   
                    obj.pLines = splitLine(obj,'pLines',splitpnt);         %split the line
                    clearGraphics(obj,{'mylines','splitpnt'}) 
                    obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);
                end
                clearGraphics(obj,{'splitpnt'})                            %clear any split graphics
                resetLines(obj);                                           %reset line state
                spline = gd_getpline(obj.Axes,promptxt1,'mylines');        %get line to split
            end
            resetMenu(obj,false)
        end

%%
        function delLinePoint(obj,~,~)
            %callback function to delete points in a line
            resetMenu(obj)
            promptxt1 = sprintf('Delete points in a line\nSelect line, right click on any line to quit');
            promptxt2 = sprintf('Delete points in a line\nLeft click to select point, right click to quit');

            [edline,~] = gd_getpline(obj.Axes,promptxt1,'mylines');       %get line to delete point on
            while ~isempty(edline)
                obj.Axes = gd_plotpoints(obj.Axes,edline,'delpnt',1);      %plot points for line to be edited                
                [delpnt,H] = gd_getpoint(obj.Axes,promptxt2,'delpnt');     %get point to delete
                while ~isempty(delpnt)
                    obj.pLines = deleteApoint(obj,'pLines',delpnt);        %delete point in line
                    delete(H)
                    resetPoints(obj);
                    [delpnt,H] = gd_getpoint(obj.Axes,promptxt2,'delpnt'); %get another point and loop
                end
                clearGraphics(obj,{'delpnt'})                              %clear any edit graphics
                resetLines(obj);                                           %reset line state
                [edline,~] = gd_getpline(obj.Axes,promptxt1,'mylines');   %get another line to edit
            end
            clearGraphics(obj,{'mylines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);     %2= plot as lines
            resetMenu(obj,false)
        end

%%
        function deleteLine(obj,~, ~)
            %callback function to delete lines
            resetMenu(obj)
            promptxt = sprintf('Delete line\nSelect line to Delete, right click on any line to quit');
            [deline,H] = gd_getpline(obj.Axes,promptxt,'mylines');         %get line to delete
            while ~isempty(deline)                
                [obj.pLines,isdel] = deleteAline(obj,'pLines',deline);       %delete the line
                if isdel, delete(H); end
                resetLines(obj);
                pause(1)
                [deline,H] = gd_getpline(obj.Axes,promptxt,'mylines');     %get line to delete
            end
            resetMenu(obj,false)
        end

%%
        function Resample(obj,~,~)
            %resample the lines at user specified interval
            cint = PLinterface.setInterval();
            obj.pLines = resampleLines(obj,cint);
            clearGraphics(obj,{'mylines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);     %2= plot as lines
        end

%%
        function Smooth(obj,~,~)
            %smooth the lines using moving average of sgolay methods
            inp = PLinterface.setSmoothingInputs();
            obj.pLines = smoothLines(obj,'pLines',inp);
            clearGraphics(obj,{'mylines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);     %2= plot as lines            
        end

%%
        function Redraw(obj,~,~)
            %redraw the current set of points and lines
            ax = obj.Axes;            
            clearGraphics(obj,{'mypoints','mylines','mytext'});

            if ~isempty(obj.pLines)   
                ax = gd_plotpoints(ax,obj.pLines,'mylines',2); %set line          
            end
  
            if ~isempty(obj.Points) 
                gd_plotpoints(ax,obj.Points,'mypoints',1);     %set points
            end
        end  

%%
        function toggleView(obj,~,~)
            %switch a line of the lines and points defined with the start point
            %of each line emphasised with a circle marker.
            ax = obj.Axes;            

            hline = findobj(ax,'Tag','mylines');
            if isempty(hline)           %toggle lines on      
                ax = gd_plotpoints(ax,obj.pLines,'mylines',2);   %set line          
            else                        %toggle line and points off
                clearGraphics(obj,{'mylines'});
            end
  
            hpnts = findobj(ax,'Tag','mypoints');
            if isempty(hpnts)           %toggle points on
                gd_plotpoints(ax,obj.Points,'mypoints',1);  %set points
            else
                clearGraphics(obj,{'mypoints'});
            end

            hpnts = findobj(ax,'Tag','mytext');
            if isempty(hpnts)           %toggle points on
                gd_plotpoints(ax,obj.Points,'mytext',3);  %set points
            else
                clearGraphics(obj,{'mytext'});
            end         
        end       

%%
        function Undo(obj,~,~)
            %function to clear current points and linesand reset to saved
            %points and lines
            obj.Points = obj.outPoints;
            obj.pLines = obj.outLines;
            clearGraphics(obj,{'mypoints','mylines'});
            toggleView(obj,[],[]);
        end

%%
        function ClearAll(obj,~,~)
            %clear all existing user graphics
            obj.Points = [];
            obj.pLines = [];
            clearGraphics(obj,{'mypoints','mylines','mytext'});
        end

%%
        function Distance(obj,~,~)
            %show distance between two points
            prompt1 = sprintf('Distance\nSelect first point');   
            prompt2 = sprintf('Distance\nSelect second point');
            prompts = {prompt1,prompt2};         
            [distpnts,H] = setPoints(obj,prompts,'distpnts');
            while ~isempty(distpnts)
                distance = obj.lineLength(distpnts(1),distpnts(2));
                getdialog(sprintf('Distance = %.2f',distance));
                delete(H)
                [distpnts,H] = setPoints(obj,prompts,'distpnts');
            end
            clearGraphics(obj,{'distpnts'});
        end

%%
        function Save(obj,~, ~)
            %save the current set of points and lines
            obj.outPoints = obj.Points;
            obj.outLines = obj.pLines;
        end

%%
        function SaveExit(obj,~, ~)
            %set waitfor trigger in waitForFigure
            obj.Figure.Tag = 'SaveExit';
        end

%%
        function Quit(obj,~, ~)
            %set waitfor trigger in waitForFigure
            obj.Figure.Tag = 'Quit';
        end

%%
        function Exit(obj,~, ~)
            %set waitfor trigger in waitForFigure
            obj.Figure.Tag = 'Exit';  %no data to save in some use cases
        end

%%
        function obj = waitForFigure(obj)
            %control exiting and closing of figure
            ok = 0;
            while ok<1
                waitfor(obj.Figure,'Tag')

                if ~ishandle(obj.Figure) %this handles the user deleting figure window
                    obj.outPoints = []; obj.outLines = [];

                elseif strcmp(obj.Figure.Tag,'SaveExit')
                    obj.outPoints = obj.Points;
                    obj.outLines = obj.pLines;

                elseif strcmp(obj.Figure.Tag,'Quit')
                    if ~isempty(obj.Points) || ~isempty(obj.pLines)
                        answer = questdlg('Save input?','PLquit','Yes','No','No');
                        if strcmp(answer,'Yes')
                            obj.outPoints = obj.Points;
                            obj.outLines = obj.pLines;
                        else
                            obj.outPoints = [];  obj.outLines = [];
                        end                                                       
                    end
                elseif strcmp(obj.Figure.Tag,'Exit')
                    obj.outPoints = []; obj.outLines = [];
                end
                ok = 1;
            end
        end

%%
        function [points,H] = setPoints(obj,promptxt,tagname)
            %define a set number of points and return in pnt struct:
            % promptxt is a cell array with prompts for each call and
            % determines the number of points returned
            npnts = length(promptxt);
            points(1,npnts) = struct('x',[],'y',[]);
            H(1,npnts) = gobjects;
            for i=1:npnts
                [Pi,Hi] = gd_setpoint(obj.Axes,promptxt{i},tagname,obj.isXYZ);
                if isempty(Pi) && i==1, points = []; return; end  %user cancels on first point
                points(i) = Pi;  H(i) = Hi;
            end
        end

%%
        function [points,H] = getPoints(obj,promptxt,tagname,H)
            %find a set number of points and return in pnt struct:
            % promptxt is a cell array with prompts for each call and
            % determines the number of points returned

            %get the start and adjacent points for the insertion of new points
            npnts = length(promptxt);
            points(1,npnts) = struct('x',[],'y',[]);
            for i=1:npnts
                [pnt,H]  = gd_getpoint(obj.Axes,promptxt{i},tagname);
                if isempty(pnt)
                    points = []; return;
                else
                    points(i) = pnt;
                    H.UserData = int32(0);   %reset the selected point
                end
            end

            delete(H)
            hold on           %plot the insertion end points
            H = plot(obj.Axes,[points(:).x],[points(:).y],'or','MarkerSize',4,...
                                              'MarkerFaceColor','w','Tag',tagname);                                                  
            hold off   
        end

%%
        function outype = setInLines(obj,inlines)
            %plot imported points and lines as active graphics
            [points,out_type] = gd_lines2points(inlines); 
   
            if any(isnan([points(:).x]))
                obj.pLines = points;   
                obj.outLines = points;                   
                obj.Axes = gd_plotpoints(obj.Axes,points,'mylines',2);     %2= plot as lines
                outype.lines = out_type;
            else
                if length(points)>5000
                    getdialog(sprintf('Large number of points (N=%d)\nLoading linework my take some time',length(points)));
                end
                obj.Points = points;   
                obj.outPoints = points;                    
                obj.Axes = gd_plotpoints(obj.Axes,points,'mypoints',1);    %1= plot as points
                outype.points = out_type;
            end
        end

        %%
        function outype = getInLines(obj,inlines)
            %use input values to determine what to do
            if isempty(inlines) 
                outype = struct('points',2,'lines',2);  
            elseif isnumeric(inlines) && isscalar(inlines)   
                %handle call to function with no lines
                outype = struct('points',inlines,'lines',inlines);   
            else                            %plot imported points and lines
                if isfield(inlines,'x')     %single input type
                    outype = setInLines(obj,inlines);
                else                        %lines and points stucts
                    fnames = fieldnames(inlines);
                    for i=1:length(fnames)
                        inpoints = inlines.(fnames{i});
                        outype = setInLines(obj,inpoints);
                    end
                end
            end
        end
%%
%--------------------------------------------------------------------------
% Utility functions to implement various actions
%--------------------------------------------------------------------------
        function points = editApoint(obj,type,edpoint,newpnt,tagname)
            %delete point defined by delpnt if newpnt is empty, otherwise edit
            %point to new value as defined in newpnt
            points = obj.(type);
            %idp = [points(:).x]==edpoint.x & [points(:).y]==edpoint.y;
            idp = find([points(:).x]==edpoint.x & [points(:).y]==edpoint.y,1,'first'); 
            if isempty(idp), return; end
            points(idp).x = newpnt.x;
            points(idp).y = newpnt.y;

            h_pnts = findobj(obj.Axes,'Tag',tagname);
            idx = [h_pnts(:).XData]==edpoint.x & [h_pnts(:).YData]==edpoint.y; 
            h_pnts(idx).XData = newpnt.x;
            h_pnts(idx).YData = newpnt.y;                 
        end

%%
        function points = deleteApoint(obj,type,delpoint)
            %delete point defined by delpnt if newpnt is empty, otherwise edit
            %point to new value as defined in newpnt
            points = obj.(type);
            %idp = [points(:).x]==delpoint.x & [points(:).y]==delpoint.y;
            idp = find([points(:).x]==delpoint.x & [points(:).y]==delpoint.y,1,'first'); 
            if isempty(idp), return; end
            answer = questdlg('Confirm deletion','Delete point','Yes','No','Yes');
            if strcmp(answer,'Yes')
                points(idp) = [];
                if strcmp(type,'pLines')
                    h_lns = findobj(obj.Axes,'Tag','mylines');
                    for i=1:length(h_lns)
                        idx = [h_lns(i).XData]==delpoint.x & [h_lns(i).YData]==delpoint.y;
                        if any(idx)
                            h_lns(i).XData(idx) =[]; 
                            h_lns(i).YData(idx) =[];
                        end
                    end
                else
                    h_pnts = findobj(obj.Axes,'Tag','mypoints');
                    idx = [h_pnts(:).XData]==delpoint.x & [h_pnts(:).YData]==delpoint.y;
                    delete([h_pnts(idx)]);  %remove any existing points
                end
            end
        end

%%
        function plines = editAline(obj,type,edpoint,newpnt,hline)
            %delete point defined by edpoint and edit point to new value 
            %as defined in newpnt
            plines = obj.(type);
            %idp = [plines(:).x]==edpoint.x & [plines(:).y]==edpoint.y; 
            idp = find([plines(:).x]==edpoint.x & [plines(:).y]==edpoint.y,1,'first'); 
            if isempty(idp), return; end
            plines(idp).x = newpnt.x;
            plines(idp).y = newpnt.y;

            %hline = findobj(obj.Axes,'Tag',tagname);
            idx = [hline.XData]==edpoint.x & [hline.YData]==edpoint.y; 
            hline.XData(idx) = newpnt.x;
            hline.YData(idx) = newpnt.y;                 
        end

        %%
        function points = extendAline(obj,type,newpnts,endpnt)
            %find end to add points to and extend existing line
            points = obj.(type);
            idp = find([points(:).x]==endpnt.x & [points(:).y]==endpnt.y,1,'first'); 
            if isempty(idp), return; end

            %check direction of new points relative to existing line
            isstart = idp==1 || isnan(points(idp-1).x);
            newpnts = obj.checkDirection(endpnt,newpnts,isstart);
            %insert the new points (cases ensure NaNs are maintained)
            if idp==1
                points = [newpnts,points];         %in front of all lines        
            elseif isnan(points(idp-1).x)
                points = [points(1:idp-1),newpnts,points(idp:end)];
            else
                points = [points(1:idp),newpnts,points(idp+1:end)];
            end
        end

%%
        function points = insertPoints(obj,type,inspnts,newpnts)
            %insert additional points, 'newpnts', between the selected points,
            %'pnt', in the digitised points vector,'points'.
            points = obj.(type);
            idpos1 = find([points(:).x]==inspnts(1).x & [points(:).y]==inspnts(1).y,1,'first'); 
            idpos2 = find([points(:).x]==inspnts(2).x & [points(:).y]==inspnts(2).y,1,'first'); 
            if isempty(idpos1) || isempty(idpos2), return; end

            %polygons can have the same point twice
            idx = abs(idpos1-idpos2)==1;  %finds the adjacent point, assuming more than 3 points
            if length(idpos1)>1
                idpos1 = idpos1(idx);
            elseif length(idpos2)>1
                idpos2 = idpos2(idx);
            end
            idpos = [idpos1,idpos2];
            if idpos(2)<idpos(1), idpos = fliplr(idpos); end
        
            %check direction of new points relative to existing line
            isstart = idpos(1)==1 || isnan(points(idpos(1)-1).x);
            newpnts = obj.checkDirection(points(idpos(1)),newpnts,isstart);
        
            %insert the new points in the correct interval for the line
            %to be in the right order they need to be created in the same direction
            %as the initial line is being created - this is NOT checked.
            points = [points(1:idpos(1)),newpnts,points(idpos(2):end)];
        end

%%
        function points = joinLines(obj,type,joinpnts)
            %join two lines based on selected end points
            points = obj.(type);
            msg = {'Cannot join line to itself','Creating a polygon'};
            idi = find([points(:).x]==joinpnts(1).x & [points(:).y]==joinpnts(1).y,1,'first'); 
            idj = find([points(:).x]==joinpnts(2).x & [points(:).y]==joinpnts(2).y,1,'first'); 
             if isempty(idi) || isempty(idj),return; end

            %check not trying to join a line to itself
            if idi==idj,titleWarning(obj,msg{1}); return; end

            %handle line being made into a polygon            
            Lineidi = gd_findline(points,joinpnts(1));
            Lineidj = gd_findline(points,joinpnts(2));
            if Lineidi==Lineidj
                idp = sort([idi,idj]);                               
                points = [points(1:idp(2)),points(idp(1)),points(idp(2)+1:end)];
                titleWarning(obj,msg{2}) 
                return; 
            end  

            %get second line and check direction of points relative to first line
            idN = [0,find(isnan([points(:).x]))]; 
            idLj =idN(Lineidj)+1:idN(Lineidj+1);
            newpnts = points(idLj(1:end-1));
            if isnan(points(idj+1).x)     %if join point is at end of line
                newpnts = fliplr(newpnts);%reverse direction of line
            end
        
            %remove the existing second line from the points set
            obj.pLines(idLj) = [];
        
            %to join the two lines, treat newpnts as a new points to be added
            points = extendAline(obj,type,newpnts,joinpnts(1));
        end

%%
        function points = splitLine(obj,type,initpnt)
            %split a line at the defined point
            points = obj.(type);
            idi = find([points(:).x]==initpnt.x & [points(:).y]==initpnt.y,1,'first'); 
            if isempty(idi), return; end
            points = [points(1:idi),obj.NaNpnts,points(idi+1:end)];
        end

%%
        function [plines,isdel] = deleteAline(obj,type,deline)
            %delete a line from a set of lines
            plines = obj.(type);
            idl = gd_findline(plines, deline(1));
            if idl<1, return; end           %line not found
            answer = questdlg('Confirm deletion','Delete line','Yes','No','Yes');
            if strcmp(answer,'Yes')
                cplines = gd_plines2cplines(plines);
                cplines(idl) = [];          %delete from lines cell array
                plines = gd_cplines2plines(cplines);
                isdel = true;
            else
                isdel = false;
            end
        end

%%
        function plines = smoothLines(obj,type,inp)
            %smooth the centre-line with option to use or revert to original
            clearGraphics(obj,{'smoothlines'});  %remove any existing smoothed lines
            plines = obj.(type);
            if inp.idm==0, method = 'movmean'; else, method = 'sgolay'; end
            blines = gd_points2lines(plines,2);
            blines = gd_smoothlines(blines,method,inp.win,inp.deg,inp.npnts);   
            hold on
            plot(obj.Axes,blines.x,blines.y,'-.g','LineWidth',1,'Tag','smoothlines')
            hold off

            %option to use or revert to original
            answer = questdlg('Accept new line or retain previous version?',...
                              'Sections','Accept','Reject','Reject');
            if strcmp(answer,'Accept')
                plines = gd_lines2points(blines); 
            else
                clearGraphics(obj,{'smoothlines'});  %remove smoothed lines
            end
        end

%%
        function plines = resampleLines(obj,cint)
            %resample the contour lines at intervals of cint
            %idN = [0;find(isnan(obj.pLines(:).x))];
            cplines = gd_plines2cplines(obj.pLines);
            nlines = length(cplines);
            newlines = [];  
            for i=1:nlines                
                cline = cplines{1,i};                         %extract line   
                cline = gd_points2lines(cline(1:end-1),1);    %convert to matrix omit trailing NaN
                clength = sum(vecnorm(diff(cline),2,2));      %cline is a column vector [Nx2]                
                if clength>cint                                   %must be > cint                    
                    cpoints = round(clength/cint);                %number of points in new line
                    newcline = curvspace(cline,cpoints);
                    if size(newcline,1)>2                         %trap single point lines
                        newlines = [newlines;newcline;[NaN,NaN]]; %#ok<AGROW>
                    end
                end                
            end
            plines = gd_lines2points(newlines);
        end

%%
        function labelLines(obj,cplines)
            %add labels to lines in plines
            clearGraphics(obj,{'mytext'});
            if ~iscell(cplines), cplines = gd_plines2cplines(cplines); end
            for i=1:length(cplines)
                gd_plotpoints(obj.Axes,cplines{i},num2str(i),3); %set labels
            end
        end

%%
        function calltable = getCallBacks(obj,submenus)
            %get the default menus for Point,Line and Figure and append
            %any user defined menus
            varnames = {'Parent','Callback','Label'};
            
            %default Figure menu variables
            stext = ["Redraw";"Undo";"Clear all";"Distance";"Save";"Save & Exit";"Quit"];
            scall = {@obj.Redraw; @obj.Undo; @obj.ClearAll; @obj.Distance; ...
                                  @obj.Save; @obj.SaveExit; @obj.Quit}; 
            nrec = length(stext);
            spart = repmat("Figure",nrec,1);
            calltable = table(spart,scall,stext,'VariableNames',varnames);

            %default point menu variables
            stext = ["Add";"Edit";"Delete"];
            scall = {@obj.addPoints; @obj.editPoint; @obj.deletePoint}; 
            nrec = length(stext);
            spart = repmat("Point",nrec,1);
            call{1} = table(spart,scall,stext,'VariableNames',varnames);

            %default line menu variables
            stext = ["Add";"Edit";"Extend";"Insert";"Join";"Split";"Delete"];
            scall = {@obj.addLine; @obj.editLine; @obj.Extend; @obj.Insert;...
                                     @obj.Join; @obj.Split; []}; 
            nrec = length(stext);
            spart = repmat("Line",nrec,1);
            call{2} = table(spart,scall,stext,'VariableNames',varnames);

            %default line submenu variables
            stext = ["Point";"Line"];
            scall = {@obj.delLinePoint; @obj.deleteLine};
            nrec = length(stext);
            spart = repmat("Delete",nrec,1);
            call{3} = table(spart,scall,stext,'VariableNames',varnames);

            for i=1:length(call)
                calltable = [calltable;call{i}]; %#ok<AGROW> 
            end
            %add any user defined menu options
            calltable = [calltable;submenus];
        end
%--------------------------------------------------------------------------
% Reset and clear graphics functions
%--------------------------------------------------------------------------
        function resetMenu(obj,ispoints)
            %switch the visbility of the graphics menu on and off
            if nargin<2, ispoints = false; end

            hm = findobj(obj.Figure,'Tag','mainPLmenu');
            if isempty(hm), return; end %user closes figure during operation
            if strcmp(hm(1).Visible,'on')
                [hm(:).Visible] = deal('off');
            else
                [hm(:).Visible] = deal('on');
                 obj.Axes.Title.String = 'Select menu option';
                 if ispoints
                     resetPoints(obj);
                 else
                     resetLines(obj)
                 end
            end
        end

%%
        function resetPoints(obj)
            %gd_getpoint sets point UserData and color when clicked on. Reset in 
            %case user clicked on points without making an action selection
            h_pnts = findobj(obj.Axes,'Tag','mypoints');
            if isempty(h_pnts), return; end
            idx = [h_pnts.UserData]>0;
            if any(idx)
                [h_pnts(idx).UserData] = deal(int32(0));
                [h_pnts(idx).Color] = deal([0,0,0]); 
            end
        end

%%
        function resetLines(obj)
            %gd_getpline sets point UserData and color when clicked on. Reset in 
            %case user clicked on points without making an action selection
            h_lines = findobj(obj.Axes,'Tag','mylines');
            if isempty(h_lines), return; end
            idx = [h_lines.UserData]>0;
            if any(idx)
                [h_lines(idx).UserData] = deal(int32(0));
                [h_lines(idx).Color] = deal([1,0,0]);
            end    
        end

%%
        function clearGraphics(obj,types)
            % delete one or line types based on Tag names
            for i=1:length(types)
                h_lines = findobj(obj.Axes,'Tag',types{i});
                delete(h_lines)
            end
        end

        %%
        function titleWarning(obj,msgtxt)
            %use title to indicate an error in the process
            obj.Axes.Title.String = msgtxt;  
            obj.Axes.Title.Color = 'r';
            pause(2);
            obj.Axes.Title.Color = 'k';
        end
    end
%--------------------------------------------------------------------------
% Static utility functions
%--------------------------------------------------------------------------
    methods (Static)
        function newpnts = checkDirection(endpnt,newpnts,isstart)  
            %check direction of new points relative to existing line
            xlen = [endpnt.x-newpnts(1).x,endpnt.x-newpnts(end).x];
            ylen = [endpnt.y-newpnts(1).y,endpnt.y-newpnts(end).y];
            [~,ide] = min(hypot(xlen,ylen)); %index of nearest point

            if (ide>1 && ~isstart) || (ide==1 && isstart)
                %flip lines if first new point is further away and is not
                %being added to the start of a line, or is the closest and
                %is being added to the start of a line
                newpnts = fliplr(newpnts);
            end
        end

%%
        function distance = lineLength(pnt1,pnt2)
            %find length of line between 2 points
            xlen = pnt1.x-pnt2.x;
            ylen = pnt1.y-pnt2.y;
            distance = hypot(xlen,ylen);
        end

%%
        function [isNear,idP,distances] = isPointNearLine(Points,Point,tol)
            % Calculate the distance from the point to points on line
            % also returns sorted indices and distances to all Points
            xlen = [Points(:).x]-Point.x;
            ylen = [Points(:).y]-Point.y;
            [distances,idP] = sort(hypot(xlen,ylen)); %indices of nearest points   
            isNear = any(distances < tol); % Threshold for proximity
        end

%%
        function zlevel = setLevel()
            %prompt user to set the level for the contour to be extracted
            promptxt = sprintf('Elevation of shoreline contour\n(or NaN for NaN mask):');
            inp = inputdlg({promptxt},'Boundary',1,{'NaN'});
            if isempty(inp), zlevel = []; return; end  %user cancelled
            zlevel = str2double(inp{1});
        end

%%
        function cint = setInterval()
            %prompt user to set point spacing interval for ampling along a
            %line.
            promptxt = sprintf('Sampling interval (m)');
            inp = inputdlg({promptxt},'Boundary',1,{'100'});
            if isempty(inp), cint = []; return; end  %user cancelled
            cint = str2double(inp{1});
        end

%%
        function inp = setSmoothingInputs()
            %prompt user for the smoothing input parameters
            %to define the method, window size, degree and mimimum number
            %of points to smooth
            promptxt = {'Method (0-moving av, 1-smooth)','Window size',...
                'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
            defaults = {'1','10','4','10'};
            input = inputdlg(promptxt,'Boundary',1,defaults);
            if isempty(input), inp = []; return; end      %user cancelled
            inp.idm= logical(str2double(input{1}));
            inp.win = str2num(input{2}); %#ok<ST2NM> allow input of vector [1x2]
            inp.deg = str2double(input{3});  
            inp.npnts = str2double(input{4});        
        end
    end
%--------------------------------------------------------------------------
% Sample function to illustrate defining additional submenus
%--------------------------------------------------------------------------
    methods (Access=protected)
       function calltable = setSubMenus(obj)
            %user defined menus to be appended to the Figure menu 
            varnames = {'Parent','Callback','Label'};
            %XXXX menu variables
            stext = ["A";"B";"C"];                                         %< Edit
            scall = {@obj.Save;@obj.SaveExit;@obj.Quit}; 
            nrec = length(stext);
            spart = repmat("Test1",nrec,1);                                %< Edit
            calltable = table(spart,scall,stext,'VariableNames',varnames);

            %YYYY point menu variables                                     %< Delete section if not required
            stext =  ["D";"E";"F"];                                        %< Edit
            scall = {@obj.addPoints;@obj.editPoint;@obj.deletePoint}; 
            nrec = length(stext);
            spart = repmat("Test2",nrec,1);                                %< Edit
            call{1} = table(spart,scall,stext,'VariableNames',varnames);

            for i=1:length(call)
                calltable = [calltable;call{i}]; %#ok<AGROW> 
            end
        end
    end
end