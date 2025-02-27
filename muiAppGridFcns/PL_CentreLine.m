classdef PL_CentreLine < PLinterface
%
%-------class help---------------------------------------------------------
% NAME
%   PL_CentreLine.m
% PURPOSE
%   Class to extract valley/channel centre-line
% USAGE
%  
% SEE ALSO
%   inherits PLinterface
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%    
    properties (Transient)
        %set Abstract properties for PLinterface
        outPoints = [] %not used but has to be set        
        isXYZ = false  %not used but has to be set
        outLines = []
        Grid           %imported grid for use in callback functions
        Props          %imported props for use in callback functions
    end
       
    methods
        function obj = PL_CentreLine(figtitle,tag,position)          
            %constructor code: 
            if nargin<3, position = []; end
            obj = setFigure(obj,figtitle,tag,position);
        end 
    end
%% 
    methods (Static)  
        function [lines,props] = Figure(grid,promptxt,inlines,props,isdel)
            %
            %-------function help------------------------------------------
            % PURPOSE
            %   Figure to interactively edit points or lines
            % USAGE
            %   lines = PL_CentreLine.setFigure(grid,promptxt,inlines,props,isdel);
            % INPUTS
            %   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
            %   promptxt- character string used for initial prompt in title
            %   inlines - can be a struct of points and lines or a struct
            %             (or table) of x,y vectors to be edited or the 
            %             format of output if no lines are being input -
            %             see Outputs for details
            %   props - struct for maximum water level and depth exponent (maxwl,dexp,cint)
            %   isdel - logical flag true to delete figure on completion - optional, 
            %           default is false
            % OUTPUTS
            %   lines - if lines are input, format is the same as the input, otherwise
            %           if nlines is scalar (default it type 2 if nlines empty)
            %            outype=0: array of structs with x, y and z fields defining selected points,
            %            outype=1: Nx2 or Nx3 array.
            %            outype=2: struct with x, y (and z) vector fields
            %            outype=3: table with x, y (and z) vector fields
            %            lines = [] if user closes figure, or no points defined
            %   props - as input with any updates
            % NOTES
            %   NB: if the figure window is closed the function returns lines=[] and
            %       not the input points
            %--------------------------------------------------------------
            %
            if nargin<5, isdel = false; end
            figtitle = sprintf('Extract centre-line');
            tag = 'PlotFig'; %used for collective deletes of a group
            position = [0,0.03,1,0.93];
            obj = PL_CentreLine(figtitle,tag,position);

            %plot grid and initialise axes (needed for context menus)
            initialiseGrid(obj,grid,props)
            obj.Axes.Title.String = promptxt; %initial prompt

            %define menu to be used
            mtext = {'Centreline','Line','Figure'};
            mcall = repmat({[]},1,length(mtext));
            menu = struct('label',mtext,'callback',mcall);
            %default menus are defined for Point, Line and Figure include:
            % Line: 'Add','Edit','Extend','Insert','Join','Split','Delete'
            % Figure: 'View','Undo','Save','Save & Exit','Quit'
            %and are defined in PLinterface. getCallBacks
            %Bespoke options are then added by creating tables for each
            %addtional menu option in the function setSubMenus        
            submenus = setSubMenus(obj);
            %set the menus and submenus            
            obj = setMenu(obj,menu,submenus);

            %handle input of existing lines
            outype = getInLines(obj,inlines);

            obj.Figure.Visible = 'on';
            obj = waitForFigure(obj);

            if isempty(obj.outLines)
                lines = [];  
            else                
                lines = gd_points2lines(obj.outLines,outype.lines);
            end
            props = obj.Props;

            %delete figure if isdel has been set by call.
            if isdel
                delete(obj.Figure)
            end
        end
    end

%--------------------------------------------------------------------------
% Additional callback functions
%--------------------------------------------------------------------------
    methods (Access=protected)
        function SetLine(obj,~,~)
            %add a centreline of a channel to the existing lines
            prompt1 = sprintf('Centre-line\nSelect first point of centre-line');   
            prompt2 = sprintf('Centre-line\nSelect second point of centre-line');
            prompts = {prompt1,prompt2};
            
            [points,Hp] = setPoints(obj,prompts,'endpoints');
            if ~isempty(points) 
                pline = getCentreLine(obj,points);
                if isempty(pline), return; end
                delete(Hp)                                                 %remove set points
                Hl = checkPlot(obj,pline);
                obj.pLines = setCentreLine(obj,pline);  
                delete(Hl)                                 %remove check line
            end

            clearGraphics(obj,{'endpoints','newcline','clines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);     %2= plot as lines (editable)
        end

%%
        function Resample(obj,~,~)
            %resample the lines at user specified interval 
            %NB overloads function in PLinterface to update Props
            cint = PLinterface.setInterval();
            obj.pLines = resampleLines(obj,cint);
            clearGraphics(obj,{'mylines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);     %2= plot as lines
            if isprop(obj,'Props')
                obj.Props.cint = cint;
            end
        end  

%%
%--------------------------------------------------------------------------
%  Utility functions to implement various actions
%-------------------------------------------------------------------------- 
        function pline = getCentreLine(obj,points)
            %create a centreline of a channel using function a_star to trace the
            %shortest path between start and end points whilst finding the deepest
            %points (ie a thalweg)
            grid = obj.Grid; props = obj.Props;
            %index of nearest grid point to selected start end end points    
            start = dsearchn(grid.xy,[points(1).x,points(1).y]); 
            goal = dsearchn(grid.xy,[points(2).x,points(2).y]);

            hwb = progressbar([],'Computing centre-line');
            %find the shortest path taking account of the cost (depths)
            %Z(Z>maxwl) = 0;
            costfnc = @(z) -(min(z,[],'all')-z).^props.dexp; %weighted inverse depth to favour staying deep
            thalweg = a_star(grid.water, costfnc(grid.Z), start, goal);
            [idy,idx] = ind2sub(size(grid.Z),thalweg); %convert indices to x,y subscripts
            progressbar(hwb);

            line.x = [flipud(grid.x(idx));NaN]; %return points in order from start point
            line.y = [flipud(grid.y(idy));NaN]; %as column vectors terminated with NaNs
            pline = gd_lines2points(line);      %return as a pline
        end

%%
        function plines = setCentreLine(obj,newline)
            %allow user to accept or reject the new addition to the centre-line
            props = obj.Props;
            plines = obj.pLines;
            nline = newline(1:end-1);                     %remove trailing NaNs
            nline = gd_points2lines(nline,1);
            ansr = questdlg('Accept the centreline?','Centre-line','Yes','No','Yes');
            if strcmp(ansr,'Yes')
                %convert output to specified point pacing
                clength = sum(vecnorm(diff(nline),2,2));  %nline is a column vector [Nx2]
                cpoints = round(clength/props.cint);      %number of points in new line
                newpoints = curvspace(nline,cpoints);     %curvespace uses [Nx2]
                newline = [newpoints;[NaN,NaN]];          %restore trailing NaNs 
                pline = gd_lines2points(newline);
                plines = [plines,pline];
            end           
        end

%%
        function H = checkPlot(obj,cline)
            %plot base map of initial grid selection 
            % cline is xy single stuct line
            clearGraphics(obj,{'mylines'});
            points = [cline(1).x,cline(1).y;cline(end-1).x,cline(end-1).y];
            hold on
            hp = plot(obj.Axes,points(:,1),points(:,2),'ok','MarkerSize',8,...
                                   'MarkerFaceColor','w','Tag','endpoints');
            hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            hl = plot(obj.Axes,[cline(:).x],[cline(:).y],'-r','LineWidth',2,...
                                           'DisplayName','newcline');
            if ~isempty(obj.pLines)
                gd_plotpoints(obj.Axes,obj.pLines,'clines',5);      %5= plot as centre-lines
            end
            hold off
            H = [hl;hp];
            title('Centre-line between defined start and end points')
        end

%%
        function initialiseGrid(obj,grid,props)
            %plot grid and add varaiables needed for extracting centre-lines
            [X,Y] = meshgrid(grid.x,grid.y);
            N = numel(X);
            grid.xy = [reshape(X,[N,1]),reshape(Y,[N,1])];
            grid.Z = grid.z';
            %accessible map (water) and use -Z as the cost map
            water = true(size(grid.Z));
            water(isnan(grid.Z) | grid.Z>props.maxwl) = false;
            grid.water = water;
            obj.Grid= grid;       
            obj.Grid.z(~water') = NaN;
            obj.Props = props;
            obj.Axes = gd_plotgrid(obj.Figure,obj.Grid);
        end 

%%
        function calltable = setSubMenus(obj)  
            %user defined menus to be appended to the Figure menu 
            varnames = {'Parent','Callback','Label'};
            %Boundary menu variables
            stext = ["Set";"Resample";"Smooth";"Reset"];
            scall = {@obj.SetLine; @obj.Resample; @obj.Smooth; @obj.Reset}; 
            nrec = length(stext);
            spart = repmat("Centreline",nrec,1);
            calltable = table(spart,scall,stext,'VariableNames',varnames);
        end
    end
%--------------------------------------------------------------------------
% Static utility functions
%--------------------------------------------------------------------------
    methods (Static, Access=protected)
        %Static methods in PLinterface
        % checkDirection
        % isPointNearLine
        % setLevel
        % setInterval
    end
end