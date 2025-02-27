classdef PL_PlotSections < PLinterface
%
%-------class help---------------------------------------------------------
% NAME
%   PL_SectionLines.m
% PURPOSE
%   Class to create to define and plot cross-sections on a grid
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
        Plots          %handles to section plots
        nsections = 0  %number of sections drawn
    end
       
    methods
        function obj = PL_PlotSections(figtitle,tag,position)          
            %constructor code: 
            if nargin<3, position = []; end
            obj = setFigure(obj,figtitle,tag,position);
        end 
    end
%% 
    methods (Static)  
        function Figure(grid,promptxt,isdel)
            %
            %-------function help------------------------------------------
            % PURPOSE
            %   Figure to interactively edit points or lines
            % USAGE
            %   lines = PL_SectionLines.Figure(grid,promptxt,clines,inlines,isdel);
            % INPUTS
            %   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
            %   cobj - instance of EDBimport class with grid or geoimages
            %   promptxt- character string used for initial prompt in title
            %   inlines - can be a struct of points and lines or a struct
            %             (or table) of x,y vectors to be edited or the 
            %             format of output if no lines are being input -
            %             see Outputs for details
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
            %   All lines are terminated with a NaN (as per shapefiles),
            %   whereas points are not.
            %--------------------------------------------------------------
            %
            if nargin<3, isdel = false; end
            figtitle = sprintf('Plot Sections');
            tag = 'PlotFig'; %used for collective deletes of a group
            position = [0,0.03,1,0.93];
            obj = PL_PlotSections(figtitle,tag,position);

            %plot grid and initialise axes (needed for context menus)
            obj.Axes = gd_plotgrid(obj.Figure,grid);          
            obj.Axes.Title.String = promptxt; %initial prompt
            obj.Grid = grid; %stored for use in callback functions

            %define menu to be used
            mtext = {'Define sections','Plot sections','Figure'};
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
            
            obj.Figure.Visible = 'on';
            obj = waitForFigure(obj);

            %delete figure if isdel has been set by call.
            if isdel
                delete(obj.Figure)
            end
        end
    end

%--------------------------------------------------------------------------
% Additional Menu callback functions
%--------------------------------------------------------------------------
    methods (Access=protected)
        function setSections(obj,~,~)
            %use centreline to generate a set of sections
            resetMenu(obj)
            prompt1 = sprintf('Define section\nSelect first point of section');   
            prompt2 = sprintf('Define section\nSelect second point of section');
            prompts = {prompt1,prompt2};

            [newpnts,Hp] = setPoints(obj,prompts,'spoints');
            while ~isempty(newpnts)     
                NewLine =  [newpnts,obj.NaNpnts];
                obj.pLines = [obj.pLines,NewLine];
                delete(Hp)      
                obj.nsections = obj.nsections+1;
                numtxt = num2str(obj.nsections);
                gd_plotpoints(obj.Axes,NewLine,'mylines',2,obj.isXYZ); %plot as line
                gd_plotpoints(obj.Axes,NewLine,numtxt,3); %set labels
                [newpnts,Hp] = setPoints(obj,prompts,'spoints');
            end
            % clearGraphics(obj,{'mylines','mytext'});
            % obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);  %set line  
            resetMenu(obj,false)  
        end

%%
        function deleteLine(obj,~,~)
            %order the lines to be in along-channel sequence
            resetMenu(obj);
            promptxt = sprintf('Delete line\nSelect line to Delete, right click on any line to quit\nRedraw to update centre-line');
            [deline,H] = gd_getpline(obj.Axes,promptxt,'mylines');         %get line to delete
            while ~isempty(deline)                
                [obj.pLines,~] = deleteAline(obj,deline);                  %delete the line
                delete(H)
                obj.nsections = obj.nsections-1;
                labelLines(obj,obj.pLines);                                %update line numbers
                resetLines(obj);
                deline = gd_getpline(obj.Axes,promptxt,'mylines');         %get line to delete
            end
            resetLines(obj)
            resetMenu(obj,false)       
        end

%%
        function Reset(obj,~,~)
            %Clear all the current section lines
            clearGraphics(obj,{'mylines','mytext'});
            obj.pLines = [];
            obj.nsections = 0;
        end

%%
        function plotSections(obj,~,~)
            %plot the current set of sections
            idx = length(obj.Plots)+1;
            if isempty(obj.pLines), return; end
            cplines = gd_plines2cplines(obj.pLines);
            [~,hf] = gd_plotsections(obj.Grid,cplines);
            hf.CloseRequestFcn = @obj.closePlot;
            obj.Plots(idx) = hf.Number;
        end

%%
        function clearPlots(obj,~,~)
            %delete all section plots
            delete(obj.Plots)
            obj.Plots = [];
        end

%%
        function Redraw(obj,~,~)
            %redraw the current set of points and lines
            ax = obj.Axes;            
            clearGraphics(obj,{'mylines'});

            if ~isempty(obj.pLines)   
                gd_plotpoints(ax,obj.pLines,'mylines',2); %set line          
            end
        end  
%%
%--------------------------------------------------------------------------
%  Utility functions to implement various actions
%-------------------------------------------------------------------------- 
        function closePlot(obj,src,~)
            %close figure and update index
            idx = obj.Plots==src.Number;
            obj.Plots(idx) = [];
            delete(src)
        end
        
%%
        function calltable = getCallBacks(obj,submenus)
            %get the default menus for Point,Line and Figure and append
            %any user defined menus
            % overload function from PLinterface to limit menu optiona
            varnames = {'Parent','Callback','Label'};
            
            %default Figure menu variables
            stext = ["Redraw";"Exit"];
            scall = {@obj.Redraw; @obj.Exit}; 
            nrec = length(stext);
            spart = repmat("Figure",nrec,1);
            calltable = table(spart,scall,stext,'VariableNames',varnames);

            %add any user defined menu options
            calltable = [calltable;submenus];
        end

%%
        function calltable = setSubMenus(obj)  
            %user defined menus to be appended to the Figure menu 
            varnames = {'Parent','Callback','Label'};
            %Set sections menu variables
            stext = ["Define";"Edit";"Delete";"Clear all"];
            scall = {@obj.setSections;  @obj.editLine; @obj.deleteLine; ...
                                                    @obj.Reset}; 
            nrec = length(stext);
            spart = repmat("Define sections",nrec,1);
            calltable = table(spart,scall,stext,'VariableNames',varnames);

            %Plot sections menu variables 
            stext =  ["Plot";"Clear all"];                                        
            scall = {@obj.plotSections; @obj.clearPlots}; 
            nrec = length(stext);
            spart = repmat("Plot sections",nrec,1);                                
            call{1} = table(spart,scall,stext,'VariableNames',varnames);

            for i=1:length(call)
                calltable = [calltable;call{i}]; %#ok<AGROW> 
            end
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