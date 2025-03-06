classdef PL_Boundary < PLinterface
%
%-------class help---------------------------------------------------------
% NAME
%   PL_Boundary.m
% PURPOSE
%   Class to extract contours and generate model boundaries
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
    end
       
    methods
        function obj = PL_Boundary(figtitle,tag,position)          
            %constructor code: 
            if nargin<3, position = []; end
            obj = setFigure(obj,figtitle,tag,position);
        end 
    end
%% 
    methods (Static)  
        function lines = Figure(grid,promptxt,inlines,isdel)
            %
            %-------function help------------------------------------------
            % PURPOSE
            %   Figure to interactively edit points or lines
            % USAGE
            %   lines = PL_Editor.setFigure(grid,promptxt,inlines,isdel);
            % INPUTS
            %   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
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
            % NOTES
            %   NB: if the figure window is closed the function returns lines=[] and
            %       not the input points
            %--------------------------------------------------------------
            %
            if nargin<4, isdel = false; end
            figtitle = sprintf('Extract boundary');
            tag = 'PlotFig'; %used for collective deletes of a group
            position = [0,0.03,1,0.93];
            obj = PL_Boundary(figtitle,tag,position);

            %plot grid and initialise axes (needed for context menus)
            obj.Axes = gd_plotgrid(obj.Figure,grid);
            obj.Axes.Title.String = promptxt; %initial prompt

            %define menu to be used
            mtext = {'Boundary','Line','Figure'};
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
            %making figure visible after call to gd_getcontour seems to
            %cause a problem so set visible moved before
            obj.Figure.Visible = 'on';

            %handle input of existing lines
            outype = getInLines(obj,inlines);
            obj.Grid = grid; %assign as property for use in callbacks

            if isempty(obj.pLines)   %no lines imported
                zlevel = PLinterface.setLevel();
                if isempty(zlevel)
                    lines = []; delete(obj.Figure); return; 
                end                
                blines =  gd_getcontour(grid,zlevel,false);
                obj.pLines = gd_lines2points(blines);
                obj.outLines = obj.pLines; %?????????check that this suits workflow
                obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2); %2= plot as lines
            end
            
            %wait for user to close the PL figure
            obj = waitForFigure(obj);

            if isempty(obj.outLines)
                lines = [];  
            else                
                lines = gd_points2lines(obj.outLines,outype.lines);
            end

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
        function Reset(obj,~,~)
            %resample contour to create new boundary
            zlevel = PLinterface.setLevel();
            if isempty(zlevel), return; end
            blines =  gd_getcontour(obj.Grid,zlevel,false);
            obj.pLines = gd_lines2points(blines);
            clearGraphics(obj,{'mylines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);     %2= plot as lines
        end

%%
        function calltable = setSubMenus(obj)  
            %user defined menus to be appended to the Figure menu 
            varnames = {'Parent','Callback','Label'};
            %Boundary menu variables
            stext = ["Resample";"Smooth";"Reset"];
            scall = {@obj.Resample; @obj.Smooth; @obj.Reset}; 
            nrec = length(stext);
            spart = repmat("Boundary",nrec,1);
            calltable = table(spart,scall,stext,'VariableNames',varnames);
        end
    end
%--------------------------------------------------------------------------
% Static utility functions
%--------------------------------------------------------------------------
    methods (Static, Access=protected)
        %Static methods in PLinterface
        % checkDirection
        % lineLength
        % isPointNearLine
        % setLevel
        % setInterval
    end
end