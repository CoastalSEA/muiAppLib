classdef PL_Editor < PLinterface
%
%-------class help---------------------------------------------------------
% NAME
%   PL_Editor.m
% PURPOSE
%   Class to edit points and lines using the methods in PLintereface
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
        outPoints = []
        outLines = []
        isXYZ = false
    end
       
    methods
        function obj = PL_Editor(figtitle,tag,position)          
            %constructor code: 
            if nargin<3, position = []; end
            obj = setFigure(obj,figtitle,tag,position);
        end 
    end
%% 
    methods (Static)  
        function [lines,points] = Figure(grid,promptxt,inlines,isxyz,isdel)
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
            %   isxyz - logical flag true to input z values
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
            %   points - same format options as lines. vectors define xy(z)
            %            points (no NaN separators)
            % NOTES
            %   NB: if the figure window is closed the function returns lines=[] and
            %       not the input points
            %   All lines are terminated with a NaN (as per shapefiles),
            %   whereas points are not.
            %--------------------------------------------------------------
            %
            if nargin<5, isdel = false; end
            figtitle = sprintf('Edit lines');
            tag = 'PlotFig'; %used for collective deletes of a group
            position = [0,0.03,1,0.93];
            obj = PL_Editor(figtitle,tag,position);
            obj.isXYZ = isxyz; 
            %plot grid and initialise axes (needed for context menus)
            obj.Axes = gd_plotgrid(obj.Figure,grid);
            obj.Axes.Title.String = promptxt; %initial prompt

            %define menu to be used
            mtext = {'Line','Point','Figure'};
            mcall = repmat({[]},1,length(mtext));
            menu = struct('label',mtext,'callback',mcall);
            %default menus are defined for Point, Line and Figure include:
            % Point: 'Add','Edit','Delete'
            % Line: 'Add','Edit','Extend','Insert','Join','Split','Delete'
            % Figure: 'View','Undo','Save','Save & Exit','Quit'
            %and are defined in PLinterface. getCallBacks
            %
            %Bespoke options are then added by creating tables for each
            %addtional menu option in the function setSubMenus        
            submenus = setSubMenus(obj);
            %set the menus and submenus            
            obj = setMenu(obj,menu,submenus);

            %handle input of existing lines
            outype = getInLines(obj,inlines);

            obj.Figure.Visible = 'on';
            obj = waitForFigure(obj);

            %convert format of output if required
            if isempty(obj.outPoints)
                points = [];     
            else
                points = gd_points2lines(obj.outPoints,outype.points);           
            end

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

%%
    methods (Access=protected)
        function calltable = setSubMenus(obj) %#ok<MANU> 
            calltable = [];
            %see PLinterface for sample setSubMenus functions
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
        % findNearestPoint
        % setLevel
        % setInterval
    end
end
