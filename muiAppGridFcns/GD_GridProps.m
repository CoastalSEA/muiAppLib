classdef GD_GridProps < muiPropertyUI    
%
%-------class help---------------------------------------------------------
% NAME
%   GD_GridProps.m
% PURPOSE
%   Class for grid and time step parameters used in ChannelForm model
% USAGE
%   obj = GD_GridProps.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'X-axis definition [x0 xN]',...
                          'No. of intervals in the x direction',...
                          'Y-axis definition [y0 yN]',...
                          'No. of intervals in the y direction',...
                          'Vertical resolution for hypsometry'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        XaxisLimits = [0 100000]   %definition of the x-axis [x0 xN]         
        Xint = int32(500)          %no of intervals in the x direction
        YaxisLimits = [-5000 5000] %definition of the y-axis [y0 yN]
        Yint = int32(500)          %no of intervals in the y direction
        histint = 0.1;             %vertical resolution for hypsometry histogram
    end    

%%   
    methods (Access=protected)
        function obj = GD_GridProps(mobj)          
            %constructor code:            
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
            
            %to use non-numeric entries then one can either pre-assign 
            %the values in the class properties defintion, above, or 
            %specify the PropertyType as a cell array here in the class 
            %constructor, e.g.:
            % obj.PropertyType = [{'datetime','string','logical'},...
            %                                       repmat({'double'},1,8)];
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'GD_GridProps';        
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = GD_GridProps(mobj);         
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end   
%%
        function obj = getGridProps(mobj)
            %enable external function to access GD_GridProps without
            %setting the properties (see setInput above)
             obj = GD_GridProps(mobj);  
        end
    end
%%        
	methods     
        function [x,y,delx,dely] = getGridDimensions(obj)
            %define x and y dimensions based on GridProps settings
            Lx = diff(obj.XaxisLimits);   %length of model domain (m)        
            Ly = diff(obj.YaxisLimits);   %width of model domain (m)  
            %define grid intervals
            delx = double(Lx/obj.Xint);            
            dely = double(Ly/obj.Yint);
            %set grid dimensions
            x = (obj.XaxisLimits(1):delx:obj.XaxisLimits(2))';
            y = (obj.YaxisLimits(1):dely:obj.YaxisLimits(2))';
        end
%%
        function obj = setGridDimensions(obj,x,y)
            %use x and y dimensions to update the grid parameters
            obj.XaxisLimits = [x(1),x(end)];  %set grid limits
            obj.YaxisLimits = [y(1),y(end)];
            %get grid intervals
            Lx = abs(diff(obj.XaxisLimits));   %length of model domain (m)        
            Ly = abs(diff(obj.YaxisLimits));   %width of model domain (m)
            delx = abs(x(2)-x(1));          
            dely = abs(y(2)-y(1));
            %set number of grid intervals
            obj.Xint = Lx/delx;
            obj.Yint = Ly/dely;            
        end
    end
%%
    methods 
        function obj = setGridProperties(obj)
            %add nrec to limit length of props UI (default=12)                           
            obj = editProperties(obj);  
        end
    end
end