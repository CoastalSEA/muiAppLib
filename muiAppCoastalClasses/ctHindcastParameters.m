classdef ctHindcastParameters < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   ctHindcastParameters.m
% PURPOSE
%   Class to handle input data required by WindWaveModel
% SEE ALSO
%   ctWindData.m and WindWaveModel.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2017
%--------------------------------------------------------------------------
%   
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Direction exponent (2-10)',...
                          'Fetch method (0=Donelan, 1=SPM)',...
                          'Height of wind speed measurement (m)',...
                          'Average water depth over fetch (m)',...
                          'Water depth at site (m)'};
        %abstract properties in PropertyInterface for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
       DirExp = 2             %Direction exponent
       FetchMethod = 0        %Fetch method (0=Donelan, 1=SPM)
       WindSpeedElev = 10     %Height of wind speed measurement (m)'
       AvFetchDepth           %Average water depth over fetch (m)'
       SiteDepth              %Water depth at site (m)
    end    

%%   
    methods (Access=protected)
        function obj = ctHindcastParameters(mobj)  
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI fcn
            
            %to use non-numeric entries then one can either pre-assign 
            %the values in properties defintion, above, or specity the 
            %PropertyType as a cell array, e.g.:
            % obj.PropertyType = [{'datetime','string','logical'},...
            %                                       repmat({'double'},1,8)];
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'ctHindcastParameters';          
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = ctHindcastParameters(mobj);             
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
    end
%%        
        %add other functions to operate on properties as required
    methods    
        function inp = inputParameters(obj)
            %construct inp struct used in models to define site properties
           inp.n = obj.DirExp;           %Direction exponent
           inp.Fopt = obj.FetchMethod;   %Fetch method (0=Donelan, 1=SPM)
           inp.zw = obj.WindSpeedElev;   %Height of wind speed measurement (m)'
           inp.df = obj.AvFetchDepth;    %Average water depth over fetch (m)'
           inp.ds = obj.SiteDepth;       %Water depth at site (m)
        end   
    end    
end