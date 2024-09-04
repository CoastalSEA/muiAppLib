classdef ctStructureInput < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   ctStructureInput.m
% PURPOSE
%   Class to handle input data required by overtopping model
% SEE ALSO
%   otop_Q.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2017
%--------------------------------------------------------------------------
%     
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Crest level (mOD)','Crest width (m)',...
                          'Upper wall slope (1:uws)','Berm level (mOD)',...
                          'Berm width (m)','Lower wall slope (1:lws)',...
                          'Toe level (mOD)','Wall roughness (-)'};
        %abstract properties in PropertyInterface for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        CrestLevel = 0;           %structure crest level (mOD)
        CrestWidth = 0;           %structure crest width (m)
        UpperWallSlope = 1;       %upper wall slope (1:uws)
        BermLevel = 0;            %berm level (mOD)
        BermWidth = 0;            %berm width (m)
        LowerWallSlope = 1;       %lower wall slope (1:lws)
        ToeLevel = 0;             %struture toe level (mOD)
        WallRoughness = 1;        %roughness coeficient (-)
    end    

%%   
    methods (Access=protected)
        function obj = ctStructureInput(mobj) 
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
            classname = 'ctStructureInput'; 
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = ctStructureInput(mobj);             
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
        function otopstruct = getStructure(obj)
            %get the structure properties as a struct
            ot.cl  = obj.CrestLevel;       %structure crest level (mOD)
            ot.cw  = obj.CrestWidth;       %structure crest width (m)
            ot.uws = obj.UpperWallSlope;   %upper wall slope (1:uws)
            ot.bl  = obj.BermLevel;        %berm level (mOD)
            ot.bw  = obj.BermWidth;        %berm width (m)
            ot.lws = obj.LowerWallSlope;   %lower wall slope (1:lws)
            ot.tl  = obj.ToeLevel;         %struture toe level (mOD)
            ot.r   = obj.WallRoughness;    %roughness coeficient (-)
            otopstruct = ot;
        end
    end
end