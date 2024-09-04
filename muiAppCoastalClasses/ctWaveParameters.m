classdef ctWaveParameters < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   ctWaveParameters.m
% PURPOSE
%   Class to handle input data required by InWaveModel
% SEE ALSO
%   ctWaveData.m and ctWaveModel.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2017
%--------------------------------------------------------------------------
%  
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Offshore bed level (mOD - so probably negative) ',...
                          'Offshore contour angle (degTN)',...
                          'Friction coefficient (value 0-1)',...
                          'Shoreline angle (degTN)',...
                          'Beach crest or HW level (mOD)',...
                          'Upper beach slope (1:s - enter value for s)',...
                          'Bed level 1km out from SWL (mOD) [or y,z]',...
                          'Nearshore bed level (NaN for surf zone depth)',... 
                          'Breaker height model for Hb (0,1,2, or 3 - see manual)',...
                          'Wave breaking model for Hsb (0,1,2, or 3 - see manual)',...
                          'Sediment grain size (m)',...
                          'CERC Drift coefficient, Kc'}
        %abstract properties in PropertyInterface for tab display (not used)
        TabDisplay   %structure defines how the property table is displayed 
    end    
    
    properties
        OffshoreBedLevel           %offshore bed level (OD-water depth) (mOD)
        OffshoreAngle              %bed contour angle at offshore point (degTN)
        FrictionCoefficient        %friction coefficient (default=1)
        ShorelineAngle             %angle of contours from north (degrees TN)
        BeachCrestLevel            %beach crest level (mOD)
        UpperBeachSlope            %bed slope (1:bs)
        BedLevelat1km              %bed level 1km out from SWL (mOD)
        InshoreBedLevel = NaN      %inshore bed level (mOD) - overrides use of surf depth
        Hb2Hsoption = 1            %flag to indicate method for estimating breaker height
                                   %0=0.78Hsi; 1=SPM on a slope; 2=SPM in front of toe; 3= as 2 for overtopping
        WaveBreakingModel = 1      %flag to indicate which wave breaking model to use
                                   %0=no breaking, 1=SPM breaking on a slope, 
                                   %2=Le Roux, 3=Hsb (Tucker/Holmes)
        GrainSize = 0.001          %sediment grains size (m)
        DriftCoefficient = 0.0006  %only used in CERC drift calcs (default=0.0006)       
    end
%%   
    methods (Access=protected)
        function obj = ctWaveParameters(mobj)       
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);     %muiPropertyUI function
            
            %to use non-numeric entries then one can either pre-assign 
            %the values in the class properties defintion, above, or specity the 
            %PropertyType as a cell array, e.g.:
            % obj.PropertyType = [{'datetime','string','logical'},...
            %                                       repmat({'double'},1,8)];
        end 
    end  
%%
    methods (Static)        
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'ctWaveParameters'; 
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = ctWaveParameters(mobj);             
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
            inp.z0 = obj.OffshoreBedLevel;    %offshore bed level (OD-water depth) (mOD)
            inp.offtheta = obj.OffshoreAngle; %bed contour angle at offshore point (degTN)
            inp.Kf = obj.FrictionCoefficient; %friction coefficient (default=1)
            inp.intheta = obj.ShorelineAngle; %angle of contours from north (degrees TN)
            inp.zBC = obj.BeachCrestLevel;    %beach crest level (mOD)
            inp.ubs = obj.UpperBeachSlope;    %bed slope (1:bs)
            inp.z1km = obj.BedLevelat1km;     %bed level 1km out from SWL (mOD)
            inp.zi = obj.InshoreBedLevel;     %inshore bed level (mOD) - overrides use of surf depth
            inp.hboption = obj.Hb2Hsoption;      %flag to indicate method for estimating breaker height
            inp.hsbflag = obj.WaveBreakingModel; %flag to indicate which wave breaking model to use
            inp.d50 = obj.GrainSize;          %sediment grains size (m)
            inp.Kc = obj.DriftCoefficient;    %only used in CERC drift calcs (default=0.0006)
        end   
    end
end    
    