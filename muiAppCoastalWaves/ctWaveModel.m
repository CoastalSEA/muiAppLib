classdef ctWaveModel < waveModels
%
%-------class help------------------------------------------------------
% NAME
%   ctWaveModel.m
% PURPOSE
%    Class for inshore and offshore wave models to be run in CoastalTools 
%    and other ModelUI apps. Inherits abstract class waveModels which is
%    subclass of muiDataSet.
% SEE ALSO
%   muiDataSet, waveModels
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
        InshoreBedSlope   %bed slope as single value or vector array
    end
    
    properties (Hidden)
        ModelType        %model used for the particular instance
    end
    
    methods (Access={?muiDataSet,?muiStats,?ctWaveData})
        function obj = ctWaveModel()                    
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj,isin)
            %function to run the wave refraction model. isin is true for
            %inshore waves and false for deepwater waves
           
            if nargin<2
                %determine if inshore or offshore to be used                
                sel = questdlg('Use offshore or inshore model?','Wave model',...
                               'Offshore','Inshore','Inshore');
                if strcmp(sel,'Inshore')
                    isin = true;
                else
                    isin = false;
                end
            end
            
            obj = ctWaveModel;                          

            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in CoastalTools
            %NB - water level data not checked because optional
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            %get the timeseries input data and site parameters
            [tsdst,meta] = obj.getInputData(mobj);            
            if isempty(tsdst), obj = []; return; end   %user cancelled data selection
            setRunParam(obj,mobj,meta.caserecs{:}) %assign run parameters
            dsp = modelDSproperties(obj,isin,meta.iselvar);

            %get the input site parameters as a class instance
            site_params = mobj.Inputs.ctWaveParameters;  
            
            if isin
                inp = inputParameters(site_params); %convert class to struct
                inp.g = mobj.Constants.Gravity;     %add gravity
                if numel(tsdst)>1
                    warndlg('Plane bed refraction does not work for compound sea states')
                    obj = []; return;
                end
                [Hs,Dir,depi,bs] = hs_surf(tsdst,inp);      
                if isempty(Hs)
                    obj = []; return
                elseif meta.iselvar %variables selected (non-standard names) so include Tp
                    results = {Hs,tsdst.Tp,Dir,tsdst.swl,depi};
                else    %default naming convention
                    results = {Hs,Dir,tsdst.swl,depi};
                end
                %inshore bed level can be array or single value so not added to table           
                obj.InshoreBedSlope = bs; 
                obj.ModelType = 'Inwave_model';
            else
                results = callRefraction(obj,tsdst,site_params);
                obj.ModelType = 'Offwave_model';
            end
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            time = tsdst.RowNames;
            dst = dstable(results{:},'RowNames',time,'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model            
            dst.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                         obj.ModelType);
            dst.MetaData = meta.inptxt;            
            %save results
            setDataSetRecord(obj,mobj.Cases,dst,'model');
            getdialog('Run complete');
        end
    end

%%
    methods     
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 

%%    
    methods (Access = private)
        function results = callRefraction(~,tsdst,site_params)
            %parse inputs for call to refraction function
            z0 = site_params.OffshoreBedLevel;
            theta0 = site_params.OffshoreAngle;
            Kf = site_params.FrictionCoefficient;
            
            prompt = {'Angle of bed contour at refraction point(s) (degTN)',...
                      'Bed level at refraction point(s) (mOD)'};
            default = {num2str(theta0),num2str(-100)};
            answer = inputdlg(prompt,'Deepwater refraction',1,default);
            if isempty(answer), return; end            
            thetaT = str2double(answer{1});
            zT = str2double(answer{2});
            dep0 = tsdst.swl-z0;   %source water depth from swl to bed level 
            depT = tsdst.swl-zT;   %target water depth from swl to bed level 
            
            %this is a general refraction model so 'isshore' is false
            [Hso,Diro] = refraction(tsdst.Hs,tsdst.Tp,tsdst.Dir,[dep0,depT],...
                                                [theta0,thetaT],Kf,false); 
            results = {Hso,tsdst.Tp,Diro};
        end
        
%%
        function dsp = modelDSproperties(~,isin,isel) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            if isin && isel                               %inshore waves & non-default variables
                dsp.Variables = struct(...                      
                    'Name',{'Hs','Tp','Dir','swl','depi'},...
                    'Description',{'Inshore wave height','Peak wave period',...
                                   'Inshore wave direction',...
                                   'Still water level','Inshore depth'},...
                    'Unit',{'m','s','deg','mOD','m'},...
                    'Label',{'Wave height (m)','Wave period (s)','Wave direction (deg)',...
                               'Water level (mOD)','Water depth (m)'},...
                    'QCflag',repmat({'model'},1,5)); 
            elseif isin                                  %inshore waves
                dsp.Variables = struct(...                      
                    'Name',{'Hs','Dir','swl','depi'},...
                    'Description',{'Inshore wave height','Inshore wave direction',...
                               'Still water level','Inshore depth'},...
                    'Unit',{'m','deg','mOD','m'},...
                    'Label',{'Wave height (m)','Wave direction (deg)',...
                               'Water level (mOD)','Water depth (m)'},...
                    'QCflag',repmat({'model'},1,4)); 
            else                                         %deepwater waves
                dsp.Variables = struct(...                      
                    'Name',{'Hs','Tp','Dir'},...
                    'Description',{'Deepwater wave height','Peak period','Deepwater wave direction'},...
                    'Unit',{'m','s','deg'},...
                    'Label',{'Wave height (m)','Peak wave period (s)',...
                                            'Wave direction (deg)'},...
                    'QCflag',repmat({'model'},1,3));                 
            end
            %
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dsp.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});  
        end
        
    end    
end    