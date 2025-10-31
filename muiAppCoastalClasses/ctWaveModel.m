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
            [tsdst,meta] = getInputData(obj,mobj);            
            if isempty(tsdst), return; end   %user cancelled data selection
            setRunParam(obj,mobj,meta.caserecs{:}) %assign run parameters
            dsp = modelDSproperties(obj,isin,meta.iselvar);

            %get the input site parameters as a class instance
            site_params = mobj.Inputs.ctWaveParameters;  
            
            if isin
                inp = inputParameters(site_params); %convert class to struct
                inp.g = mobj.Constants.Gravity;     %add gravity
                [Hsi,Diri,depi,bs] = hs_surf(tsdst,inp);      
                if meta.iselvar %variables selected (non-standard names) so include Tp
                    results = {Hsi,tsdst.Tp,Diri,tsdst.swl,depi};
                else    %default naming convention
                    results = {Hsi,Diri,tsdst.swl,depi};
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
        % function [tsdst,caserec] = getWaveModelDataset(obj,mobj,type,addnames,caserec)
        %     %prompt user to select a model wave dataset and add Tp if inshore
        %     %
        %     muicat = mobj.Cases;
        %     if nargin<4
        %         addnames = {'Tp'};  %default is to add Tp
        %     end
        %     %
        %     if nargin<5   %no caserec to prompt for selection
        %         [wvobj,wvdst,ok] = selectClassInstance(obj,'ModelType',type);
        %         if ok<1, tsdst = []; caserec = []; return; end
        %         caserec = caseRec(muicat,wvobj.CaseIndex);
        %     else
        %         wvobj = getCase(muicat,caserec);
        %         wvdst = wvobj.Data.Dataset;
        %     end            
        % 
        %     %if inshore wave dataset add variables requested otherwise just
        %     %copy dstable
        %     tsdst = copy(wvdst);
        %     if strcmp(type,'Inwave_model')
        %         inpwavecid = wvobj.RunParam.ctWaveData.caseid;  %source dataset (offshore)
        %         inpwaverec = caseRec(muicat,inpwavecid);        %case record
        %         inpdst = getDataset(muicat,inpwaverec,1);       %dst used to create inshore waves
        %         dstnames = inpdst.VariableNames;                %source variables
        %         varnames = tsdst.VariableNames;                 %inshore wave variables
        %         for i=1:length(addnames)    
        %             if any(strcmp(varnames,addnames{i}))
        %                 %variable already included
        %                 continue;
        %             elseif any(strcmp(dstnames,addnames{i})) && ...
        %                                     ~isempty(inpdst.(addnames{i}))
        %                 %variable to be added exists in source dataset
        %                 tsdst = addvars(tsdst,inpdst.(addnames{i}),...
        %                                    'NewVariableNames',addnames{i});
        %             else
        %                 warndlg(sprintf('Variable %s not found so not added to wave dataset',...
        %                                                    addnames{i}));
        %             end
        %         end
        %     end
        % end
%%        
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 
%%    
    methods (Access = private)
        % function [tsdst,inputxt,isel] = getInputData(obj,mobj)
        %     %prompt user to select wave and water level data and return in
        %     %input dstable of data and metadata for inputs used
        %     tsdst = []; inputxt = []; isel = false;
        %     muicat = mobj.Cases;
        %     promptxt = 'Select input wave data set:';           
        %     [wv_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
        %                    'CaseClass',{'ctWaveData','muiUserModel'},'ListSize',[300,100]);                                    
        %     if ok<1, return; end
        %     wvdst = getDataset(muicat,wv_crec,1);
        %     wvtime = wvdst.RowNames;
        %     inputxt = sprintf('%s used for offshore waves',wvdst.Description);
        % 
        %     %check whether default variable names are not used and selection needed
        %     varnames = wvdst.VariableNames;
        %     if ~any(strcmp(varnames,'Hs'))
        %         wvdst = extract_wave_data(wvdst);
        %         if isempty(wvdst), return; end
        %         isel = true;
        %     end
        % 
        %     promptxt = 'Select input water level data set (Cancel to use SWL=0):';           
        %     [wl_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
        %                         'CaseClass',{'ctWaterLevelData','ctTidalAnalysis',...
        %                                       'muiUserModel'},...
        %                         'ListSize',[300,100]); 
        % 
        %     swl = zeros(size(wvtime));               
        %     if ok<1 || isempty(wl_crec)
        %         getdialog('Using SWL=0');
        %         inputxt = sprintf('%s, 0mOD used for water level',inputxt);
        %         wl_crec = 0;     %assign a null value if no water level data available
        %     else
        %         wldst = getDataset(muicat,wl_crec,1); 
        %         %check that there is water level data for period of interest
        %         [idst,idnd] = ts2_endpoints_in_ts1(wvdst,wldst);
        %         if isempty(idst)
        %             getdialog('Data do not overlap. Using SWL=0');
        %             inputxt = sprintf('%s, 0mOD used for water level',inputxt);
        %             wl_crec = 0; %assign a null value if no water level data available
        %         else 
        %             %select a variable from the water level dataset
        %             varnames = wldst.VariableNames;
        %             idx = 1;
        %             if length(varnames)>1
        %                 [idx,ok] = listdlg('Name','WL options', ...
        %                     'PromptString','Select a variable:', ...
        %                     'SelectionMode','single','ListSize',[200,100],...
        %                     'ListString',varnames);
        %                 if ok<1, idx = 1; end
        %             end
        %             wldata = wldst.(varnames{idx});
        %             wltime = wldst.RowNames;
        %             swltime = wvtime(idst:idnd);                                        
        %             %now interpolate water levels onto wave height times
        %             swl(idst:idnd,1) = interp1(wltime,wldata,swltime,'linear','extrap');
        %             swl(isnan(swl)) = 0;
        %             inputxt = sprintf('%s, %s used for water levels',...
        %                                         inputxt,wldst.Description);
        %         end
        %     end
        %     tsdst = addvars(wvdst,swl,'NewVariableNames','swl');
        %     %assign the run parameters to the model instance
        %     if wl_crec==0
        %         setRunParam(obj,mobj,wv_crec);
        %     else
        %         setRunParam(obj,mobj,wv_crec,wl_crec); %input caserecs passed as varargin     
        %     end
        % end
%%
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
            if isin && isel                               %inshore waves
                dsp.Variables = struct(...                      
                    'Name',{'Hsi','Tp','Diri','swl','depi'},...
                    'Description',{'Inshore wave height','Peak wave period',...
                                   'Inshore wave direction',...
                                   'Still water level','Inshore depth'},...
                    'Unit',{'m','s','deg','mOD','m'},...
                    'Label',{'Wave height (m)','Wave period (s)','Wave direction (deg)',...
                               'Water level (mOD)','Water depth (m)'},...
                    'QCflag',repmat({'model'},1,5)); 
            elseif isin                                  %inshore waves
                dsp.Variables = struct(...                      
                    'Name',{'Hsi','Diri','swl','depi'},...
                    'Description',{'Inshore wave height','Inshore wave direction',...
                               'Still water level','Inshore depth'},...
                    'Unit',{'m','deg','mOD','m'},...
                    'Label',{'Wave height (m)','Wave direction (deg)',...
                               'Water level (mOD)','Water depth (m)'},...
                    'QCflag',repmat({'model'},1,4)); 
            else                                         %deepwater waves
                dsp.Variables = struct(...                      
                    'Name',{'Hso','Tp','Diro'},...
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