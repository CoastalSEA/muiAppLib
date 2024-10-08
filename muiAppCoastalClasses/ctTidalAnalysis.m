classdef ctTidalAnalysis < muiDataSet
%
%-------class help---------------------------------------------------------
% NAME
%   TidalAnalysis.m
% PURPOSE
%   Class to hold tidal analysis methods and output. Implements calling of
%   utide to to analyses water level data and store the constituents and to
%   use the constituents to generate a tidal record.
% NOTE
%   Calls 
%   Daniel Codiga (2017). UTide Unified Tidal Analysis and Prediction 
%   Functions. MATLAB Central File Exchange.
%   https://www.mathworks.com/matlabcentral/fileexchange/46523-utide-unified-tidal-analysis-and-prediction-functions
% SEE ALSO
%   muiDataSet.m, ctWaterLevelData.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%    
    properties (Hidden)      
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties: 
        Constituents      %tidal constituents generated by utide
    end  
    
    methods %(Access={?muiDataSet,?muiStats,?ct_data_cleanup})
        function obj = ctTidalAnalysis()
            %constructor to initialise object
        end
    end           
    %%   
    methods (Static)               
%% -- model calculation functions --->           
        function obj = runTidalAnalysis(mobj)
            %check to see if there are existing models and add record

            obj = ctTidalAnalysis();         %instantiate object
            dsp = modelDSproperties(obj);
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in InWave
            if ~isValidModel(mobj, metaclass(obj).Name)
                warndlg('No data available. Load water level data to do tidal analysis');
                return;
            end      
            muicat = mobj.Cases;
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------       
            promptxt = 'Select water level data set:';           
            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                                'CaseClass',{'ctWaterLevelData'},...
                                'ListSize',[300,100]);                                        
            if ok<1, obj = []; return; end

            %assign the run parameters to the model instance
            setRunParam(obj,mobj,caserec); 

            %select a variable from the water level dataset
            wldst = getDataset(muicat,caserec,1);            
            varnames = wldst.VariableNames;
            idx = 1;
            if length(varnames)>1
                [idx,ok] = listdlg('Name','WL options', ...
                    'PromptString','Select a variable:', ...
                    'SelectionMode','single','ListSize',[200,100],...
                    'ListString',varnames);
                if ok<1, obj = []; return; end
            end
            
            %select portion of record to use
            range = wldst.RowRange;
            defaults = {char(range{1}),char(range{2})};
            promptxt = {'Start date', 'End date'};
            values = inputdlg(promptxt,'Subsample record',1,defaults);
            if isempty(values) %user cancelled
                obj = [];
                return;
            else
                %offset ensures selected range is extracted  
                fmt = range{1}.Format;
                startime = datetime(values{1},'InputFormat',fmt)-minutes(1);  
                endtime = datetime(values{2},'InputFormat',fmt)+minutes(1);  
                if startime>endtime
                    warndlg('Start time must preceed End time')
                    return; 
                end
            end  
            timeidx = isbetween(wldst.RowNames,startime,endtime);
            time = wldst.RowNames(timeidx);
            data = wldst.(varnames{idx})(timeidx);

            %setup options for utide
            hw = waitbar(0,'Processing. Please wait');
            numtime = convertTo(time,'datenum'); %format required by ut_solv
            waitbar(0.2,hw);
            lat = [];
            cnstit = 'auto';
            %see utide manual for full list of analysis options
            options = {'NodsatNone','White','NoDiagn'};
            
            coef = ut_solv(numtime,data,[],lat,cnstit,options{:});
            predicted = ut_reconstr(numtime,coef);
            diff = data-predicted;
            plotUTideLevels(obj,time,data,predicted,diff)
            delete(hw)

            answer = questdlg('Save results?','Tidal analysis','Yes','No','Yes');
            if strcmp(answer,'No'), obj = []; return; end
            
            %assign the model results to the dstable
            dst = dstable(predicted,'RowNames',time,'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = sprintf('Tides using: %s',wldst.Description);
            %save results
            obj.Constituents = coef;
            setDataSetRecord(obj,muicat,dst,'model');            
            getdialog('Run complete');
        end
%%
        function obj = getTidalData(mobj)
            %check to see if there are existing models and add record
            muicat = mobj.Cases;
            if ~any(strcmp(muicat.Catalogue.CaseClass,'ctTidalAnalysis'))
                msgtxt = 'No constituents available. Run tidal analysis on water level data';
                warndlg(msgtxt)
                return
            end
            
            promptxt =  'Select tidal constituent data set:';
            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                                'CaseClass',{'ctTidalAnalysis'},...
                                'ListSize',[300,100]);
            if ok<1, return; end
            cobj = getCase(muicat,caserec);
            dst = cobj.Data.Dataset;
            coef = cobj.Constituents;

            %get the time period to reconstruct
            startdate = convertTo(dst.RowNames(1),'datenum');%format required by uigetdate
            enddate = convertTo(dst.RowNames(end),'datenum');             
            stnd = {startdate,enddate};
            title = {'From:','To:'};

            datestnd = NaT(2,1);
            ok = 0;
            while ok==0
                for i=1:2
                    %uigetdate is from Matlab Forum (copyright Elmar Tarajan)
                    %input is serial date number or character string
                    ti = datetime(uigetdate(stnd{i},title{i}),...
                                           'ConvertFrom','datenum'); 
                    datestnd(i) = ti;                
                end
                %
                if datestnd(1)<datestnd(2)
                    ok = 1; 
                else
                    hw = warndlg('Start time must preceed End time');
                    waitfor(hw)
                end
            end

            interval = inputdlg('Record time interval (mins)?',...
                                            'Tidal analysis',1,{'60'});
            if isempty(interval), return; end

            time = (datestnd(1):minutes(str2double(interval{1})):datestnd(2))';            
            numtime = convertTo(time,'datenum'); %format required by ut_reconstr
            
            hw = waitbar(0,'Processing. Please wait');
            %reconstruct the tidal record, time input is datenum
            %see utide manual for full list of analysis options
            options = {};
            predicted = ut_reconstr(numtime,coef,options{:});
            delete(hw)
            if isempty(predicted), return; end

            %assign the model results to the dstable
            obj = ctTidalAnalysis();            %instantiate object
            dsp = modelDSproperties(obj);
            time.Format = dst.RowNames.Format; %force time format to same as source
            dst = dstable(predicted,'RowNames',time,'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = sprintf('Tides using: %s',dst.Description);
            %save results
            
            obj.Constituents = coef;
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
%%
        function plotUTideLevels(~,time,data,predicted,diff)
            %plot the input data, the predicted data and the differences
            figure('Name','Tidal Analysis','Tag','PlotFig')
            subplot(3,1,1)
            plot(time,data)
            ylabel('Input water levels')
            
            subplot(3,1,2)
            yyaxis right
            plot(time,predicted)
            ylabel('Predicted water levels')
            
            subplot(3,1,3)
            plot(time,diff)
            ylabel('Residual differences (m)')
        end
    end
%%    
    methods (Access = private)
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                      
                'Name',{'TLOD'},...
                'Description',{'Tidal level'},...
                'Unit',{'mOD'},...
                'Label',{'Tidal level (mOD)'},...
                'QCflag',{'model'}); 
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