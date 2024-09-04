classdef ctWindWaveModel < muiDataSet                        
%
%-------class help------------------------------------------------------
% NAME
%   Model_template.m
% PURPOSE
%   Class to handle WindWaveModel methods and output for use in
%   CoastalTools and other ModellUI apps
%
% SEE ALSO
%   muiDataSet
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:     
    end
    
    methods (Access={?muiDataSet,?muiStats}) 
        function obj = ctWindWaveModel()                      
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = ctWindWaveModel;                            
            dsp = modelDSproperties(obj);
            
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
%--------------------------------------------------------------------------
% Model code  <<INSERT MODEL CODE or CALL MODEL>>
%--------------------------------------------------------------------------
            %get the timeseries input data and site parameters
            [dst,inputxt] = getInputData(obj,mobj);
            if isempty(dst), return; end   %user cancelled data selection
            %get the input site parameters as a class instance
            site_params = mobj.Inputs.ctHindcastParameters;  
            inp = inputParameters(site_params); %convert class to struct using short names 
                                                %method is in ctHindcastParameters

            [F_dir,F_len] = readfetchfile();
            if isempty(F_dir), return; end %no fetch data
            Fch = eff_fetch(F_len, F_dir, dst.Dir, inp.n, inp.Fopt);

            [Hs,Tp,Tz] = tma_spectrum(dst.U,inp.zw,Fch,inp.df,inp.ds);
            if isempty(Hs), return; end %no wave data
            
            results = {Hs,Tp,Tz,dst.Dir};
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(results{:},'RowNames',dst.RowNames,'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = inputxt;
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
        function [tsdst,inputxt] = getInputData(obj,mobj)
            %prompt user to select wind data and return in
            %input dstable of data and metadata for inputs used
            tsdst = []; inputxt = [];
            muicat = mobj.Cases;
            promptxt = 'Select input wind data set:';
            [wd_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                            'CaseClass',{'ctWindData'},'ListSize',[300,100]);                                                                        
            if ok<1, return; end
            windst = getDataset(muicat,wd_crec,1);
            inputxt = sprintf('%s used for offshore waves',windst.Description);
            
            varnames = windst.VariableNames;            
            idu = 1;
            if length(varnames)>1
                promptxt = 'Select wind speed variable:';
                [idu,ok] = listdlg('Name','Wind variables', ...
                            'PromptString',promptxt,'ListSize',[200,100],...
                            'SelectionMode','single','ListString',varnames);
                if ok<1, return; end
            end
            idd = 1;
            if length(varnames)>1
                promptxt = 'Select wind direction variable:';
                [idd,ok] = listdlg('Name','Wind variables', ...
                            'PromptString',promptxt,'ListSize',[200,100],...
                            'SelectionMode','single','ListString',varnames);
                if ok<1, return; end
            end
            
            tsdst = getDSTable(windst,[],[idu,idd]);
            tsdst.VariableNames = {'U','Dir'};
            %assign the run parameters to the model instance
            setRunParam(obj,mobj,wd_crec); %input caserecs passed as varargin  
        end
%%
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                       % <<Edit metadata to suit model
                'Name',{'Hs','Tp','Tz','Dir'},...
                'Description',{'Significant wave height','Peak wave period',...
                           'Mean zero crossing period','Wave direction'},...
                'Unit',{'m','s','s','deg'},...
                'Label',{'Wave height (m)','Wave period (s)',...
                           'Wave period (s)','Wave direction (deg)'},...
                'QCflag',repmat({'model'},1,4));  
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