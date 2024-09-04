classdef ctModelData < muiDataSet                    
%
%-------class help------------------------------------------------------===
% NAME
%   ctModelData.m
% PURPOSE
%   Class to hold wave data
% USAGE
%   obj = ctModelData() 
% SEE ALSO
%   inherits muiDataSet and uses dstable and dscatalogue
%   format files used to load data of varying formats (variables and file format)
%
% Author: Ian Townend
% CoastalSEA (c) Apr 2024
%--------------------------------------------------------------------------
%    
    properties  
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        % importing data requires muiDataSet propertiesm DataFormats and
        % FileSpec to be defined in class constructor.
        %Additional properties:  
    end
    
    methods 
        function obj = ctModelData()              
            %class constructor
            %initialise list of available input file formats. Format is:
            %{'label 1','function name 1';'label 2','function name 2'; etc}
            obj.DataFormats = {'BlueKenu timeseries data','bluekenue_ts_format'};
            %define file specification, format is: {multiselect,file extension types}
            obj.FileSpec = {'on','*.txt;*.csv;*.ts*'};
        end
%%
        function tabPlot(obj,src)
            %generate plot for display on Q-Plot tab
            tabDefaultPlot(obj,src);
        end 
    end
end