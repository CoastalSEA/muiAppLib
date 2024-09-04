classdef ctShorelineData < muiDataSet
%
%-------class help---------------------------------------------------------
% NAME
%   ctShorelineData.m
% PURPOSE
%   Class to hold shoreline position data
% USAGE
%   obj = ctShorelineData()
% SEE ALSO
%   inherits muiDataSet and uses dstable and dscatalogue
%   format files used to load data of varying formats (variables and file format)
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%     
    properties  
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:  
    end
        
    methods  
        function obj = ctShorelineData()
            %constructor to initialise object

            %initialise list of available input file formats. Format is:
            %{'label 1','function name 1';'label 2','function name 2'; etc}
            obj.DataFormats = {'Date-Time-Position','shoreline_data_format';...
                               'Vector-Date-Position','shoreline_data_format'};   
            %define file specification, format is: {multiselect,file extension types}
            obj.FileSpec = {'on','*.txt;*.csv'};               
        end
%%
        function tabPlot(obj,src)
            %generate plot for display on Q-Plot tab
            tabDefaultPlot(obj,src);
        end         
    end  
%--------------------------------------------------------------------------
%   functions to read data from file and load as a timeseries
%   default function for loadData, addData and dataQC are in muiDataSet
%   format defintions are in the format files
%--------------------------------------------------------------------------
end