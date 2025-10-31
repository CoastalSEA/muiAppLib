classdef ctWindData < muiDataSet
%
%-------class help---------------------------------------------------------
% NAME
%   ctWindData.m
% PURPOSE
%   Class to hold wind data
% USAGE
%   obj = ctWindData()
% SEE ALSO
%   inherits muiDataSet and uses dstable and dscatalogue
%   format files used to load data of varying formats (variables and file format)
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%  
    properties  
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:  
        zW   %height of wind measurement
    end
    
    methods  
        function obj = ctWindData()
            %class constructor
            %initialise list of available input file formats. Format is:
            %{'label 1','function name 1';'label 2','function name 2'; etc}
            obj.DataFormats = {'Date-Record format','wind_daterec_format';...
                               'CEDA MIDAS format','wind_midas_format';...
                               'Hong Kong data format','wind_hk_format'};
            %define file specification, format is: {multiselect,file extension types}
            obj.FileSpec = {'on','*.txt;*.csv'};               
        end

%%
        function tabPlot(obj,src)
            %generate plot for display on Q-Plot tab
            tabDefaultPlot(obj,src);
        end  

%%
        function zW = setHeight(~)
            %prompt user to define the height of measurement
            promptxt = {'Height of measurement above msl (m)'};
            inpt = inputdlg(promptxt,'Wind height',1,{'10'});
            if isempty(inpt)
                zW = NaN; 
            else
                zW = str2double(inpt{1});
            end 
        end
    end
end