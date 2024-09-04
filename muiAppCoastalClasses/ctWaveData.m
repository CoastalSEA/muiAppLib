classdef ctWaveData < muiDataSet                    
%
%-------class help------------------------------------------------------===
% NAME
%   ctWaveData.m
% PURPOSE
%   Class to hold wave data
% USAGE
%   obj = ctWaveData() 
% SEE ALSO
%   inherits muiDataSet and uses dstable and dscatalogue
%   format files used to load data of varying formats (variables and file format)
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%    
    properties  
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        % importing data requires muiDataSet propertiesm DataFormats and
        % FileSpec to be defined in class constructor.
        %Additional properties:  
    end
    
    methods 
        function obj = ctWaveData()              
            %class constructor
            %initialise list of available input file formats. Format is:
            %{'label 1','function name 1';'label 2','function name 2'; etc}
            obj.DataFormats = {'Date-Record format','wave_daterec_format';...
                               'Channel Coastal Observatory format','wave_cco_format';...
                               'Cefas Wavenet format','wavenet_format';...                               
                               'ShoreCast data format','wave_scast_format';...
                               'CCO directional wave spectra','wave_cco_spectra'};
            %define file specification, format is: {multiselect,file extension types}
            obj.FileSpec = {'on','*.txt;*.csv;*.spt'};
        end
%%
        function tabPlot(obj,src)
            %generate plot for display on Q-Plot tab
            tabDefaultPlot(obj,src);
        end 
    end
end