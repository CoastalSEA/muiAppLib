function output = wl_scast_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wl_scast_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   ShoreCast data format (vector time and date)
% USAGE
%   output = wl_scast_format(funcall,varargin)
% USAGE
%   obj = dataimport_format_template(obj,funcall)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
%
% Author: Ian Townend
% CoastalSEA (c)June 2017
%--------------------------------------------------------------------------
%
    switch funcall
        case 'getData'
          output = getData(varargin{:});
        case 'dataQC'
            output = dataQC(varargin{1});
        case 'getPlot'
            output = 0;
    end
end
%%
%--------------------------------------------------------------------------
% getData
%--------------------------------------------------------------------------
function dst = getData(~,filename)
    %read and load a data set from a file
    [data,~] = readInputData(filename);             
    if isempty(data), dst = []; return; end
    
    %set metadata
    dsp = setDSproperties;

    % assign date and time
    myDatetime = data{1};
    myDatetime.Format = dsp.Row.Format;
    
    % information on data location
%     Latitude = [];
%     Longitude = [];
    
    varData = table(data{1,3:end});
    varData = standardizeMissing(varData,[99,99.9,99.99,999,9999]);

    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
%     dst.Dimensions.Position = [Latitude,Longitude];    
end
%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    dataSpec = '%u %u %u %u %u %u %f'; 
    nhead = 1;     %number of header lines
    [data,header] = readinputfile(filename,nhead,dataSpec); %see muifunctions
    myDatetime = datetime(data{1},data{2},data{3},data{4},data{5},data{6});
    data =[{myDatetime},{myDatetime},data(7)];
end
%%
%--------------------------------------------------------------------------
% dataDSproperties
%--------------------------------------------------------------------------
function dsp = setDSproperties()
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'WLOD'},...                   
        'Description',{'Water Level'},...
        'Unit',{'mAD'},...
        'Label',{'Water levels (mAD)'},...
        'QCflag',{'raw'}); 
    dsp.Row = struct(...
        'Name',{'Time'},...
        'Description',{'Time'},...
        'Unit',{'h'},...
        'Label',{'Time'},...
        'Format',{'dd-MM-yyyy HH:mm:ss'});        
    dsp.Dimensions = struct(...    
        'Name',{'Position'},...
        'Description',{'Latitude and Longitude'},...
        'Unit',{'deg'},...
        'Label',{'Latitude and Longitude'},...
        'Format',{''});          
end
%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
%%
function output = dataQC(obj)
    %quality control a dataset
    datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    dst = obj.Data.(datasetname);      %selected dstable
    dst = wl_dataQC(dst);        %generic data quality control for water levels
    output = {dst,datasetname};
end