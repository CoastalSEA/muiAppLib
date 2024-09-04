function output = wavenet_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wavenet_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   cefas WaveNet data
% USAGE
%   output = wavenet_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   ouput - updated class object
% NOTES
%   cefas WaveNet data format
    % https://www.cefas.co.uk/data-and-publications/wavenet/
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
    [data,~] = readInputData(filename);             
    if isempty(data), dst = []; return; end
    
    %set metadata
    dsp = setDSproperties;
    
    % concatenate date and time
    mdat = data{1};
    mtim = data{2};
    mdat.Format = dsp.Row.Format;
    mtim.Format = dsp.Row.Format;
    myDatetime = mdat + timeofday(mtim);

    % information on data location
%     Latitude  = [];
%     Longitude = [];
    
    %re-order to Hs, Tp, Tz, Dir, spread, temp
    varData = table(data{1,6},data{1,5},data{1,3},data{1,4},data{1,7},data{1,8});
    varData = standardizeMissing(varData,[99,99.9,99.99,999,9999]);
    
    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
%     dst.Dimensions.Position = [Latitude,Longitude];    
end
%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    dataSpec ='%{yyyy-MM-dd}D %{HH:mm:ss}D %f %d %f %f %f %d';
    nhead = 1;     %number of header lines
    [data,header] = readinputfile(filename,nhead,dataSpec); %see muifunctions
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
        'Name',{'Hs','Tp','Tz','Dir','Spr','SST'},...              
        'Description',{'Significant wave height','Peak period',...
                   'Mean zero-crossing period','Wave direction',...
                   'Spread','Sea surface temperature'},...
        'Unit',{'m','s','s','deg','deg','deg C'},...
        'Label',{'Wave height (m)','Wave period (s)',...
                   'Wave period (s)','Wave direction (deg)',...
                   'Wave spread (deg)','Sea surface temperature (deg)'},...
        'QCflag',repmat({'raw'},1,6)); 
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
function output = dataQC(obj)
    %quality control a user data timeseries    
    datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    dst = obj.Data.(datasetname);      %selected dstable
    dst = wave_dataQC(dst);            %generic data quality control for waves
    output = {dst,datasetname};
end