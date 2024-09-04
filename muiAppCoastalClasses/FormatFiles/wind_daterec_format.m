function output = wind_daterec_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%  wind_daterec_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   Date-Record data format
% USAGE
%   output = wind_daterec_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - updated class object
% NOTES
%   Date-Record data format
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

    % default information on data location
%     Latitude = 0;
%     Longitude = 0; 
    % concatenate date and time
    mdat = data{1};
    mtim = data{2};
    mdat.Format = dsp.Row.Format;
    mtim.Format = dsp.Row.Format;
    myDatetime = mdat + timeofday(mtim);
    
    %check for duplicate records
    [UniqueTime,iu] = unique(myDatetime);
    if length(UniqueTime)~=length(myDatetime)
        myDatetime = UniqueTime;                
    end
    
    %check for missing data
    varData = table(data{1,3:end});
    varData = standardizeMissing(varData,[99,99.9,99.99,999,9999]);
    varData = varData(iu,:);

    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp);
%     dst.Dimensions.Position = [Latitude,Longitude];    
end
%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    dataSpec = []; 
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
        'Name',{'AvSpeed','MaxSpeed','Dir'},...              
        'Description',{'Mean wind speed','Maximum wind speed','Mean wind direction'},...
        'Unit',{'m/s','m/s','deg'},...
        'Label',{'Wind speed (m/s)','Wind speed (m/s)','Wind direction (deg)'},...
        'QCflag',repmat({'raw'},1,3)); 
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
function output = dataQC(obj)
    %quality control a user data timeseries
    % datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    % dst = obj.Data.(datasetname);      %selected dstable
    warndlg('No quality control defined for this format');
    output = [];    %if no QC implemented in dataQC
end