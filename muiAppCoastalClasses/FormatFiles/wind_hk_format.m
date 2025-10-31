function output = wind_hk_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%  wind_hk_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   Hong Kong wind data format
% USAGE
%   output = wind_hk_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
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
function dst = getData(obj,filename) 
    %read and load a data set from a file
    [data,~] = readInputData(filename);             
    if isempty(data), dst = []; return; end
    
    %set metadata
    dsp = setDSproperties;

    % unpack data for time variables
    mdat = data{1};
    mtim = data{2};
    
    idx = mtim==24;
    mdat(idx) = mdat(idx)+1;
    mtim(idx) = 0;
    
    mdat = datetime(mdat,'ConvertFrom','yyyymmdd');
    mtim = hours(mtim);
    %{HH}D 
    % concatenate date and time
    myDatetime = mdat + mtim;
    myDatetime.Format = dsp.Row.Format;
    %remove text string flags
    data(:,3:6) = cellfun(@str2double,data(:,3:6),'UniformOutput',false);
    %reorder to be speed direction speed direction
    temp = data(:,3);
    data(:,3) = data(:,4);
    data(:,4) = temp;
    temp = data(:,5);
    data(:,5) = data(:,6);
    data(:,6) = temp;
    %check for missing data
    varData = table(data{1,3:end});
    varData = standardizeMissing(varData,[99,99.9,99.99,999,9999]);
    % information on data location (could be in header)
    % Latitude = [];
    % Longitude = [];

    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
    dst.MetaData.zW = setHeight(obj);
    % dst.Dimensions.Position = [Latitude,Longitude];    
end

%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    dataSpec = '%d %d %s %s %s %s'; 
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
        'Name',{'Speed10min','Dir10min','Speed1hr','Dir1hr'},...    
        'Description',{'Mean wind speed 10min','Mean wind direction 10min',...
                   'Mean wind speed 1hr','Mean wind direction 1hr'},...
        'Unit',{'m/s','deg','m/s','deg'},...
        'Label',{'Wind speed (m/s)','Wind direction (deg)',...
                   'Wind speed (m/s)','Wind direction (deg)'},...
        'QCflag',repmat({'raw'},1,4));
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
    %quality control a dataset
    % datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    % dst = obj.Data.(datasetname);      %selected dstable
    warndlg('No quality control defined for this format');
    output = [];    %if no QC implemented in dataQC
end