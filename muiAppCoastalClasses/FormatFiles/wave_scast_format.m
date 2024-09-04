function output = wave_scast_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wave_scast_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   ShoreCast data format (vector time and date)
% USAGE
%   output = wave_scast_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - updated class object
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

    % default information on data location
%     Latitude =0;
%     Longitude = 0;
    myDatetime = datetime(data{1},data{2},data{3},data{4},data{5},data{6});
    % concatenate date and time
    myDatetime.Format = dsp.Row.Format;

    [UniqueTime,iu] = unique(myDatetime);
    if length(UniqueTime)~=length(myDatetime)
        myDatetime = UniqueTime;                
    end
    varData = [{data{7}(iu)},{data{8}(iu)},{data{11}(iu)}];    
    varData = table(varData{1,:});
    varData = standardizeMissing(varData,[99,99.9,99.99,999,9999]);

    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
%     dst.Dimensions.Position = [Latitude,Longitude];    
end
%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    dataSpec = '%u %u %u %u %u %u %f %f %f %f %f'; 
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
        'Name',{'Hs','Tp','Dir'},...                   % <<Edit metadata to suit model
        'Description',{'Significant wave height','Peak period','Wave direction'},...
        'Unit',{'m','s','deg'},...
        'Label',{'Wave height (m)','Wave period (s)','Wave direction (deg)'},...
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
    %quality control a user wave data timeseries    
    datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    dst = obj.Data.(datasetname);      %selected dstable
    dst = wave_dataQC(dst);            %generic data quality control for waves
    output = {dst,datasetname};
end