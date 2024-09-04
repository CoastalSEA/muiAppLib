function output = shoreline_data_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   shoreline_data_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   Shoreline position data
% USAGE
%   output = shoreline_data_format(funcall,varargin)
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
    selection = obj.DataFormats{obj.idFormat,1};
    %set metadata
    dsp = setDSproperties;
    
    [data,~] = readInputData(filename,selection,dsp);             
    if isempty(data), dst = []; return; end

    % concatenate date and time
    mdat = data{1};
    mtim = data{2};
    mdat.Format = dsp.Row.Format;
    mtim.Format = dsp.Row.Format;
    myDatetime = mdat + timeofday(mtim);
    
    % information on data location
%     Latitude = [];
%     Longitude = [];
    
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
function [data,header] = readInputData(filename,selection,dsp)
    %read wshoreline data (read format is file specific).
        
    switch selection
        case 'Date-Time-Position'
            %File format defined in header line
            %Date Time Shoreline position
            %eg:  '%{dd/MM/yyyy}D %{HH:mm:ss}D %q';
            dataSpec = [];           
        case 'Vector-Date-Position'
            %header line for vector data format (eg ShoreCast data)
            dataSpec = '%u %u %u %u %u %u %f';   
    end
    nhead = 1;
    [data,header] = readinputfile(filename,nhead,dataSpec); %see muifunctions
    %concatenate the date vector into a datetime
    if strcmp(selection, 'Vector-Date-Position')
        myDatetime = datetime(data{1},data{2},data{3},...
                                         data{4},data{5},data{6});
        myDatetime.Format = dsp.Row.Format;
        data =[{myDatetime},{myDatetime},data(7)];
    end    
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
        'Name',{'Chainage'},...                   % <<Edit metadata to suit model
        'Description',{'Shoreline position'},...
        'Unit',{'m'},...
        'Label',{'Shoreline (m)'},...
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
function output = dataQC(obj)
    %quality control a dataset
    % datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    % dst = obj.Data.(datasetname);      %selected dstable
    warndlg('No quality control defined for this format');
    output = [];    %if no QC implemented in dataQC
end