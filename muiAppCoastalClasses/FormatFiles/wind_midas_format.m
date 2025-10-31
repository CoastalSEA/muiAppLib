function output = wind_midas_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%  wind_midas_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   CEDA MIDAS data format
% USAGE
%   output = wind_midas_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   CEDA MIDAS data format
%   https://catalogue.ceda.ac.uk/
%   data.ceda.ac.uk/badc/ukmo-midas-open/data/uk-mean-wind-obs/dataset-version-201908/
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
    %use readTable with auto detect import options  
    data = readTSinputFile(obj,filename);           
    if isempty(data), dst = []; return; end
    
    %set metadata
    dsp = setDSproperties;

    % default information on data location
    % Latitude = 0;
    % Longitude = 0; 

    %extract required subset of data        
    myDatetime = data.ob_end_time;
    myDatetime.Format = dsp.Row.Format;
    
    %convert knots to m/s
    factor = 0.51444;
    vars = {'mean_wind_speed','mean_wind_dir','max_gust_speed','max_gust_dir'};
    varData = data(:,vars);
    varData{:,1} = varData{:,1}*factor;
    varData{:,3} = varData{:,3}*factor;

    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp);
    dst.MetaData.zW = setHeight(obj);
    % dst.Dimensions.Position = [Latitude,Longitude];    
end

%% ------------------------------------------------------------------------
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
        'Name',{'meanwindspeed','meanwinddir','maxgustspeed','maxgustdir'},...               
        'Description',{'Mean wind speed','Mean wind direction',...
                   'Speed of maximum wind gust','Direction of maximum wind gust'},...
        'Unit',{'m/s','degTN','m/s','degTN'},...
        'Label',{'Wind speed (m/s)','Wind direction (degTN)',...
                   'Wind speed (m/s)','Wind direction (degTN)'},...
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
