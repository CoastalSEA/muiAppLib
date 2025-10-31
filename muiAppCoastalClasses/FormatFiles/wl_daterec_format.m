function output = wl_daterec_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wl_daterec_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   Date-Record data format
% USAGE
%   output = wl_daterec_format(funcall,varargin)
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
    %%% NB: this no longer uses the file header %%%
    data = readTSinputFile(obj,filename);   %data returned as a table
    if isempty(data), dst = []; return; end

    if ischar(data{1,1}{1})
        varData = data(:,4:end);    %format with id) date time data
        % concatenate date and time
        mdat = datetime(data{:,2},'InputFormat','yyyy/MM/dd');
        mtim = data{:,3};
        myDatetime = mdat + mtim;        
    elseif isdatetime(data{1,1})
        varData = data(:,2:end);    %format with datetime data
        %extract datetime
        myDatetime = data{:,1};
    else
        dst = [];
        warndlg('Unrecognised file format')
        return;
    end
    varData = addDatum(obj,varData,false); 
    %extract required subset of data
    %set metadata
    nvar = width(varData);
    dsp = setDSproperties(nvar);         
    myDatetime.Format = dsp.Row.Format;

    % information on data location
    % Latitude = [];
    % Longitude = [];        

    %load the table varData into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
    % dst.Dimensions.Position = [Latitude,Longitude];    
end
%%
%--------------------------------------------------------------------------
% dataDSproperties
%--------------------------------------------------------------------------
function dsp = setDSproperties(nvar)
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'WLOD','WLCD','WLResidual','WLFlag'},...            
        'Description',{'Water Level to OD','Water Level to CD','Residual','Flag'},...
        'Unit',{'mAD','mCD','m','-'},...
        'Label',{'Water Levels (mAD)','Water Levels (mCD)','Residual (m)','Flag'},...
        'QCflag',repmat({'raw'},1,4)); 
    dsp.Row = struct(...
        'Name',{'Time'},...
        'Description',{'Time'},...
        'Unit',{'h'},...
        'Label',{'Time'},...
        'Format',{'dd-MM-yyyy HH:mm:ss'});        
    dsp.Dimensions = struct(...    
        'Name',{''},...
        'Description',{''},...
        'Unit',{''},...
        'Label',{''},...
        'Format',{''});    
    %remove excess variables
    if nvar==2
        dsp.Variables(3:4) = [];
    elseif nvar==3
        dsp.Variables(4) = [];
    end   
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