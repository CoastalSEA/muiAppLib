function output = wl_gesla_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wl_gesla_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   GESLA water level data
% USAGE
%   obj = wl_gesla_format(funcall,varargin)
% INPUTS
%   obj - instance of WaterLevelData class object
%   funcall - function being called  
%             getData loads the dataset into the class object
%                - calls readInputData to read data from file
%             dataQC handles quality control of data
% OUTPUT
%   output - updated class object
% NOTES
%   Global Extreme Sea Level Analysis data format
%   https://gesla.org/
%   for details of format see: https://gesla787883612.wordpress.com/format/
%   updated to GESLA v5.0 on 15 Feb 2022
%
% Author: Ian Townend
% CoastalSEA (c)Jan 2021
%--------------------------------------------------------------------------
%
    switch funcall
        case 'getData'
            output = getData(varargin{:});
        case 'dataQC'
            output = dataQC(varargin{:});
        case 'getPlot'
            output = 0;
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
        'Name',{'WL2D','WLflag'},... 
        'Description',{'Water Level to Datum','Flag'},...
        'Unit',{'mAD','-'},...
        'Label',{'Water Level (mAD)','Flag'},...
        'QCflag',repmat({'raw'},1,2)); 
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
% getData
%--------------------------------------------------------------------------
function dst = getData(~,filename)
    %read data from file (function is at end of file)
    [data,header] = readInputData(filename);  
    if isempty(data), dst = []; return; end
    %set metadata
    dsp = setDSproperties;
    
    % concatenate date and time
    mdat = data{1};
    mtim = data{2};
    mdat.Format = dsp.Row.Format;
    mtim.Format = dsp.Row.Format;
    rdata = mdat + timeofday(mtim);
    
    % information on data location
    temp = split(header{11});           %depends on GESLA format
    Latitude = str2double(temp{3});
    temp = split(header{12});           %depends on GESLA format
    Longitude = str2double(temp{3});
    
    %capture datum info as metadata
    mdata = sprintf('Datum used: %s',header{18}(21:end)); %depends on GESLA format
    
    %check that datetime values are unique
    [UniqueTime,iu] = unique(rdata);
    if length(UniqueTime)~=length(rdata)
        rdata = UniqueTime;                
    end    
    
    varData = table(data{1,3:end-1});
    varData = standardizeMissing(varData,[99,99.9,99.99,999,9999]);
    varData = varData(iu,:);
    
    %load the results into a dstable  
    dst = dstable(varData,'RowNames',rdata,'DSproperties',dsp); 
    dst.Dimensions.Position = [Latitude,Longitude];
    dst.MetaData = mdata;
end
%%
function [data,header] = readInputData(filename)
    %read wave data (read format is file specific).
    dataSpec = '%{yyyy/MM/dd}D %{HH:mm:ss}D %f %u %u';
    nhead = 41;  %see: https://gesla787883612.wordpress.com/format/
    [data,header] = readinputfile(filename,nhead,dataSpec); %see muifunctions
end
%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
function output = dataQC(obj)
    %quality control a user data timeseries
    %the obj passed to the function is the ts collection  
%   Quality-control (QC) flags 
%   0 - no quality control
%   1 - correct value 
%   2 - interpolated value
%   3 - doubtful value
%   4 - isolated spike or wrong value
%   5 - missing value
    datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    dst = obj.Data.(datasetname);      %selected dstable
    
    flag = dst.WLflag;
    idx = flag>1;
    dst.WL2D(idx) = NaN;
    
    WL2D = dst.WL2D;
    figure('Name','QC Plot','Tag','PlotFig');
    plot(dst,'WL2D');                  
    maxWL = ceil(max(WL2D));
    minWL = floor(min(WL2D));
    %get input parameters from user
    prompt = {'Maximum water level (mOD):','Minimum water level (mOD):'};
    title = 'Define limiting water levels';
    numlines = 1;
    defaultvalues = {num2str(maxWL),num2str(minWL)};
    useInp=inputdlg(prompt,title,numlines,defaultvalues);
    if isempty(useInp), return; end %user cancelled
    maxWL = str2double(useInp{1});
    minWL = str2double(useInp{2});

    idx = idx | (WL2D>maxWL);
    idx = idx | (WL2D<minWL);

    WL2D(idx) = NaN;

    hw = waitbar(0, 'Loading data. Please wait');            

    dst.WL2D = WL2D;            
    dst.VariableQCflags(1:2) = repmat({'qc'},1,2); 

    output = {dst,datasetname};    
    waitbar(1); 
    close(hw);
end