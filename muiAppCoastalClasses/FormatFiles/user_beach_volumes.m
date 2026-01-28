function output = user_beach_volumes(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   user_beach_volumes.m
% PURPOSE
%   Functions to define metadata, read and load data from file to
%   import beach volumes, adding the results to dstable
%   and a record in a dscatlogue (as a property of muiCatalogue)
% USAGE
%   obj = user_beach_volumes(obj,funcall)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   Data format is date rows and column profile volumes
%   Header has 2 rows. first row defines data format, second row the
%   profile numbers.
%   Data is added to class muiUserData
%
% Author: Ian Townend
% CoastalSEA (c)Feb 2026
%--------------------------------------------------------------------------
%
    switch funcall
        %standard calls from muiDataSet - do not change if data class 
        %inherits from muiDataSet. The function getPlot is called from the
        %Abstract method tabPlot. The class definition can use tabDefaultPlot
        %define plot function in the class file, or call getPlot
        case 'getData'
          output = getData(varargin{:});
        case 'dataQC'
            output = dataQC(varargin{1});  
        case 'getPlot'
            output = 0; %uses the default tab plot in muiDataSet, else
            %output = getPlot(varargin{:});
    end
end
%%
%--------------------------------------------------------------------------
% getData
%--------------------------------------------------------------------------
function dst = getData(obj,filename) %#ok<INUSD>
    %read and load a data set from a file
    [data,header] = readInputData(filename);             
    if isempty(data), dst = []; return; end
    
    %set metadata
    dsp = setDSproperties;

    %code to parse input data and assign to varData  
    myDatetime = data{1};                            %Read data column
    varData = [data{2:end}];                              %
                                                    
    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
    pid = split(string(header{2}));
    dst.Dimensions.Pid = str2double(pid(2:end-1));  %omit text and return
end
%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    dataSpec = []; 
    nhead = 2;     %number of header lines
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
        'Name',{'Volume'},...  
        'Description',{'Beach volume'},...
        'Unit',{'m^3'},...
        'Label',{'Beach volume (m^3)'},...
        'QCflag',{'raw'}); 
    dsp.Row = struct(...
        'Name',{'Time'},...
        'Description',{'Time'},...
        'Unit',{'h'},...
        'Label',{'Time'},...
        'Format',{'dd/MM/yyyy'});        
    dsp.Dimensions = struct(...    
        'Name',{'Pid'},...
        'Description',{'Profile Number'},...
        'Unit',{'-'},...
        'Label',{'Profile No.'},...
        'Format',{''});         
end
%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
function output = dataQC(obj)                        % <<Add any quality control to be applied (optional)
    %quality control a dataset
    % datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    % dst = obj.Data.(datasetname);      %selected dstable
    warndlg('No quality control defined for this format');
    output = [];    %if no QC implemented in dataQC
end
%%
%--------------------------------------------------------------------------
% getPlot
%--------------------------------------------------------------------------
function ok = getPlot(obj,src)                       % <<Add code for bespoke Q-Plot is required (optional)
    %generate a plot on the src graphical object handle    
    ok = 0;  %ok=0 if no plot implemented in getPlot
    %return some other value if a plot is implemented here
end



