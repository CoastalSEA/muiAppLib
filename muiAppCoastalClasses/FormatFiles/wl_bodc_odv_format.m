function output = wl_bodc_odv_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wl_bodc_odv_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   BODC water level data ODV format
% USAGE
%   obj = wl_bodc_odv_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   British Oceanographic Data Centre: www.bodc.ac.uk
%
% Author: Ian Townend
% CoastalSEA (c)Feb 2021
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
            output = 0;
    end
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
        'Name',{'WLOD','WLCD','WLResidual'},... 
        'Description',{'Water Level to OD','Water Level to CD','Residual'},...
        'Unit',{'mOD','mCD','m'},...
        'Label',{'Water Level (mOD)','Water Level (mCD)','Residual (m)'},...
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
    if nvar==2  %only water level data
        dsp.Variables(3) = [];
    end   
end
%%
%--------------------------------------------------------------------------
% getData
%--------------------------------------------------------------------------
function dst = getData(obj,filename)
    %read and load a data set from a file
    [data,header] = readInputData(filename); %data returned as cell array              
    if isempty(data), dst = []; return; end
       
    %extract variables
    varData = table(data{3:end});
    varData = addDatum(obj,varData,true); 
    varData = standardizeMissing(varData,[99,99.9,99.99,999,9999,-99]);

    nvar = width(varData);
    %set metadata
    dsp = setDSproperties(nvar);   

    % concatenate date and time
    rdata = data{1};
    rdata.Format = dsp.Row.Format;

    % information on data location
    Latitude = header{6};
    Longitude = header{5};

    %load the results into a dstable  
    dst = dstable(varData,'RowNames',rdata,'DSproperties',dsp); 
    dst.Dimensions.Position = [Latitude,Longitude];
end
%%
function [data,header] = readInputData(filename)
    %read water level data (read format is file specific).
    fid = fopen(filename, 'r');
    if fid<0
        errordlg('Could not open file for reading','File write error','modal')
        header = ''; data = [];
        return;
    end
    
    %read BODC data using Ocean Data View (ODV) format
    tline = fgetl(fid); hflag = true;
    while ischar(tline) && hflag
        disp(tline)
        tline = fgetl(fid);
        if length(tline)>=11
            if strcmp(tline(1:11),'unspecified')
                hflag = false;
            end
        end
    end
    %find number of columns in file
    delimiter = sprintf('\t'); %tab
    ncols = numel(strfind(tline,delimiter)) + 1;
    %spec of first line of data which includes location information
    dataspec = ['%s %s %s %s %f %f %f %u %u', repmat([' %f',' %u'],1,floor((ncols-9)/2))];
    row1 = textscan(tline,dataspec);
    %Cruise	Station	Type	yyyy-mm-ddThh:mm:ss.sss	Longitude [degrees_east]	Latitude [degrees_north]	LOCAL_CDI_ID	EDMO_code	Bot. Depth [m]
    header = row1(1:9);
    %Chronological Julian Date [days]	QV:SEADATANET	SeaLvl_bubbler [m]	QV:SEADATANET	SeaLvl_bubbler2 [m]	QV:SEADATANET	HTSeaLvl [m]	QV:SEADATANET
    row1 = row1(10:end); %first row of data
    ncol = length(row1);
    dataSpec = repmat(['%f ','%u '],1,floor(ncol/2));
    rown = textscan(fid,dataSpec); %remaining rows
    for i=1:ncol
        if isempty(row1{1,i})
            row1{1,i} = 0;
        end
        rows{1,i} = vertcat(row1{1,i}(1),rown{1,i}(:));
    end
    %convert Julian data to date-time
    data{1,1} = datetime(rows{1,1},...
        'ConvertFrom','juliandate','Format','dd-MMM-yyyy');
    data{1,2} = datetime(rows{1,1},...
        'ConvertFrom','juliandate','Format','HH:mm:ss'); %not used
    %assign variables                            
    data{1,3} = rows{1,3};

%     data = addDatum(obj,data); 
    if isempty(data)
        warndlg('No data. Check file format selected')
    end
    fclose(fid);
end
%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
function output = dataQC(obj)
    %quality control a dataset
    datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    dst = obj.Data.(datasetname);      %selected dstable
    dst = wl_dataQC(dst);        %generic data quality control for water levels
    output = {dst,datasetname};
end



