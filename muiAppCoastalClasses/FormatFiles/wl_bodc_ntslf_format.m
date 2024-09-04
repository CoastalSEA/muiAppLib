function output = wl_bodc_ntslf_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wl_bodc_ntslf_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   BODC water level data NTSLF format
% USAGE
%   obj = wl_bodc_ntslf_format(funcall,varargin)
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
            output = dataQC(varargin{:});
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
    mdat = data{1};
    mtim = data{2};
    mdat.Format = dsp.Row.Format;
    mtim.Format = dsp.Row.Format;
    rdata = mdat + timeofday(mtim);
    
    % information on data location
    Latitude = header{1};
    Longitude = header{2};        

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
    
    %read BODC data using National Tidal & Sea Level Facility (NTSLF) format
    tline = fgetl(fid); hflag = true;
    while ischar(tline) && hflag %read header
        disp(tline)
        tline = fgetl(fid);
        if length(tline)>=11
            if strcmp(tline(1:9),'Latitude:')
                header = textscan(tline,'%s %f');
                header = header(1,2);
            elseif strcmp(tline(1:10),'Longitude:')
                temp = textscan(tline,'%s %f');
                header(1,2) = temp(1,2);
            elseif strcmp(tline(1:6),'Number') || strcmp(tline(1:7),' Number')
                hflag = false;
            end
        end   
    end
    %File format
    %Cycle Date Time Tidelevel (flag) Residual (flag) 
    %some data may be flagged so have to read line
    %by line
    tline = fgetl(fid); %remaining rows
    nvar = length(split(strtrim(tline)));
    varstr = repmat('%s ',1,nvar-3);
    intstr = '%s %{yyyy/MM/dd}D %{HH:mm:ss}D';
    dataSpec = sprintf('%s %s',intstr,varstr);
    count = 1;
    rowflag = 1;
    while ischar(tline) && ~isempty(tline)
        temp = textscan(tline,dataSpec);
        for j=1:nvar-3
            varj = char(temp{j+3});
            if isletter(varj(end))
                flag = 9;
                temp{j+3} = str2double(varj(1:end-1));
            else
                flag = 1;
                temp{j+3} = str2double(varj);
            end
        end

        rowflag = max(rowflag,flag);
        temp(1,end+1) = {rowflag};
        rows(count,:) = temp;            
        count = count+1;
        tline = fgetl(fid);
    end

    for i=1:nvar
        rown{1,i} = vertcat(rows{:,i});    
    end                            

    data = rown(2:end);              

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
    dst = wl_dataQC(dst);  %generic data quality control for water levels
    output = {dst,datasetname};
end



