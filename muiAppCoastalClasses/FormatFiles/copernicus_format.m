function output = copernicus_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   copernicus_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   Copernicus wave reanalysis data
% USAGE
%   output = copernicus_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   ouput - updated class object
% NOTES
%   Full set of variables for the Atlantic- European North West Shelf- 
%   Wave Physics Reanalysis output includes:
%   VHM0 - Sea surface wave significant height[m]
%   VHM0_SW1 - Sea surface primary swell wave significant height [m]
%   VHM0_SW2 - Sea surface secondary swell wave significant height [m]
%   VHM0_WW - Sea surface wind wave significant height [m]
%   VMDR - Sea surface wave from direction [°]
%   VMDR_SW1 - Sea surface primary swell wave from direction [°]
%   VMDR_SW2 - Sea surface secondary swell wave from direction [°]
%   VMDR_WW - Sea surface wind wave from direction[°]
%   VPED - Sea surface wave from direction at variance spectral density maximum [°]
%   VSDX - Sea surface wave stokes drift x velocity [m/s]
%   VSDY - Sea surface wave stokes drift y velocity [m/s]
%   VTM01_SW1 - Sea surface primary swell wave mean period [s]
%   VTM01_SW2 - Sea surface secondary swell wave mean period [s]
%   VTM01_WW - Sea surface wind wave mean period [s]
%   VTM02 - Sea surface wave mean period from variance spectral density second frequency moment [s]
%   VTM10 - Sea surface wave mean period from variance spectral density inverse frequency moment [s]
%   VTPK - Sea surface wave period at variance spectral density maximum [s]
% SEE ALSO
%   Atlantic- European North West Shelf- Wave Physics Reanalysis
%   https://data.marine.copernicus.eu/product/NWSHELF_REANALYSIS_WAV_004_015/description
%   files are in NetCDF format
%
% Author: Ian Townend
% CoastalSEA (c)Dec 2024
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
    %get the datetime and variables from the netcdf file
    dur = ncread(filename,'time');
    date = datetime('1970-01-01 00:00:00');
    myDatetime =  date+seconds(dur);
    myDatetime.Format = 'dd-MM-yyyy HH:mm:ss';
    
    info = ncinfo(filename);
    varnames = {info.Variables(:).Name};
    %first 3 variables are time, latitude and longitude
    nvar = length(varnames)-3;

    lat = ncread(filename,'latitude');
    long = ncread(filename,'longitude');
    
    if length(lat)>1   %more than one point in file user will need to select a position
        %not yet implemented ****
        warndlg('More than one location in file. Subselection not yet implemented')
        dst = []; return;
    end
    
    varData = cell(1,nvar);
    for i=1:nvar
        varData{i} = squeeze(ncread(filename,varnames{i+3}));
    end

    %set metadata
    dsp = setDSproperties(varnames);
    %load the results into a dstable  
    dst = dstable(varData{:},'RowNames',myDatetime,'DSproperties',dsp); 
    dst.Dimensions.Position = [lat,long];   
end
%%
%--------------------------------------------------------------------------
% dataDSproperties
%--------------------------------------------------------------------------
function dsp = setDSproperties(varnames)
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'VHM0','VHM0_SW1','VHM0_SW2','VHM0_WW','VMDR','VMDR_SW1',...
                'VMDR_SW2','VMDR_WW','VPED','VSDX','VSDY','VTM01_SW1',...
                'VTM01_SW2','VTM01_WW','VTM02','VTM10','VTPK'},...              
        'Description',{'Sea surface wave significant height',...
                       'Sea surface primary swell wave significant height',...
                       'Sea surface secondary swell wave significant height',...
                       'Sea surface wind wave significant height',...
                       'Sea surface wave from direction',...
                       'Sea surface primary swell wave from direction',...
                       'Sea surface secondary swell wave from direction',...
                       'Sea surface wind wave from direction',...
                       'Sea surface wave from direction at variance spectral density maximum',...
                       'Sea surface wave stokes drift x velocity',...
                       'Sea surface wave stokes drift y velocity',...
                       'Sea surface primary swell wave mean period',...
                       'Sea surface secondary swell wave mean period',...
                       'Sea surface wind wave mean period',...
                       'Sea surface wave mean period from variance spectral density second frequency moment',...
                       'Sea surface wave mean period from variance spectral density inverse frequency moment',...
                       'Sea surface wave period at variance spectral density maximum'},...
        'Unit',{'m','m','m','m','deg','deg','deg','deg','deg','m/s','m/s',...
                's','s','s','s','s','s'},...
        'Label',{'Wave height (m)',...
                 'Wave height (m)',...
                 'Wave height (m)',...
                 'Wave height (m)',...
                 'Wave direction (deg)',...
                 'Wave direction (deg)',...
                 'Wave direction (deg)',...
                 'Wave direction (deg)',...
                 'Wave direction (deg)',...
                 'Velocity (m/s)',...
                 'Velocity (m/s)',...
                 'Wave period (s)',...
                 'Wave period (s)',...
                 'Wave period (s)',...
                 'Wave period (s)',...
                 'Wave period (s)',...
                 'Wave period (s)'},...
        'QCflag',repmat({'raw'},1,17)); 
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
    
    %now remove the variables that are not being loaded
    idx = ismember({dsp.Variables(:).Name},varnames);
    dsp.Variables(~idx) = [];

    %rename variables to be compatible with wave model which uses Hs, Tp
    %and Dir
    %replaced with clean-up function to subsample and rename variables.
%     modelvars = {'Hs','Dir','Tp'};
%     idx = ismember({dsp.Variables(:).Name},{'VHM0','VMDR','VTPK'});
%     [dsp.Variables(idx).Name] = modelvars{:};
end
%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
function output = dataQC(obj)
    %quality control a user data timeseries    
    datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    dst = obj.Data.(datasetname);      %selected dstable
    dst = wave_dataQC(dst);            %generic data quality control for waves
    output = {dst,datasetname};
end
