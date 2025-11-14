function output = wave_cco_spectra(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wave_cco_spectra.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   CCO wave data format
% USAGE
%   obj = wave_cco_spectra(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   Channel Coastal Observatory (CCO) data
%   https://www.channelcoast.org/
%
% Author: Ian Townend
% CoastalSEA (c) March 2023
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
        case 'setDSproperties'
            [output.dspec,output.dsprop] = setDSproperties();
    end
end

%%
%--------------------------------------------------------------------------
% getData using readinputfile function
%--------------------------------------------------------------------------
function dst = getData(~,filename)
    %read and load a data set from a file
    [data,header] = readInputData(filename);             
    if isempty(data), dst = []; return; end

    if length(data{1})~=64
        getdialog(sprintf('Incorrect record length in file:\n%s',filename))
        dst = []; return;
    end

    %set metadata
    [dspectra,dsparams] = setDSproperties;
    
    % concatenate date and time from filename
    ids = regexp(filename,'}');
    idm = regexp(filename,'T');
    idm = idm(idm>ids);
    ide = regexp(filename,'Z.');
    mdat = datetime(filename(ids+1:idm-1));
    mdat.Format = dspectra.Row.Format;
    mtim = filename(idm+1:ide-1);
    mtim(3) = ':';
    mtim = datetime(mtim,'InputFormat','HH:mm');
    myDatetime = mdat + timeofday(mtim);

    idp = regexp(filename,'/');
    if isempty(idp)
        idp = regexp(filename,'\');
    end
    Location = filename(idp(end)+1:ids-1);

    %extract spectral data
    varData = data(2:end);    
    Smax = str2double(header{4});

    varData{1} = varData{1}*Smax;  %convert relative to absolute spectral energy
    varData = cellfun(@transpose,varData,'UniformOutput',false);

    %load the results into a dstable  
    dst.Spectra = dstable(varData{:},'RowNames',myDatetime,'DSproperties',dspectra); 
    dst.Spectra.Dimensions.freq = data{1};
    dst.Spectra.Description = Location;

    %add header information to a dstable
    header = cellfun(@str2double,header,'UniformOutput',false);
    input = header([2,3,4,6]);    %extract required variables
    input{1} = input{1}/100;      %convert wave height from cm to m
    dst.Properties =  dstable(input{:},'RowNames',myDatetime,'DSproperties',dsparams);
    dst.Properties.Description = Location;
end

%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    dataSpec = '%f,%f,%f,%f,%f,%f';
    nhead = 12;
    [data,header] = readinputfile(filename,nhead,dataSpec); %see muifunctions
end

%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
function output = dataQC(obj) %#ok<INUSD> 
    %quality control a dataset
    % datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    % dst = obj.Data.(datasetname);      %selected dstable
    warndlg('No quality control defined for this format');
    output = [];    %if no QC implemented in dataQC
end

%%
%--------------------------------------------------------------------------
% dataDSproperties
%--------------------------------------------------------------------------
function [dspec,dsprop] = setDSproperties()
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dspec = struct('Variables',[],'Row',[],'Dimensions',[]); dsprop = dspec;    
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dspec.Variables = struct(...
        'Name',{'S','Dir','Spr','Skew','Kurt'},...
        'Description',{'Spectral density','Direction',...
                'Directional spread','Skewness of the directional distribution',...
                'Kurtosis of the directional distribution'},...
        'Unit',{'m^2/Hz','deg','deg','-','-'},...
        'Label',{'Spectral density (m^2/Hz)','Direction (degTN)',...
                'Directional spread (deg)','Skewness of the directional distribution (-)',...
                'Kurtosis of the directional distribution (-)'},...
        'QCflag',repmat({'raw'},1,5)); 
    dspec.Row = struct(...
        'Name',{'Time'},...
        'Description',{'Time'},...
        'Unit',{'h'},...
        'Label',{'Time'},...
        'Format',{'dd-MM-yyyy HH:mm:ss'});        
    dspec.Dimensions = struct(...    
        'Name',{'freq'},...
        'Description',{'Frequency'},...
        'Unit',{'Hz'},...
        'Label',{'Frequency (Hz)'},...
        'Format',{''});   

    %struct for wave specctra header information
    % could include all parameters in header    
    %     tn: transmission index (1 to 8) 
    %     Hs: significant wave height [cm] 
    %     Tz: zero-upcross period [s] 
    %     Smax: maximum of the psd S(f) [m^2/Hz] 
    %     Tref: reference temperature [°C, centigrade] 
    %     Tsea: Sea surface temperature [°C, centigrade] 
    %     Bat: Battery status (0 = empty to 7 = full) 
    %     Av: offset of the vertical accelerometer 
    %     Ax: offset of the x-accelerometer 
    %     Ay: offset of the y-accelerometer 
    %     Ori: buoy orientation [°] 
    %     Incli: magnetic inclination [°] 
    % for now just use Hs, Tz, Smax and SST (Tsea).
    dsprop.Variables = struct(...
        'Name',{'Hs','Tz','Sp','SST'},...
        'Description',{'Significant wave height','Zero-upcrossing period ',...
                'Spectral energy at peak','Sea surface temperature'},...
        'Unit',{'m','s','m^2/Hz','degC'},...
        'Label',{'Wave height (m)','Wave period (s)',...
                'Spectral energy at peak (m^2/Hz)','Temperature (degC)'},...
        'QCflag',repmat({'raw'},1,4)); 
    dsprop.Row = struct(...
        'Name',{'Time'},...
        'Description',{'Time'},...
        'Unit',{'h'},...
        'Label',{'Time'},...
        'Format',{'dd-MM-yyyy HH:mm:ss'});        
    dsprop.Dimensions = struct(...    
        'Name',{''},...
        'Description',{''},...
        'Unit',{''},...
        'Label',{''},...
        'Format',{''});     
end

