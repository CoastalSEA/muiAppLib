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
%   output - dstables for the sprectrum data and properties held in
%            datasets named sptSpectrum and sptProperties
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
        case 'addData'
            output = addData(varargin{:});
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

    flim = [0.025,0.58]; %frequency limits
    f = [flim(1):0.005:0.1,0.11:0.01:flim(2)];     
    
    if numel(data{1})>64   %some files have more than 64 frequencies
        nvar = numel(data);
        % Compute pairwise absolute differences
        [~, idf] = min(abs(data{1}(:)' - f(:)), [], 2);
        d = cell(1,nvar);
        for j=1:nvar
            d{j} =data{j}(idf);
        end
        data = d;
    elseif numel(data{1})<64 %some have less
        msgbox(sprintf('Incorrect record length in file:\n%s',filename))
        dst = []; return;
    end

    %set metadata
    [dspectra,dsparams] = setDSproperties;
    [~, name, ext] = fileparts(filename);
    full = [name'.', ext];

    % concatenate date and time from filename
    % Regex for timestamps like "2004-01-01T00h22"    
    expr = '(\d{4}-\d{2}-\d{2})[T](\d{2})h(\d{2})';
    tokens = regexp(full, expr, 'tokens', 'once');
    if ~isempty(tokens)
        dateStr = tokens{1};
        hh = str2double(tokens{2});
        mm = str2double(tokens{3});
        myDatetime = datetime(dateStr + " " + sprintf('%02d:%02d', hh, mm), ...
                      'InputFormat','yyyy-MM-dd HH:mm');
        myDatetime.Format = 'dd-MM-yyyy HH:mm:ss';
    else
        dst = []; return;
    end

    %check and cleanup frequencies if duplicated
    ide = find(data{1}==flim(1));
    if numel(ide)>1
        data{1}(ide) = f(ide);
    end
    
    %get location
    idp = regexp(name,'}');
    Location = name(1:idp-1);

    %extract spectral data
    varData = data(2:end);    
    Smax = str2double(header{4});
    Sin = max(varData{1});
    if Sin<=1
        %some files have relative data others absolute values. If relative:
        varData{1} = varData{1}*Smax;  %convert relative to absolute spectral energy
    end
    varData = cellfun(@transpose,varData,'UniformOutput',false);

    %load the results into a dstable  
    dst.sptSpectrum = dstable(varData{:},'RowNames',myDatetime,'DSproperties',dspectra); 
    dst.sptSpectrum.Dimensions.freq = data{1};
    dst.sptSpectrum.Description = Location;

    %add header information to a dstable
    header = cellfun(@str2double,header,'UniformOutput',false);
    input = header([2,3,4,6]);    %extract required variables
    input{1} = input{1}/100;      %convert wave height from cm to m    
    dst.sptProperties =  dstable(input{:},'RowNames',myDatetime,'DSproperties',dsparams);
    dst.sptProperties.Description = Location;
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
function output = dataQC(obj)
    %quality control a dataset
    inp = inputdlg({'Limit for peak spectral density (m^2/Hz)'},'QC',1,{'100'});
    if isempty(inp); inp = {'100'}; end
    threshold = str2double(inp{1});

    hw = waitbar(0, 'Loading data. Please wait');
    dstspec = obj.Data.sptSpectrum;      %selected dstables
    dstprop = obj.Data.sptProperties;    
    
    %run checks based on wave steepness and peak spectral density
    Steep = 2*pi*dstprop.Hs./(9.81*dstprop.Tz.^2);  %Tz deep water wave steepness
    for i=1:numel(Steep)
        if Steep(i)>0.1 || dstprop.Sp(i)>threshold || dstprop.Tz(i)>30
            dstprop.Hs(i) = NaN;
            dstprop.Sp(i) = NaN;
            dstprop.Tz(i) = NaN;
            dstspec = removeRecord(dstspec,i);
        end
    end
    dstspec.VariableQCflags(1:5) = repmat({'qc'},1,5);
    dstprop.VariableQCflags(1:4) = repmat({'qc'},1,4);

    dst.sptSpectrum = dstspec;
    dst.sptProperties = dstprop;
    output = {dst};
    waitbar(1); 
    close(hw);
end

%%
function dst = removeRecord(dst,idx)
    %remove record for all variables
    dst.S(idx,:) = NaN;
    dst.Dir(idx,:) = NaN;
    dst.Spr(idx,:) = NaN;
    dst.Skew(idx,:) = NaN;
    dst.Kurt(idx,:) = NaN;
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

