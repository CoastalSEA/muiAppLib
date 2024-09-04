function output = wave_cco_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wave_cco_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   CCO wave data format
% USAGE
%   obj = wave_cco_format(funcall,varargin)
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
function dsp = setDSproperties()
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'Hs','Hmax','Tp','Tz','Dir','Spr','SST','Flag'},...
        'Description',{'Significant wave height','Maximum wave height',...
                'Peak period','Mean zero-crossing period','Wave direction',...
                'Spread','Sea surface temperature','Flag'},...
        'Unit',{'m','m','s','s','deg','deg','deg C','-'},...
        'Label',{'Wave height (m)','Wave height (m)','Wave period (s)',...
                 'Wave period (s)','Wave direction (deg)',...
                 'Wave spread (deg)','Sea surface temperature (deg)',...
                 'Error flag'},...
        'QCflag',repmat({'raw'},1,8)); 
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
% getData - option 1 using
%--------------------------------------------------------------------------
function dst = getData(obj,filename)
    %read data from file 
    %use readTable with auto detect import options  
    data = readTSinputFile(obj,filename);           
    if isempty(data), dst = []; return; end
    
    %some locales mean that first column is loaded as a cell array of char
    %hence 12 columns instead of 11
    if width(data)==12 && iscell(data{1,1})
        date = datetime(data.Var1,'InputFormat','dd-MMM-yyyy',...
                                            'Locale','en_GB'); %reformat date
        date = date+data.Var2;
        data = removevars(data,'Var1');                         %remove extra column
        data.Var2 = date;                                       %assign datetime
        varnames = {'Date_Time_GMT_','Latitude','Longitude','Flag','Hs_Hm0__m_','Hmax_m_','Tp_s_','Tz_Tm__s_','Dirp_degrees_','Spread_deg_','SST_degC_'};
        data.Properties.VariableNames = varnames;               %assign variables so that movevars works
    end    
    
    %set metadata
    dsp = setDSproperties;

    %extract required subset of data
    data = movevars(data,'Flag','After',11);  %move Flag to end of list
    varData = data(:,4:end);
    myDatetime = data{:,1};
    myDatetime.Format = dsp.Row.Format;

    % information on data location
    Lat  = data.Latitude;
    Long = data.Longitude;
    Lat(Lat==99 | Lat==999 | Lat==9999) = NaN;
    Latitude = mean(Lat,'omitnan');
    Long(Long==99 | Long==999 | Long==9999) = NaN;
    Longitude = mean(Long,'omitnan');

    %load the results into a dstable  
    dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
    dst.Dimensions.Position = [Latitude,Longitude];
end
%%
%--------------------------------------------------------------------------
% getData option 2 using readinputfile function
%--------------------------------------------------------------------------
% function dst = getData(~,filename)
%     %read and load a data set from a file
%     [data,~] = readInputData(filename);             
%     if isempty(data), dst = []; return; end
% 
%     %set metadata
%     dsp = setDSproperties;
%     
%     % concatenate date and time
%     mdat = data{1};
%     mtim = data{2};
%     mdat.Format = dsp.Row.Format;
%     mtim.Format = dsp.Row.Format;
%     myDatetime = mdat + timeofday(mtim);
% 
%     %move flag to end of list     
%     temp = data{5};
%     data(:,5) = [];
%     data{12} = temp;  
% 
%     %check that datetime values are unique
%     [UniqueTime,iu] = unique(myDatetime);
%     if length(UniqueTime)~=length(myDatetime)
%         myDatetime = UniqueTime;                
%     end
%     %check for missing data
%     varData = table(data{1,5:end});
%     varData = standardizeMissing(varData,[99,99.9,99.99,999,9999]);
%     varData = varData(iu,:);
%     
%     %load the results into a dstable  
%     dst = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp); 
% %     dst.Dimensions.Position = [Latitude,Longitude];
% end
% %%
% function [data,header] = readInputData(filename)
%     %read wind data (read format is file specific).
%     dataSpec = '%{dd-MMM-yyyy}D %{HH:mm:ss}D %f %f %d %f %f %f %f %f %f %f';
%     nhead = 1;
%     [data,header] = readinputfile(filename,nhead,dataSpec); %see muifunctions
% end
%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
function output = dataQC(obj)
    %quality control a dataset
    datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    dst = obj.Data.(datasetname);      %selected dstable
    
    maxTp = 25;
    %get input parameters from user
    prompt = {'Maximum value of Wave Peak Period (s):'};
    title = 'Define limiting peak period';
    numlines = 1;
    defaultvalues{1} = num2str(maxTp);      
    useInp=inputdlg(prompt,title,numlines,defaultvalues);
    if isempty(useInp), return; end %user cancelled
    maxTp = str2double(useInp{1});

    tol=0.1;
    idx = dst.Hs-tol > dst.Hmax;
    idx = idx | (dst.Tp>maxTp);
    %originally used Hs as 1/18 and Hmax as 1/16 but changed to
    %single wave limit of 1/7 = 0.14 for both
    Steep = dst.Hs./(1.56*dst.Tz.^2);  %denominator is gT^2/2pi
    idx = idx | (Steep>0.14);
    Steep = dst.Hmax./(1.56*dst.Tp.^2);  %denominator is gT^2/2pi
    idx = idx | (Steep>0.14);

    idx = idx | (dst.Tp<=0);   %Negative wave periods

    %CCO checks for Hs, Tz and Tp
    % Flag = 0 all data pass, 
    % Flag = 1 either Hs or Tz fail, so all data fail, 
    % Flag = 2 Tp fail + derivatives, 
    % Flag = 3 Direction fail + derivatives,  
    % Flag = 4  Spread fail + derivatives, 
    % Flag = 5 Tp fail Jump test only, (exclude on advice from Travis Morgan at CCO)
    % Flag = 7 Buoy adrift, 
    % Flag = 8 Sea Temperature fail, 
    % Flag = 9 Missing data (already assigned as NaN)
    idx = idx | (dst.Flag>0 & dst.Flag<5); 

    dst.Dir(dst.Dir<0) = dst.Dir(dst.Dir<0)+360;
    dst.Dir(dst.Dir>360) = dst.Dir(dst.Dir>360)-360;
    
    hw = waitbar(0, 'Loading data. Please wait');
    
    dst.Hs(idx) = NaN;    dst.Hmax(idx) = NaN;
    dst.Tz(idx) = NaN;    dst.Tp(idx) = NaN;
    dst.Dir(idx) = NaN;

    dst.VariableQCflags(1:5) = repmat({'qc'},1,5);
    
    output = {dst,datasetname};
    waitbar(1); 
    close(hw);
end



