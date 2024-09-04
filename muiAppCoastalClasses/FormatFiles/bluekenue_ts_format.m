function output = bluekenue_ts_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wave_cco_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   BlueKenue timeseries data format
% USAGE
%   obj = wave_cco_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   https://nrc.canada.ca/en/research-development/products-services/software-applications/blue-kenuetm-software-tool-hydraulic-modellers
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
% getData
%--------------------------------------------------------------------------
function dst = getData(~,filename)
    %load BlueKenue timeseries data based file formats ts1-ts4
    [header,data] = read_data(filename);  %read data from file
    if isempty(data), return; end
    % unpack data for time variable
    nvar = length(data);    %number of columns in data set (cells)
    nrec = length(data{1}); %number of records in each column
    switch header.filetype
        case {'ts1','ts2'}
            myDatetime = datetime(zeros(nrec,6));
            myDatetime.Format = 'dd-MMM-uuuu HH:mm:ss';                    
            myDatetime(1) = header.StartTime;
            dt = header.DeltaT-datetime('0:00:00');
            myDatetime(2:nrec) = myDatetime(1)+dt*(1:nrec-1);
            offset = 0;  %number of variables in 'data' already unpacked
        case {'ts3','ts4'}
            mdat = data{1};
            if ~isdatetime(mdat)
                myDuration = getTimeData(mdat);
                startyear = sprintf('01-Jan-%d 00:00:00',0);            
                stdat = datetime(startyear);
                myDatetime = stdat + myDuration;
                offset = 1;  %number of variables in 'data' already unpacked
            else
                % concatenate date and time
                mdat = data{1};
                mtim = data{2};
                mdat.Format = 'dd-MMM-uuuu HH:mm:ss';
                mtim.Format = 'dd-MMM-uuuu HH:mm:ss';
                myDatetime = mdat + timeofday(mtim);
                offset = 2;  %number of variables in 'data' already unpacked
            end
    end

    varData = data(offset+1:nvar);
        
    
    dsp = setDSproperties(header);
    %load the table varData into a dstable  
    dst = dstable(varData{:},'RowNames',myDatetime,'DSproperties',dsp); 
    
    % information on data location
    Xposition = header.location(1);
    Yposition = header.location(2);  
    if Xposition==Yposition
        %add small offset to ensure unique values
        %alternative is to load x and y as separate dimensions
        Yposition = Yposition+0.1;
    end
    dst.Dimensions.Position = [Xposition,Yposition];    
end
%%
function [header,data] = read_data(filename)
    %read BlueKenue timeseries data (read format is file specific).
    fid = fopen(filename, 'r');
    if fid<0
        errordlg('Could not open file for reading','File write error','modal')
        header = ''; data = [];
        return;
    end
    %read BlueKenue ts header
    tline = fgetl(fid); hflag = true;
    XYposition = [];
    h.vardesc = 'No Name';
    while ischar(tline) && hflag %read header
        disp(tline)
        tline = fgetl(fid);                
        if strcmp(tline(1),':')  %ordered by length of definition string             
            if strcmpi(tline(2:5),'Name')
                k1 = strfind(tline,'(X');
                k2 = strfind(tline,')'); 
                if ~isempty(k1) && ~isempty(k2) %coordinates included
                    temp = textscan(tline(k1+1:k2(k2>k1)-1),'%s %f %s %f');
                    XYposition(1,1) = temp{2};
                    XYposition(1,2) = temp{4};
                    h.vardesc = tline(6:k1-1);
                else                            %name only
                    h.vardesc = tline(6:end);
                end
            elseif strcmpi(tline(2:7),'DeltaT')      
                temp = textscan(tline,'%s %{HH:mm:ss.SSS}D');
                h.DeltaT = temp{2};    
            elseif strcmpi(tline(2:10),'LocationX')
                temp = textscan(tline,'%s %f');
                XYlocation(1,1) = temp{1,2};
            elseif strcmpi(tline(2:10),'LocationY')    
                temp = textscan(tline,'%s %f');
                XYlocation(1,2) = temp{1,2};
            elseif strcmpi(tline(2:10),'StartTime') 
                temp = textscan(tline,'%s %{yyyy/MM/dd}D %{HH:mm:ss.SSS}D');
                mdat = temp{2};
                mdat.Format = 'dd-MMM-uuuu HH:mm:ss';
                mtim = temp{3};
                mtim.Format = 'dd-MMM-uuuu HH:mm:ss';
                h.StartTime = mdat + timeofday(mtim);                    
            elseif strcmp(tline(2:10),'EndHeader')
                hflag = false;    
            elseif strcmpi(tline(2:9),'FileType')
                h.filetype = tline(11:13);    
            elseif strcmpi(tline(2:15),'AttributeUnits') 
                temp = textscan(tline,'%s %f %s');
                h.units = lower(temp{3});
            elseif strcmpi(tline(2:15),'DataDefinition') 
                temp = textscan(tline,'%s %s');
                if strcmpi(temp{2},'MAGDIR')
                    h.units{2} = 'degrees';   
                else
                    h.units{2} = 'm/s';
                end                    
            end
        end
    end
    % select which location definition to use
    if XYlocation(1)==0 && XYlocation(2)==0
        if ~isempty(XYposition)
            h.location = XYposition;
        else
            h.location = [0,0];
        end
    else
        h.location = XYlocation;
    end
    %read numeric data
    switch h.filetype
        case {'ts1','ts3'}
            dataSpec = '%f';
        case {'ts2','ts4'}
            dataSpec = '%f %f';
    end 

    if strcmp(h.filetype,'ts3') || strcmp(h.filetype,'ts4')
        dataSpec = getDateTimeFormat(fid,dataSpec);
    end

    data = textscan(fid,dataSpec);
    header = h;
    fclose(fid);
end  
%%
%%
function dataspec = getDateTimeFormat(fid,dataspec)
    %read the first line of data to determin date format
    ndata = numel(split(dataspec)); %number of data columns
    filepos = ftell(fid); %find current position in file
    tline = fgetl(fid);   %read first line of data
    ndatim = numel(split(tline))-ndata; %number of parts to data-time        
    if ndatim>1
        %if two parts then this should be date and time
        dataspec = ['%{yyyy/MM/dd}D %{HH:mm:ss.SSS}D ',dataspec];
    else
        %if one part this should be elapsed time in HH:mm:ss.SSS
        dataspec = ['%q ',dataspec];
    end
    fseek(fid,filepos,'bof'); %return file pointer to first line of data
end
%%
function timedata = getTimeData(timecells)
    %convert cell time strings to duration data
    hms = cell2table(split(timecells,':')); %create hh,mm,ss table
    numhms = zeros(size(hms));
    for i=1:3 %convert hh,mm,ss characters to numerical values
        var = ['Var',num2str(i)];
        alphms = string(hms.(var));
        numhms(:,i) = str2double(alphms);
    end
    timedata = duration(numhms);   %return hh:mm:ss as a duration
end
%%
%--------------------------------------------------------------------------
% dataDSproperties
%--------------------------------------------------------------------------
function dsp = setDSproperties(header)
    %define the variables in the dataset
        % unpack remaining variables
    if isempty(header.units)
        varUnit = {'-'};
    else
        varUnit = header.units;
    end            

    prompt = {'Variable name','Variable description:',...
              'Variable units','Variable Label'};
    title = 'Load BKts data';
    numlines = 1;
    vlabel = sprintf('%s (%s)',header.vardesc,varUnit{1});
    default = {'var1',header.vardesc,varUnit{1},vlabel};
    answer = inputdlg(prompt,title,numlines,default);
    
    varQC = {'raw'};
    switch header.filetype   
        case {'ts1','ts3'}
            varName = answer(1);
            varDesc = answer(2); 
            varUnit = answer(3);
            varLabels = answer(4);
        case {'ts2','ts4'}
            if strcmpi(header.units{2},'degrees')
                varName = {answer{1},'Dir'};
                varDesc = {answer{2},'Direction'}; 
                varUnit = {answer{3},'degrees'};
                varLabels = {answer(4),'Direction'};
            else
                varName = {'U','V'}; 
                varDesc = {['U ',answer{2}],['V ',answer{2}]};  
                varUnit = {answer{3},answer{3}};
                varLabels = {answer{4},answer{4}};
                varQC = {'raw','raw'};
            end
    end
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',varName,...            
        'Description',varDesc,...
        'Unit',varUnit,...
        'Label',varLabels,...
        'QCflag',varQC); 
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

