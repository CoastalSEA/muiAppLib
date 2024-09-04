function output = profiles_chainage_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   profiles_chainage_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   Profiles with data chaninage and elevation format
%   
% USAGE
%   output = profiles_chainage_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   e.g. Duck and Narrabeen data sets
%   CT_Plots assumes that the x-values are named 'Chainage'
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
function dst = getData(~,filename)
    %read and load a data set from a file
    data = readInputData(filename);             
    if isempty(data), dst = []; return; end
    
    %set metadata
    dsp = setDSproperties;

    mdat = unique(data.Date);
    mdat.Format = dsp.Row.Format;
    myDatetime = mdat + hours(12);
    
    % unpack data to specific variables   
    site = data.Site;        %location
    pid = data.ProfileID;    %profile id
    pdat = data.Date;        %date of survey for each point
    
    %clean up and assign profile names
    site_id = unique(site);
    %prof_id is used as a variable name and so can only use
    %ascii characters A-Z, a-z, 0-9 and start with a letter
    newpid = strrep(pid, '-', 'm'); %replace - with m
    newpid = strrep(newpid, '+', 'p'); %replace + with p
    newpid = strcat('P',newpid);            
    %Construct valid MATLAB identifiers from input strings
    cleanpid = matlab.lang.makeValidName(newpid,...
                     'ReplacementStyle','delete');
    prof_id = unique(cleanpid);
    nprof = length(prof_id);
    
    for k=1:length(site_id)
        %some files only have one site but this provides for multiples
        for j=1:nprof
            varData = [];
            for i=1:length(myDatetime)
                idx = find(strcmp(cleanpid,prof_id{j}) & ...
                    pdat==mdat(i) & strcmp(site,site_id{k}));
                
                %data are for a single time step and need to be
                %1xN vectors so that each data to concatenate time steps
                profData = mat2cell(data{idx,4:5},length(idx),[1,1]);
                profData = cellfun(@ctranspose,profData,'UniformOutput',false);
                %
                if isempty(varData)
                    varData = profData;
                else
                    newlen = size(profData{1},2);
                    oldlen = size(varData{1,1},2);
                    nrec = size(varData,1);
                    nvar = 2;   %chainage and elevation are the variables
                    diff = abs(newlen-oldlen);
                    pad = ones(1,diff)*NaN;
                    %
                    if newlen>oldlen
                        for it = 1:nrec
                            for jv = 1:nvar
                                varData{it,jv}(1,oldlen+1:newlen) = pad;
                            end
                        end
                    elseif newlen<oldlen
                        for jv = 1:nvar
                            profData{1,jv}(1,newlen+1:oldlen) = pad;
                        end
                    end
                    varData = vertcat(varData,profData);
                end
            end
            varData= {cell2mat(varData(:,1)),cell2mat(varData(:,2))};

            %load the results into a dstable  
            adst = dstable(varData{:},'RowNames',myDatetime,'DSproperties',dsp); 
            %assign site and profile_id as description
            adst.Description = sprintf('%s-%s',site_id{k},prof_id{j});
            %add to stsc as a structure with field names using profile id
            dst.(adst.Description) = adst;
            clear adst
        end
    end
end
%%
function data = readInputData(filename)
    %read chainage beach profile data (read format is file specific)
    %currently reads Narrabeen and Duck beach profile data files with format:
    %Site ProfileID Date Chainage Elevation Flag
    fid = fopen(filename, 'r');
    if fid<0
        errordlg('Could not open file for reading','File read error','modal')
        data = [];
        return;
    end           
    %find number of columns in file
    header = fgets(fid); %read first line of file
    delimiter = sprintf(' '); %space
    ncols = numel(strfind(header,delimiter)) + 1;     
    header = split(header,delimiter);
    %
    %read data
    if ncols==6
        %Site ProfileID Date Chainage Elevation Flag
        dataSpec = '%s %s %{yyyy-MM-dd}D %f %f %s';
    else
        errmsg = sprintf('Unknown file format in %s',filename);                                    
        errordlg(errmsg,'File read error','modal')
        data = [];
        fclose(fid);
        return;
    end
    din = textscan(fid,dataSpec); 
    data = table(din{1},din{2},din{3},din{4},din{5},din{6},...
             'VariableNames',header);  
    fclose(fid);
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
        'Name',{'Chainage','Elevation'},...                   % <<Edit metadata to suit model
        'Description',{'Chainage','Elevation'},...
        'Unit',{'m','mOD'},...
        'Label',{'Chainage (m)','Elevation (mOD)'},...
        'QCflag',{'raw','raw'}); 
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
function output = dataQC(obj)
    %quality control a dataset
    % datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    % dst = obj.Data.(datasetname);      %selected dstable
    warndlg('No quality control defined for this format');
    output = [];    %if no QC implemented in dataQC
end
