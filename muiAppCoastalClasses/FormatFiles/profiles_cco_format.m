function output = profiles_cco_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   profiles_cco_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for:
%   CCO beach profile data
% USAGE
%   output = profiles_cco_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   Channel Coastal Observatory (CCO) data
%   https://www.channelcoast.org/
% SEE ALSO
%   clean_cco_bpfiles.m (now incorparated in getTSdata)
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
    [data,~] = readInputData(filename);             
    if isempty(data), dst = []; return; end
    %use readTable with auto detect import options 
%     opts = detectImportOptions(obj.filename,'FileType','text');  %v2016b
%     data = readtable(obj.filename,opts); 

    %set metadata
    dsp = setDSproperties;

    % get date of survey    
    ids = find(filename==filesep,1,'last'); 
    fileonlyname = filename(ids+1:end);
    idx = find(fileonlyname=='_');
    f_date = str2double(fileonlyname(idx+1:idx+8));
    f_date = datetime(f_date,'ConvertFrom','yyyymmdd');
    s_date = datetime(19000101,'ConvertFrom','yyyymmdd');
    c_date = datetime('now');
    %
    while (f_date>s_date && f_date<=c_date)==0
        prompt = {'Date of survey:'};
        dlg_title = 'Input';
        num_lines = 1;
        def = {num2str(yyyymmdd(f_date))};
        f_date = inputdlg(prompt,dlg_title,num_lines,def);
        if isempty(f_date)
            fprintf('Cancel and exit program \n');
            return
        end
        f_date = datetime(str2double(f_date),'ConvertFrom','yyyymmdd');
    end 
    mdat = f_date;
    mdat.Format = dsp.Row.Format;
    myDatetime = mdat + hours(12);
    
    %profile id and feature codes
    pid = data.Reg_ID;    
    prof = data.Profile;
    pfc = upper(data.FC);
    
    SFC = SFcodes;        %function file lists CCO codes
    codeset = SFC(:,1);
    valueset = SFC(:,3);
    FCunknown = 'Codes:';
    pFC = NaN(size(pfc));
    warnmsg = [];
    for j=1:length(pfc)
        FCid = find(strcmpi(pfc(j),codeset)); %case insensitive
        if isempty(FCid)
            if strcmpi(pfc(j),'B') %legacy code now replaced by BO
                pFC(j,1) = 10;                
            elseif all(cellfun(@length,pid)==0)
                %if profile ids not valid probably because FC column missing
                pFC(j,1) = 1; %set to zz = 'Unknown'
                warnmsg = sprintf('Beach profile FC column missing in file\n%s',...
                                                             filename);
                if all(cellfun(@length,pfc)~=0) &&  all(cellfun(@length,prof)==0)
                    %FC and Profiles columns missing and Reg_ID in FC column    
                    pFC(j,1) = 1;        %set to zz = 'Unknown
                    warnmsg = sprintf('Beach profile FC & Profile columns missing in file\n%s',...
                                                             filename); 
                end
            else
                %an unknown code identified                
                pFC(j,1) = 1; %set to zz = 'Unknown'
                warnmsg = sprintf('%s %s',warnmsg,pfc{j});                
            end
        else
            pFC(j,1) = valueset{FCid};
        end
    end
    %
    if ~isempty(warnmsg)
        if ~contains(warnmsg,'Beach profile FC')
            %attempt to correct for invalid feature code
            warnmsg = sprintf('Beach profile FC column contains an invalid entry:\n ''%s'' in file %s\nSet to ZZ = Unknown',...
                                                    warnmsg,filename);
        elseif contains(warnmsg,'Beach profile FC & Profile columns missing in file')
            %attempt to correct for missing columns
            pid = pfc;  
        end
        warndlg(warnmsg);
    end
    data.FC = pFC;
    %
    if ~isempty(pid)
        %prof_id is used as a variable name and so can only use
        %ascii characters A-Z, a-z, 0-9 and start with a letter
        if all(cellfun(@length,pid)==0)
            pid = prof;       %FC code probably missing so columsn shifted 1
            if all(cellfun(@length,pid)==0)  %check if this is solution
                warndlg(sprintf('Invalid data in file: %s NOT loaded \n',...
                    filename));
                return;
            end
        end
        if any(cellfun(@isempty,pid))
            idx = cellfun(@isempty,pid);
            pid(idx) = prof(idx);
            if any(cellfun(@isempty,pid)) && all(cellfun(@length,pid)>5)
                 warndlg(sprintf('Invalid data in file: %s NOT loaded \n',...
                        filename));
             return;
            end
         end
        
        newpid = strrep(pid, '-', 'm'); %replace - with m
        newpid = strrep(newpid, '+', 'p'); %replace + with p
        newpid = strcat('P',newpid);
        %Construct valid MATLAB identifiers from input strings
        cleanpid = matlab.lang.makeValidName(newpid,...
            'ReplacementStyle','delete');
        [prof_id,~] = unique(cleanpid);
        nprof = length(prof_id);
        
        %create a matrix of profile id and data of survey for each variable
        for j=1:nprof
            idx = find(strcmp(cleanpid,prof_id{j}));
            % information on data location
            Position(1) = data{idx(1),1};
            Position(2) = data{idx(1),2};
%             %some CCO data profiles have a 0,0 point at the end of the profile          
%             if data{idx,3}(end,1)==0 || data{idx,4}(end,1)==0
%                 data{idx,3}(1,end) = NaN;
%                 data{idx,4}(1,end) = NaN;
%             end
            %data to be loaded are for a single time step and need to be 
            %1xN vectors so that each data set is loaded as variable in tsc
            varData = mat2cell(data{idx,1:5},length(idx),[1,1,1,1,1]);
            varData = cellfun(@ctranspose,varData,'UniformOutput',false);
            
            %load the results into a dstable  
            adst = dstable(varData{:},'RowNames',myDatetime,'DSproperties',dsp); 
            %assign baseline point
            adst.UserData = Position;
%             adst.Dimensions.E = Position(1);
%             adst.Dimensions.N = Position(2);
            %assign site and profile_id as description
            adst.Description = prof_id{j};
            %add to stsc as a structure with field names using profile id
            dst.(adst.Description) = adst;
            clear adst
        end
    else
        warndlg(sprintf('Profile IDs not read: %s NOT loaded \n',...
            filename));
        dst = [];
    end  
end
%%
function [data,header] = readInputData(filename)
    %read wind data (read format is file specific).
    fid = fopen(filename, 'r');
    if fid<0
        errordlg('Could not open file for reading','File read error','modal')
        data = [];
        return;
    end            
    %find number of columns in file
    header = fgets(fid); %read first line of file
    delimiter = sprintf('\t'); %tab    
    ncols = numel(strfind(header,delimiter)) + 1;
    header = split(header,delimiter);
    header = matlab.lang.makeValidName(header,'ReplacementStyle','delete');    
    %
    %read data
    if ncols==7
        %Easting Northing Elevation_OD Chainage FC Profile Reg_ID
        nnum = 4;
        dataSpec = '%f %f %f %f %s %s %s';
    elseif ncols==5
        %Easting Northing Elevation_OD FC Name
        % nnum = 5;
        % dataSpec = '%f %f %f %s %s';
        errmsg = sprintf('Chainage missing in %s',...
                            filename);
        errordlg(errmsg,'Data error','modal')
        data = [];
        fclose(fid);
        return;
    elseif ncols==4
        %Easting Northing Elevation_OD FC
        %dataSpec = '%f %f %f %s';        
        errmsg = sprintf('Chainage and profile ID missing in %s',...
                            filename);
        errordlg(errmsg,'Data error','modal')
        data = [];
        fclose(fid);
        return;
    else
        errmsg = sprintf('Unknown file format in %s',filename);                                    
        errordlg(errmsg,'File read error','modal')
        data = [];
        fclose(fid);
        return;
    end
    data = textscan(fid,dataSpec); 
    datalen = cellfun(@length,data,'UniformOutput',false);
    isclean = sum(diff(cell2mat(datalen)))==0;
    if ~isclean
        frewind(fid);  %return to the beginning of the file
        data = cleanBPdata(fid);  
        if isempty(data)
            errmsg = sprintf('Unable to clean file %s',filename);                                    
            errordlg(errmsg,'File read error','modal')
            return;
        else
             warndlg(sprintf('Had to clean file:\n%s',filename));
        end   
    end  
    temp = num2cell(horzcat(data{:,1:nnum}));
    data = horzcat(temp,data{:,nnum+1:end});  
    data = cell2table(data,'VariableNames',header);
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
        'Name',{'Eastings','Northings','Elevation','Chainage','FeatureCode'},...                   % <<Edit metadata to suit model
        'Description',{'Eastings','Northings','Elevation','Chainage','FeatureCode'},...
        'Unit',{'m','m','mOD','m','-'},...
        'Label',{'Eastings (m)','Northings (m)','Elevation (mOD)','Chainage (m)','FeatureCode'},...
        'QCflag',repmat({'raw'},1,5));
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
%%
function data = cleanBPdata(fid)
    %remove additional field that occurs in some rows of CCO beach profile
    %data. See also: clean_CCO_bpfiles.m
    data = []; %#ok<NASGU>
    inpdata = textscan(fid,'%s','Delimiter','\n');
    txtdata = inpdata{1};
    delimiter = sprintf('\t'); %tab 
    %find number of columns in header
    ncols = numel(split(txtdata{1},delimiter));
    data = cell(1,ncols);
    for i=2:length(txtdata)
        %check whether each row as correct number of columns
        rowtxt = split(txtdata{i});
        ndata = numel(rowtxt);
        if ndata~=ncols
            %check that error is likely to be extra column at ncols
            if ndata==ncols+1
                %remove extra column - penultimate column of input row
                rowtxt(ncols) = [];
                ndata = numel(rowtxt);
                if ndata~=ncols, return; end
            else
                %a different error                
                return;
            end
        end
        %
        for j=1:ncols
            if j<5
                data{j}(i-1,1) = str2double(rowtxt{j});
            else
                data{j}{i-1,1} = rowtxt{j};
            end
        end
    end
end