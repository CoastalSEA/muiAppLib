function extract_synoptic_from_zip(outDir,tol)
%
%-------function help------------------------------------------------------
% NAME
%   extract_synoptic_from_zip.m
% PURPOSE
%   Extracts files closest to synoptic times inorde to subsample a set of
%   files. e.g. spt files from a wave buoy at half hour intervals to 
%   synoptic hours.
% USAGE
%   extract_synoptic_from_zip(outDir,tol);
% INPUTS
%   outDir - name of folder to write selected files to
%            OR
%            numeric vector of years (based on years in zip file names)
%   tol - tolerance relative to synoptic times to use for selection.
%   Optional, default is 0 s and no tolerance is used.
% OUTPUT
%   selected files written to outDir
% NOTES
%   used for Channel Coastal Observatory (CCO) data. https://www.channelcoast.org/
% SEE ASLO  
%   clean_synoptic_file_names(outDir) to convert files with timestamps of
%   "2004-01-01 00h22" to standard format of "2004-01-01T00h22" 
%
% Author: Ian Townend
% CoastalSEA (c)Feb 2021
%--------------------------------------------------------------------------
%
% extract_synoptic_from_zip  Extracts files closest to synoptic times
%   extract_synoptic_from_zip('input.zip', 'synoptic_out')
%   
    if nargin<2, tol = 0; end
    if ~isduration(tol), tol = seconds(tol); end   %make tol a duration in seconds 
   
    if isnumeric(outDir)
        [fname,path,nfiles] = getfiles('MultiSelect','on',...
                    'FileType',{'*.zip'},'PromptText','Select zip file:');
        if nfiles==0, return; end
        if numel(outDir)~=nfiles
            warndlg('Number of files selected does not match year inputs')
            return;
        end
        outDir = cellstr(num2str(outDir(:)));        
    else
        [fname,path,nfiles] = getfiles('MultiSelect','off',...
                    'FileType',{'*.zip'},'PromptText','Select zip file:');
        if nfiles==0, return; end
        fname = {fname}; outDir = {outDir};
    end

    for i=1:numel(fname)
        zipFile = [path,fname{i}];
        outDiri = [path,outDir{i}];

        if ~isfolder(outDiri)
            mkdir(outDiri);
        end
    
        % Synoptic hours
        synopticHours = [0 3 6 9 12 15 18 21];
    
        % List contents of ZIP
        listing = unzip(zipFile, tempname);  % extract to temp folder
    
        % Prepare storage
        fileInfo = struct('fname', {}, 'dt', {});
    
        % Regex for timestamps like "2004-01-01 00h22"    
        % expr = '(\d{4}-\d{2}-\d{2})[ _](\d{2})h(\d{2})';
        % Regex for timestamps like "2004-01-01T00h22"    
        expr = '(\d{4}-\d{2}-\d{2})[T](\d{2})h(\d{2})';
    
        % Parse timestamps
        for i = 1:numel(listing)
            [~, name, ext] = fileparts(listing{i});
            full = [name'.', ext];
    
            tokens = regexp(full, expr, 'tokens', 'once');
            if ~isempty(tokens)
                dateStr = tokens{1};
                hh = str2double(tokens{2});
                mm = str2double(tokens{3});
                dt = datetime(dateStr + " " + sprintf('%02d:%02d', hh, mm), ...
                              'InputFormat','yyyy-MM-dd HH:mm');
                fileInfo(end+1).fname = listing{i};
                fileInfo(end).dt = dt;
            end
        end
    
        if isempty(fileInfo)
            warning('No timestamped files found in ZIP.');
            return
        end
    
        % Group by date
        allDates = unique([fileInfo.dt].');
        allDates = unique(dateshift(allDates, 'start', 'day'));
    
        selected = {};
    
        for d = allDates.'
            % Files for this day
            idx = isbetween([fileInfo.dt], d, d+days(1));
            todays = fileInfo(idx);
    
            if isempty(todays)
                continue
            end
    
            for h = synopticHours
                target = d + hours(h);
    
                % Compute absolute time difference
                diffs = abs([todays.dt] - target);
    
                [mindiff, k] = min(diffs);
                if tol==0 || mindiff<tol
                    selected{end+1} = todays(k).fname; %#ok<AGROW>
                end
            end
        end
    
        % Copy selected files to output directory
        for i = 1:numel(selected)
            [~, base, ext] = fileparts(selected{i});
            copyfile(selected{i}, fullfile(outDiri, [base,'.', ext]));
        end
    
        fprintf('Extracted %d synoptic-nearest files to %s\n', numel(selected), outDiri);
    end
end
