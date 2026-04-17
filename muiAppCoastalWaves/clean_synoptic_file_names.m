function clean_synoptic_file_names(outDir)
    %some spectra files have a non-standard format (eg Russington in 2004)
    %correct filename to standard format replacing space with a T    
    %This replaces timestamps like "2004-01-01 00h22"    
    %with timestamps like "2004-01-01T00h22" 
    % see also: extract_synoptic_from_zip.m
    [fname,path,nfiles] = getfiles('MultiSelect','on',...
                'FileType',{'*.spt'},'PromptText','Select zip file:');
    if nfiles==0, return; end

    outDir = [path,outDir];    
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    for i=1:nfiles
        ids = regexp(fname{i},'}');
        newfname = [fname{i}(1:ids+10),'T',fname{i}(ids+12:end)];
        copyfile([path,fname{i}], fullfile(outDir,newfname));
    end
    fprintf('Cleaned %d synoptic-nearest files to %s\n', nfiles, outDir);
end