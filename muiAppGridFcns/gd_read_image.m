function im = gd_read_image(filename)
%
%-------function help------------------------------------------------------
% NAME
%   gd_read_image.m
% PURPOSE
%   read cdata for an image from an ASCII text file, with the position info
% USAGE
%    im = gd_read_image()
% INPUTS
%   filename - full path and file name to tif file to read
% OUTPUT
%   im - image struct array with limits for x and y and colormap used
%        struct('XData',xLim,'YData',yLim,'CData',C.CData,'CMap',cmap,'CLim',cLimits)   
% SEE ALSO
%   called from gd_geoimage_format to load image from file
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2025
%--------------------------------------------------------------------------
%
    %read header
    fid = fopen(filename,'r');
    txt = fgetl(fid);
    if ~contains(txt,'muiGrid image')
        warndlg('File is not a muiGrid image'); return;
    end
    %read Location, XData, YData, CLim
    for i=1:5
        aline = fgetl(fid);
        vartxt = split(aline,':');
        im.(vartxt{1}) = strip(vartxt{2});
        if i>2
            im.(vartxt{1}) = str2num(im.(vartxt{1})); %#ok<ST2NM> vector
        end
    end
    fgetl(fid);
    nlines = findTableLines(fid,{'CData','CMap'});
    opts = detectImportOptions(filename);
    opts.DataLines = [nlines(1)+1,nlines(2)-1];
    cdata = readmatrix(filename,opts);
    nrows = size(cdata,1);
    for i=1:nrows
        CData(i,:) = str2num(cdata{i,:}); %#ok<AGROW,ST2NM> 
    end
    im.CData = CData;
    opts.DataLines = [nlines(2)+1,inf];
    cmap = readmatrix(filename,opts);
    nrows = size(cmap,1);
    for i=1:nrows
        CMap(i,:) = str2num(cmap{i,:}); %#ok<ST2NM,AGROW> 
    end
    im.CMap = CMap;
    fclose(fid);
end

%%
function nlines = findTableLines(fid,propnames)
    %find the start and end lines of a property table
    frewind(fid);
    nlines = zeros(size(propnames));
    nline = 1;
    tline = fgetl(fid);
    while ~feof(fid)
        if ischar(tline) && any(matches(propnames,tline(1:end)))
            nlines(matches(propnames,tline(1:end))) = nline;
        end
        nline = nline+1;
        tline = fgetl(fid);
    end     
end

