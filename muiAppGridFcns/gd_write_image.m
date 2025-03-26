function gd_write_image(filename,desc,im)
%
%-------function help------------------------------------------------------
% NAME
%   gd_write_image.m
% PURPOSE
%   write cdata from an image to an ASCII text file, with the position info
% USAGE
%    gd_write_image(filename,desc,im)
% INPUTS
%   filename - full path and file name to tif file to write
%   desc - description of image
%   im - image struct array with limits for x and y and colormap used
%        struct('Location', desc,'XData',xLim,'YData',yLim,'CData',...
%                                      C.CData,'CMap',cmap,'CLim',cLimits)
% OUTPUT
%   image file written as ASCII text with '.txt' extension
% SEE ALSO
%   called from GDinterface.gridImage
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2025
%--------------------------------------------------------------------------
%
    fid = fopen(filename,'w'); 
    %write header
    projdate = cellstr(datetime('now'),'dd-MM-yyyy');
    fprintf(fid,'muiGrid image file\nDate: %s\n',projdate{1});
    fprintf(fid,'Location: %s\n',desc);
    fprintf(fid,'XData: %g %g\n',im.XData(1),im.XData(2));
    fprintf(fid,'YData: %g %g\n',im.YData(1),im.YData(2));
    fprintf(fid,'CLim: %g %g\n',im.CLim(1),im.CLim(2));
    fprintf(fid,'CData\n'); 
    %write data
    [rows, cols] = size(im.CData);
    for i = 1:rows
        for j = 1:cols
            fprintf(fid,'%g ',im.CData(i, j));
        end
        fprintf(fid, '\n');
    end
    %write color map
    fprintf(fid,'CMap\n');
    [rows, cols] = size(im.CMap);
    for i = 1:rows
        for j = 1:cols
            fprintf(fid,'%g ',im.CMap(i, j));
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end