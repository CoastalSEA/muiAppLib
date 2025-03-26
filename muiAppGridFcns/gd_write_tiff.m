function gd_write_tiff(filename,cdata,im)
%
%-------function help------------------------------------------------------
% NAME
%   gd_write_tiff.m
% PURPOSE
%   write the data from an image to a tif file including the position info
% USAGE
%    gd_write_tiff(filename,cdata,im)
% INPUTS
%   filename - full path and file name to tif file to write
%   cdata - image data MxNx3 RGB array of image
%   im - image struct array with limits for x and y and colormap used
% OUTPUT
%   tif file
% SEE ALSO
%   called from GDinterface.gridImage. gd_geoimage_format.m has a function
%   to read the file and produce an image struct
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2025
%--------------------------------------------------------------------------
%

    nrows = size(cdata, 1);
    ncols = size(cdata, 2);
    dx = (im.XData(2)-im.XData(1))/ncols;
    dy = (im.YData(2)-im.YData(1))/nrows;
    
    %colormap has to be 256x3
    % Create an array of indices for the original and new colormaps
    idx = linspace(1, size(im.CMap, 1), size(im.CMap, 1));
    newIndices = linspace(1, size(im.CMap, 1), 256);    
    % Interpolate each color channel (R, G, B) separately
    CMap = zeros(256, 3);
    for i = 1:3
        CMap(:, i) = interp1(idx, im.CMap(:, i), newIndices);
    end

    t = Tiff(filename, 'w');    
    % Set the tags for the TIFF file
    tagstruct.ImageLength = nrows;
    tagstruct.ImageWidth = ncols;
    tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = size(cdata, 3);
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    tagstruct.XPosition = im.XData(1)-dx/2;
    tagstruct.YPosition = im.YData(1)-dy/2;
    tagstruct.XResolution = dx;
    tagstruct.YResolution = dy;
    tagstruct.ColorMap = CMap;
    tagstruct.DocumentName = sprintf('%f %f',im.CLim(1),im.CLim(2));
    % Set the tags
    t.setTag(tagstruct);    
    % Write the image data
    t.write(cdata);
    % Close the Tiff object
    t.close();    
end