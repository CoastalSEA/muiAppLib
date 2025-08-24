function gd_write_image(filename,im)
%
%-------function help------------------------------------------------------
% NAME
%   gd_write_image.m
% PURPOSE
%   write image to jpg or tif file
% USAGE
%    gd_write_image(filename,im)
% INPUTS
%   filename - full path and file name to tif file to write
%   im - image struct array with limits for x and y and colormap used
%        struct('Location', desc,'XData',xLim,'YData',yLim,'CData',...
%                                      C.CData,'CMap',cmap,'CLim',cLimits)
% OUTPUT
%   image file written as image file depending on the file extension
% SEE ALSO
%   called from GDinterface.gridImage
%
% Author: Ian Townend
% CoastalSEA (c) Aug 2025
%--------------------------------------------------------------------------
%
    hf = figure('Visible','off');
    ax = axes(hf);
    axis tight
    axis off
    ax.Units = 'pixels';
    h_im = imagesc(ax,'XData',im.XData,'YData',im.YData,'CData',im.CData);
    if any(isnan(im.CData))
        set(h_im, 'AlphaData', 1-isnan(im.CData)); %set Nan values to be transparent 
    end
    colormap(im.CMap);
    clim(im.CLim);
    frame = getframe(hf,ax.Position);
    img = frame2im(frame);
    imwrite(img,filename);
    delete(hf)
end