function output = gd_geoimage_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   gd_geoimage_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for
%   image data
% USAGE
%   output = gd_geoimage_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   This file loads a geoimage created from a grid which is saved in a
%   dstable as a struct because imshow uses XLim and Ylin rather than x,y
%   vectors
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    switch funcall
        %standard calls from muiDataSet - do not change if data class 
        %inherits from muiDataSet. The function getPlot is called from the
        %Abstract method tabPlot. The class definition can use tabDefaultPlot
        %define plot function in the class file, or call getPlot
        case 'getFormat'
            output = getFormat(varargin{:});
        case 'getData'
          output = getData(varargin{:});
        case 'setData'
          output = setData(varargin{:});  
        case 'dataQC'
            output = dataQC(varargin{1});  
        case 'getPlot'
            %output = 0; if using the default tab plot in muiDataSet, else
            output = getPlot(varargin{:});
    end
end

%%
%--------------------------------------------------------------------------
% getFormat
%--------------------------------------------------------------------------
function obj = getFormat(obj,formatfile)
    %return the file import format settings
    obj.DataFormats = {'muiUserData',formatfile,'data'};
    obj.idFormat = 1;
    obj.FileSpec = {'on','*.mat;*.tif;*.tiff;*.txt;'};
end

%%
%--------------------------------------------------------------------------
% getData
%--------------------------------------------------------------------------
function newdst = getData(obj,filename,metatxt) %#ok<INUSD>
    %read and load a data set from a file
    dsp = setDSproperties;                 %set metadata
    [~,location,ext] = fileparts(filename);  
    if strcmp(ext,'.mat')
        imdata = {load(filename)};    
        %load the results into a dstable
        dst = dstable(imdata{1}.im,'RowNames',{location},'DSproperties',dsp); 
    elseif strcmp(ext,'.txt')
        im = gd_read_image(filename);
        %load the results into a dstable
        dst = dstable(im,'RowNames',{location},'DSproperties',dsp); 
    else
        im = readTiff(filename,ext(2:end));
        %load the results into a dstable
        dst = dstable(im,'RowNames',{location},'DSproperties',dsp); 
    end
    dst.Description = location;
    dst.Source = filename;
    dst.MetaData = metatxt;

    dst.UserData = [];         %unused
    newdst.GeoImage = dst;        %GeoImage is the dataset name for this format
end

%%
%--------------------------------------------------------------------------
% getData
%--------------------------------------------------------------------------
function newdst = setData(obj,imobj,metatxt) %#ok<INUSD>
    %read and load a data set from a file
    %imobj = struct('XData',[],'YData',[],'CData',[],'CMap',[],'CLim',[]);
    dsp = setDSproperties;                 %set metadata
    %load the results into a dstable - Image is the dataset name for this format
    dst = dstable(imobj,'RowNames',{metatxt},'DSproperties',dsp);  
    dst.Description = metatxt;
    dst.Source = 'GeoImage of grid plot';
    dst.MetaData = sprintf('Georeferenced image for a grid of %s',metatxt); 
    dst.UserData = [];         %unused
    newdst.GeoImage = dst;        %GeoImage is the dataset name for this format
end

%%
%--------------------------------------------------------------------------
% setDSproperties
%--------------------------------------------------------------------------
function dsp = setDSproperties()
    %define the variables and metadata properties for the dataset
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique
    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...                      
        'Name',{'geoimage'},...
        'Description',{'GeoImage'},...
        'Unit',{''},...
        'Label',{'GeoImage'},...
        'QCflag',repmat({'data'},1,1)); 
    dsp.Row = struct(...
        'Name',{'Location'},...
        'Description',{'Location'},...
        'Unit',{'-'},...
        'Label',{'Location'},...
        'Format',{''});         
    dsp.Dimensions = struct(...    
        'Name',{''},...
        'Description',{''},...
        'Unit',{''},...
        'Label',{''},...
        'Format',{''});   
end

%%
%--------------------------------------------------------------------------
% getPlot
%--------------------------------------------------------------------------
function ok = getPlot(obj,src,dsetname)
    %generate a plot on the src graphical object handle
    ok = [];
    tabcb  = @(src,evdat)tabPlot(obj,src);
    tabfigureplot(obj,src,tabcb,false);
    ax = gca;
    %get data and variable id
    dst = obj.Data.(dsetname);
    if isempty(dst), return; end
    %test for array of allowed data types for a color image
    img = dst.geoimage;     %image object
    him = imagesc(ax,'XData',img.XData,'YData',img.YData,'CData',img.CData);
    set(gca, 'YDir', 'normal'); % Correct the Y direction
    if any(isnan(img.CData))
        set(him, 'AlphaData', 1-isnan(img.CData)); %set Nan values to be transparent
    end
    colormap(img.CMap);
    clim(img.CLim);
    shading interp
    axis equal tight
    cb = colorbar;
    cb.Label.String = 'Elevation (mAD)';
    xlabel('Eastings (m)'); 
    ylabel('Northings (m)');    
    title(sprintf('%s(%s)',dst.Description,dsetname));

    % Draw a rectangle around the image
    hold on; % Keep the image displayed while adding the rectangle
    addBox(img)
    hold off;
end

%%
%--------------------------------------------------------------------------
% dataQC
%--------------------------------------------------------------------------
function output = dataQC(obj)                    %#ok<INUSD> 
    %quality control a dataset
    % datasetname = getDataSetName(obj); %prompts user to select dataset if more than one
    % dst = obj.Data.(datasetname);      %selected dstable
    warndlg('No quality control defined for this format');
    output = [];    %if no QC implemented in dataQC
end

%%
%--------------------------------------------------------------------------
% Utility functions
%--------------------------------------------------------------------------
function addBox(img)
    %add enclosing box to image
    mnmxX = minmax(img.XData);
    mnmxY = minmax(img.YData);
    rectangle('position',[img.XData(1),img.YData(1),diff(mnmxX),diff(mnmxY)],...
                                             'EdgeColor', 'k','linewidth',0.5);
    % Expand the axis limits to show full frame
    xlim(xlim()+[-.001,.001]*diff(mnmxX)) % add 1% to the x axis limits
    ylim(ylim()+[-.001,.001]*diff(mnmxY)) % add 1% to the y axis limits
end

 %%
 function im = readTiff(filename,ext)
    %read data from a tiff file
    cdata = flipud(imread(filename,ext));
    info = imfinfo(filename,ext);
    cmap = info.Colormap;
    climits = str2num(info.DocumentName); %#ok<ST2NM> vector
    nrows = info.Height;
    ncols = info.Width;
    dx = info.XResolution;              %dx
    dy = info.YResolution; 
    xLim = [info.XPosition,info.XPosition+ncols*dx];
    yLim = [info.YPosition,info.YPosition+nrows*dy];
    im = struct('XData',xLim,'YData',yLim,'CData',cdata,'CMap',cmap,'CLim',climits);
 end