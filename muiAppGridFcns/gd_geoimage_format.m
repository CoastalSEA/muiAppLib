function output = gd_geoimage_format(funcall,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   gd_geoimage_format.m
% PURPOSE
%   Functions to define metadata, read and load data from file for estuary
%   image data
% USAGE
%   output = gd_geoimage_format(funcall,varargin)
% INPUTS
%   funcall - function being called
%   varargin - function specific input (filename,class instance,dsp,src, etc)
% OUTPUT
%   output - function specific output
% NOTES
%   This file loads a geoimage created from a grid
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
    obj.DataFormats = {'muiUserData',formatfile,'geoimage'};
    obj.idFormat = 1;
    obj.FileSpec = {'on','*.mat;'};
end
%%
%--------------------------------------------------------------------------
% getData
%--------------------------------------------------------------------------
function newdst = getData(obj,filename,metatxt) %#ok<INUSD>
    %read and load a data set from a file
    dsp = setDSproperties;                 %set metadata
    [~,location,~] = fileparts(filename);   
    imdata = {load(filename)};
    %load the results into a dstable - Image is the dataset name for this format
    dst = dstable(imdata,'RowNames',{location},'DSproperties',dsp);  
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
    im = dst.geoimage;     %image object
%     h_im = imagesc(ax,'XData',im.XData,'YData',im.YData,'CData',im.CData);
%     set(h_im, 'AlphaData', 1-isnan(im.CData)); %set Nan values to be transparent
h_im = imshow(im.CData, 'XData',im.XData,'YData',im.YData); 
    axis equal tight
%     cb = colorbar;
%     cb.Label.String = 'Elevation (mAD)';
    xlabel('Eastings (m)'); 
    ylabel('Northings (m)');    
    title(dst.Description);
    ax.Color = [0.96,0.96,0.96];  %needs to be set after plot
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