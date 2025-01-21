function [hypsdst,histdst]  = gd_basin_hypsometry(grid,wl,histint,limits,isHW)
%
%-------function help------------------------------------------------------
% NAME
%   gd_basin_hypsometry.m
% PURPOSE
%   compute area and volume hypsometry from gridded elevation data
% USAGE
%   [hypsdst,histdst]  = gd_basin_hypsometry(grid,wl,histint,limits,isHW)
% INPUTS
%   grid - struct of x,y,z,t values that define grid 
%   wl - struct of zhw, zmt and zlw with scalar values at the mouth, or a
%        vector with the same dimension as the x vector, or CF_HydroData instance
%   histint - vertical interval for histogram
%   limits - upper and lower limits for vertical range:
%            0=use grid; 1=prompt user; [x,y]=limits to use
%   isHW - flag to apply mask to data. true = values above HW removed (optional)
%          allows limits to be set to range of multiple grids with varying water
%          levels (eg due to slr) but masks individual grids to specific HW
% OUTPUT
%   hypsdst - zcentre, zsurf, zvol, zhist as a dstable with dimensions of z
%            z: vertical elevation of interval mid-point
%            zsurf: surface plan area for each interval
%            zvol: cumulative volumes below each interval
%            zhist: plan area for each interval (zsurf is cumulative total)
%   histdst - elevation histogram, SArea, as a dstable with dimension of x and z
% NOTES
%   This function is only valid for cartesian grids of constant grid size
% SEE ALSO
%   used in GDinterface, along with gd_gross_properties and
%   gd_section_properties
%
% Author: Ian Townend
% CoastalSEA (c) July 2022
%--------------------------------------------------------------------------
%
    if nargin<5
        isHW = false;
    end
    
    ich = gd_basin_indices(grid);  %account for offset to mouth and head (if defined)
    % fprintf('xM = %.1f ich = %d-%d\n',grid.xM,ich(1),ich(end)) %debug check
    if isscalar(wl.zhw)
        zhw = ones(size(grid.x,1))*wl.zhw;
    else
        zhw = wl.zhw;
    end
    %grid intervals
    gdim = gd_dimensions(grid);

    % range for histogram data - general
    uplim = max(zhw(ich))+histint;
    lowlim = min(min(grid.z))-histint;    
    if (isscalar(limits) && limits==1) || length(limits)==2
        %use limits provided or prompt user to define limits
        [uplim,lowlim] = setLimits(limits,uplim,lowlim);
    end
        
    zedge = lowlim:histint:uplim;
    %centre of each bin
    zcentre = movsum(zedge,[0,1])/2;
    zcentre(end) = [];
    
    zhist = zeros(gdim.xint,length(zcentre)); %NB uses full grid

    for ix=ich
        %to account for varying upper and lower surfaces compute surface 
        %and volume for each increment along length        
        zed = grid.z(ix,:);        
        ibot = find(zedge>min(zed)-histint,1,'first');         
        itop = find(zedge<zhw(ix)+histint,1,'last');
        if ibot<itop
            subzedge = zedge(ibot:itop); %sub-sample that covers vertical range at z(ix)
            zhistc = histcounts(zed,subzedge); %bin counts for each interval defined by zedge
            zhist(ix,ibot:itop-1) = zhistc*gdim.delx*gdim.dely; %scale occurences by grid cell area
            if isHW
                %apply mask to zedge above zhw
                zhist(ix,subzedge>zhw(ix)+histint) = 0;
            end
        else
            zhist(ix,ibot:itop-1) = 0;
        end
    end
    
    %cumulative values from head along x-axis
    zhistsumx = sum(zhist,1);    %sum(bin counts x cell area) over all-x to give total area at each z
    zsurf = cumsum(zhistsumx);   %cumulative surface area at or below each elevation z
    zvol = cumtrapz(histint,zsurf); %cumulative volume below each elevation z
    hyps  = [{zsurf},{zvol},{zhistsumx}]; %formated for dstable input   

    %write cumulative basin hypsometry to a dstable 
    hdsp = setDSproperties_hypsometry;
    hypsdst = dstable(hyps{:},'RowNames',grid.t,'DSproperties',hdsp);    
    hypsdst.Dimensions.Z = zcentre; %dstable holds as a column vector
    hypsdst.UserData.histint = histint;      %save vertical interval used
    hypsdst.UserData.limits = [uplim,lowlim];%save z-limits used
    
    %write x-z histogram of hypsometry surfacae areas to a dstable
    sz = num2cell(size(zhist));
    zhist = reshape(zhist,1,sz{:});
    hgdsp = setDSproperties_histogram();
    histdst = dstable(zhist,'RowNames',grid.t,'DSproperties',hgdsp);    
    histdst.Dimensions.X = grid.x;  %NB uses full grid
    histdst.Dimensions.Z = zcentre; %dstable holds as a column vector
    %assign metadata about grid
    if ~isfield(grid,'desc'), grid.desc = 'New grid'; end
    histdst.Source = [grid.desc,': using gd_basin_hypsometry'];
    histdst.UserData.histint = histint;      %save vertical interval used
    histdst.UserData.limits = [uplim,lowlim];%save z-limits used
    histdst.Description = grid.desc;
end
%%
function [uplim,lowlim] = setLimits(limits,uplim,lowlim)
    %set the limits for the vertical range
    if length(limits)==2
        uplim = limits(1);
        lowlim = limits(2);
    else
        promptxt = {'Upper limit','Lower limit'};
        defaults = {num2str(uplim),num2str(lowlim)};
        title = 'Set limits';
        inp = inputdlg(promptxt,title,1,defaults);
        if isempty(inp), return; end %user cancels - use existing values
        uplim = str2double(inp{1});
        lowlim = str2double(inp{2});
    end
end
%%
function dsp = setDSproperties_histogram()
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
                    'Name',{'SArea'},...                  
                    'Description',{'Surface area'},...
                    'Unit',{'m^2'},...
                    'Label',{'Surface area (m^2)'},...
                    'QCflag',{'statistic'});  
    dsp.Row = struct(...
                    'Name',{'Time'},...
                    'Description',{'Time'},...
                    'Unit',{'y'},...
                    'Label',{'Time (yr)'},...
                    'Format',{'y'});          
    dsp.Dimensions = struct(...    
                    'Name',{'X','Z'},...
                    'Description',{'Distance','Elevation'},...
                    'Unit',{'m','mAD'},...
                    'Label',{'Distance (m)','Elevation (mAD)'},...
                    'Format',{'-','-'}); 
end
%%
function dsp = setDSproperties_hypsometry()
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
                    'Name',{'SurfaceArea','Volume','SAfreq'},...                  
                    'Description',{'Surface Area','Volume',...
                            'Surface area frequency'},...
                    'Unit',{'m2','m3','-'},...
                    'Label',{'Surface area (m^2)','Volume (m^3)',...
                             'Surface area frequency'},...
                    'QCflag',{'model','model','statistic'});  
    dsp.Row = struct(...
                    'Name',{'Time'},...
                    'Description',{'Time'},...
                    'Unit',{'y'},...
                    'Label',{'Time (yr)'},...
                    'Format',{'y'});          
    dsp.Dimensions = struct(...    
                    'Name',{'Z'},...
                    'Description',{'Elevation'},...
                    'Unit',{'mAD'},...
                    'Label',{'Elevation (mAD)'},...
                    'Format',{'-'});  
end
