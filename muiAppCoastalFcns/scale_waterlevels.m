function dst = scale_waterlevels(wl,t,issave,isplot)
% 
%-------function help------------------------------------------------------
% NAME
%   scale_waterlevels.m
% PURPOSE
%   function to scale the water levels based on factors for above and
%   below zero
% USAGE
%   dst = scale_waterlevels(wl,t,issave,isplot)
%   e.g. in Derive output UI: scale_waterlevels(x,t,1,0)
% INPUT
%   wl - water levels (mOD)
%   t - datetime
%   issave - true if range is to be returned and saved (optional)
%   isplot - true if plotting options are to be called (optional)
% OUTPUT (optional)
%   dst - dstable containing
%
% Author: Ian Townend
% CoastalSEA (c)Dec 2024
%--------------------------------------------------------------------------
%
    if nargin<3           %catch call using scale_waterlevels(wl,t)
        isplot = false;
        issave = false;
    elseif nargin<4      %catch call using scale_waterlevels(wl,t,issave)
        isplot = false;         
    end
    %
    if isempty(issave)   %catch call using scale_waterlevels(wl,t,[],isplot)
        issave = false;
    end
    promptxt = {'Scale for high waters', 'Scale for low waters'};
    answer = inputdlg(promptxt,'Scale WLs',1,{'1','1'});
    if isempty(answer), return; end  %user cancelled
    upper = str2double(answer{1});
    lower = str2double(answer{2}); 

    wl(wl>0) = wl(wl>0)*upper;
    wl(wl<=0) = wl(wl<=0)*lower;

    if isplot
        %option to produce various plots    
        figure; plot(t,wl); ylabel('Scaled water levels');
    end
    %option to save output
    if issave
        dsp = setDSproperties();
        dst = dstable(wl,'RowNames',t,'DSproperties',dsp);
        %Source and MetaData are set in muiUserModel. 
        %Put fit parameters in UserData
        dst.UserData = sprintf('[wl>0]*%.1f, [wl<=0]*%.1f',upper,lower);
    else
        dst = 'no output';
    end
end

%%
function dsp = setDSproperties()
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'WL'},... 
        'Description',{'Water levels'},...
        'Unit',{'mAD'},...
        'Label',{'Water level (mAD)'},...
        'QCflag',repmat({'derived'},1,1)); 
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