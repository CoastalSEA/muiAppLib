function dst = tidalhighlow(wl,t)    
% 
%-------function help------------------------------------------------------
% NAME
%   tidalhighlow.m
% PURPOSE
%   function to compute tidal high and low waters from water level timeseries
% USAGE
%   dst = tidalhighlow(wl,t)   
% INPUT
%   wl - water levels (mOD)
%   t - datetime
% OUTPUT (optional)
%   dst - struct (fields HW, LW) with two dstables containing high water  
%         levels (mAD) and low water levels (mAD).
% NOTES
%   function can be used in Derive Output
% SEE ALSO
%   tidalrange.m compute tidal range and associated high and low waters
%   with times based on the cycle mid point.
%
% Author: Ian Townend
% CoastalSEA (c)Feb 2026
%--------------------------------------------------------------------------
%    
    meanwl = mean(wl,1,'omitnan');
    %get the up and down-crossing of mwl 
    
    %use peakoverthrwshold based on mean wl to find high and low water
    % method: 3 = peaks that have a separation of at least tint hours
    % tint: time interval between peaks; separation is >=tint (hours) 
    tint = hours(12);
    [locs,hwl] = peaksoverthreshold(wl,meanwl,3,t,tint);
    tHW = t(locs);
    
    [locs,lwl] = peaksoverthreshold(-wl,-meanwl,3,t,tint);
    tLW = t(locs);
    
    %check plot
    % figure;plot(t,wl)
    % hold on
    % plot(tHW,hwl,'xr')
    % plot(tLW,-lwl,'or')
    % hold off

    %save output
    dsp = setDSproperties('HW');
    dst.HW = dstable(hwl(:),'RowNames',tHW,'DSproperties',dsp);
    dsp = setDSproperties('LW');
    dst.LW = dstable(-lwl(:),'RowNames',tLW,'DSproperties',dsp);
    %Source and MetaData are set in muiUserModel.
end

%%

%%
function dsp = setDSproperties(option)
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    if strcmp(option,'HW')
        dsp.Variables = struct(...
            'Name',{'HWL'},... 
            'Description',{'High water levels'},...
            'Unit',{'mAD'},...
            'Label',{'Water level (mAD)'},...
            'QCflag',{'derived'});     
    else
        dsp.Variables = struct(...
            'Name',{'LWL'},... 
            'Description',{'Low water levels'},...
            'Unit',{'mAD'},...
            'Label',{'Water level (mAD)'},...
            'QCflag',{'derived'});         
    end

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