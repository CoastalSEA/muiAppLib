function dst = tidalrange(wl,t,issave,isplot)
% 
%-------function help------------------------------------------------------
% NAME
%   tidalrange.m
% PURPOSE
%   function to compute tidal range from a water level time series
% USAGE
%   dst = tidalrange(wl,t,issave,isplot)
% INPUT
%   wl - water levels (mOD)
%   t - datetime
%   issave - true if range is to be returned and saved (optional)
%   isplot - true if plotting options are to be called (optional)
% OUTPUT (optional)
%   dst - dstable containing, tidal range (m),  high water levels (mAD), 
%         low water levels (mAD). Times are based on the central
%         up/down-crossing of mean water level. ie NOT explicit HW and LW
%         times. This provides the high and low water values associated
%         with each value of tidal range.
%
% Author: Ian Townend
% CoastalSEA (c)June 2018
%--------------------------------------------------------------------------
%
    if nargin<3           %catch call using tidalrange(wl,t)
        isplot = false;
        issave = false;
    elseif nargin<4      %catch call using tidalrange(wl,t,issave)
        isplot = false;         
    end
    %
    if isempty(issave)   %catch call using tidalrange(wl,t,[],isplot)
        issave = false;
    end
    meanwl = mean(wl,1,'omitnan');
    %get the up and down-crossing of mwl
    [upid,downid] = zero_crossing(wl,meanwl); %upid and downid are the same length
    %zero-crossing returns up-down or down-up pairs.The first up or down 
    %determine which.For a tide we need complete up-down-up of down-up-down cycles
    %ie up-to-up or down-to-down.

    umeanwl = mean(wl(upid(1):upid(end)),1,'omitnan');     %mean for exact number of cycles
    dmeanwl = mean(wl(downid(1):downid(end)),1,'omitnan'); %mean for exact number of cycles
    %check but for most tide records should be the same.
    meanwl = (umeanwl+dmeanwl)/2;

    %find the mean cycle interval
    umodedt = mode(diff(t(upid)));
    dmodedt = mode(diff(t(downid)));
    meandt  = hours((umodedt+dmodedt)/2);
    
    %find first crossing point
    isdowncross = downid(1)<upid(1);
    %find last crossing point
    isupcross = upid(end)>downid(end);
    % %find first tide with a max and a min
    % firstmax = max(wl(1:upid(1)),[],1,'omitnan');
    % %find last tide with a max and min
    % lastmin = min(wl(downid(end):end),[],1,'omitnan');
    ntide = length(upid)-1;
    if ~isdowncross && ~isupcross 
        %correct record length if starts on upcross and ends on downcross
        ntide = ntide-1;
    end
    %getRange assumes wls around 0 datum. If all wls +ve add offset
    %the 0.2 offset handles a few extremes going below zero
    if all(wl(~isnan(wl))>-0.2)
        offset = mean(wl,1,'omitnan');
        wl0 = wl-offset;
        offlag = true;
    else
        wl0 = wl;
        offlag = false;
    end
    
    %get the range after each up or down crossing    
    if isdowncross
        %starts on downcross
        [R,rt,hwl,lwl] = getRange(wl0,t,ntide,downid,upid,meanwl,meandt);
    else
        %starts on upcross
        [R,rt,hwl,lwl] = getRange(wl0,t,ntide,upid,downid,meanwl,meandt);   
    end
    
    %adjust the high and low water levels to remove offset
    if offlag
        hwl = hwl+offset;
        lwl = lwl+offset;
    end

    %report statistics of water levels and tidal range
    minwl = min(wl,[],1,'omitnan');
    maxwl = max(wl,[],1,'omitnan');
    minTR = min(R,[],1,'omitnan');
    maxTR = max(R,[],1,'omitnan');
    mnTR = mean(R,1,'omitnan');
    
    mhw = mean(hwl,1,'omitnan');
    mlw = mean(lwl,1,'omitnan');
    mhhw = mean(hwl(hwl>mhw),1,'omitnan');
    mlhw = mean(hwl(hwl<mhw),1,'omitnan');
    mhlw = mean(lwl(lwl>mlw),1,'omitnan');
    mllw = mean(lwl(lwl<mlw),1,'omitnan');
    
    varname = {'Elevation'};    
    rowname = {'Maximum water level';'Mean water level';'Minimum water level';...
               'Maximum tidal range';'Mean tidal range';'Minimum tidal range';...
               'MHHW';'MHW';'MLHW';'MHLW';'MLW';'MLLW'};
    wlstats = [maxwl;meanwl;minwl;maxTR;mnTR;minTR;mhhw;mhw;mlhw;mhlw;mlw;mllw];
    
    restable = table(wlstats,'RowNames',rowname,'VariableNames',varname);
    tablefigure('Tidal ranges','Summary tidal statistics',restable);
    
    if isplot
        %option to produce various plots    
        waterlevelfreqplots(wl,t);
    end
    %option to save output
    if issave
        dsp = setDSproperties();
        results = {R,hwl,lwl};
        dst = dstable(results{:},'RowNames',rt,'DSproperties',dsp);
        %Source and MetaData are set in muiUserModel. 
        %Put fit parameters in UserData
        dst.UserData = restable;
    else
        dst = 'no output';
    end
end
%%
function  [R,rt,hwl,lwl] = getRange(wl,t,ntide,stid,ndid,meanwl,meandt)
    %subsample waterlevels for each tidal cycle defined by stid(i) to (i+1)
    %asssign to mid-tide defined by ndid(i).
    %filters on the period between up/down crossings and number of null
    %records for a tidal cycle record to be void
    tfact = 1.2;  %allow some variance around mean modal value (meandt)
    nulls = 5;    %number of records that can be NaN in one half cycle
    R = zeros(ntide,1); hwl = R; lwl = R; rt = datetime;
    for i=1:ntide
        if hours(t(stid(i+1))-t(stid(i)))>tfact*meandt
            R(i,1) = NaN; hwl(i,1) = NaN; lwl(i,1) = NaN;
        else
            firstcyc = wl(stid(i):ndid(i));
            firstpk = max(abs(firstcyc),[],1,'omitnan');
            ishw = mean(firstcyc,'omitnan')>meanwl;
            secondcyc = wl(ndid(i):stid(i+1));
            secondpk = max(abs(secondcyc),[],1,'omitnan');
            if sum(isnan(firstcyc))>nulls || sum(isnan(secondcyc))>nulls 
                %gaps total O[1hr] for 10-15min data so peak may be missed
                R(i,1) = NaN; hwl(i,1) = NaN; lwl(i,1) = NaN;
            else
                R(i,1) = firstpk+secondpk;
                if ishw
                    hwl(i,1) = firstpk;
                    lwl(i,1) = -secondpk;
                else
                    hwl(i,1) = secondpk; 
                    lwl(i,1) = -firstpk;   
                end
            end
        end
        rt(i,1) = t(ndid(i));
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
        'Name',{'TR','HWL','LWL'},... 
        'Description',{'Tidal range','High water levels','Low water levels'},...
        'Unit',{'m','mAD','mAD'},...
        'Label',{'Tidal Range (m)','Water level (mAD)','Water level (mAD)'},...
        'QCflag',repmat({'derived'},1,3)); 
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
