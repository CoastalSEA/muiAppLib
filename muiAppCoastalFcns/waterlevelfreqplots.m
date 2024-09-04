function res = waterlevelfreqplots(wl,t)
%
%-------function help------------------------------------------------------
% NAME
%   waterlevelfreqplots.m
% PURPOSE
%   Create various water level exceedance and duration plots
% USAGE
%   res = waterlevelfreqplots(wl,t)
% INPUTS
%   wl - waterlevel data
%   t - time
% OUTPUT
%   res - dummy text so that function can be called from Derive Output UI
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%
    res = {'Plots completed'}; %cell ouput required by call from DataManip.createVar   
    ok = 1;
    while ok>0
        %allow user to generate various plots
        plotlist = {'Water level elevation frequency',...
                    'Water level spectrum',...
                    'Elevations above a threshold',...
                    'Duration of threshold exceedance',...
                    'Elevation frequency above threshold',...
                    'Rolling mean duration above a threshold'};
        [idx,ok] = listdlg('Name','Plot options', ...
            'PromptString','Select a plot:', ...
            'SelectionMode','single', ...
            'ListSize',[250,150],...
            'ListString',plotlist);
        if ok<1, return; end

        switch plotlist{idx}
            case 'Water level elevation frequency'
                wl_elev_freq(wl)
            case 'Water level spectrum'
                wl_spectrum(wl,t)       
            otherwise            
                %get threshold elevation from user
                prompt = {'Threshold elevation (mOD):'};
                dlgtitle = 'Define elevation (mOD)';
                numlines = 1;
                defaultvalues{1} = num2str(0);
                useInp=inputdlg(prompt,dlgtitle,numlines,defaultvalues);
                if isempty(useInp), return; end %user cancelled
                z0 = str2double(useInp{1});

                switch plotlist{idx}
                    case 'Elevations above a threshold'
                        wl_exceed_thr(wl,t,z0)
                    case 'Duration of threshold exceedance'   
                        thr_durations(wl,t,z0)
                    case 'Elevation frequency above threshold'
                        thr_elev_freq(wl,z0)
                    case 'Rolling mean duration above a threshold'
                        movingtime_thr_duration(wl,t,z0)
                end
        end
    end
end
%%
function wl_elev_freq(wl)
    %Plot - 'Water level elevation frequency'
    %plot histogram of probability of a given water level against elevation
    minwl = min(wl,[],1,'omitnan');
    maxwl = max(wl,[],1,'omitnan');
    zedges = floor(minwl):0.1:ceil(maxwl);
    fc = histcounts(wl,zedges,'Normalization', 'probability');
    z = zedges(1)+0.05:0.1:zedges(end)-0.05;
    figure('Name','Elevation frequency','Units','normalized',...                
           'Resize','on','HandleVisibility','on','Tag','PlotFig');            
    barh(z,fc*100);
    title('Elevation frequency')
    ylabel('Elevation (mOD)')
    xlabel('Probability of occurrence (%)')
end
%%
function wl_spectrum(wl,t)
    %Plot - 'Water level spectrum'
    %code to compute the fft of the water level signal (based on Matlab Help)
    wl(isnan(wl)) = [];  %remove Nans
    Y= fft(wl);
    Ts = seconds(mean(diff(datenum(t))));
    Fs = seconds(1)/Ts;
    L = length(wl);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = 1/Fs/24.*(0:(L/2))/L;
    figure('Name','Water level spectrum','Units','normalized',...                
           'Resize','on','HandleVisibility','on','Tag','PlotFig'); 
    plot(f,P1)
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')     
end
%%
function wl_exceed_thr(wl,t,z0)
    %Plot - 'Elevations above a threshold'
    %plot probability of exceeding the specified threshold   
    wl(wl<z0) = NaN;
    figure('Name','Elevation exceedance','Units','normalized',...                
           'Resize','on','HandleVisibility','on','Tag','PlotFig'); 
    plot(t,wl);
    title(sprintf('Elevations above %.3g (mOD) threshold',z0))
    ylabel('Elevation (mOD)')
    xlabel('Time')
end
%%
function thr_durations(wl,t,z0)
    %Plot - 'Duration of threshold exceedance' 
    [stid,edid] = zero_crossing(wl,z0);
    if isempty(stid)
        hw = warndlg('No exceedances found'); 
        waitfor(hw);
        return; 
    end
    %
    if stid(1)>edid(1)
        temp = stid;
        stid = edid;
        edid  = temp;
    end
    wetdur = t(edid)-t(stid);
    wetdur.Format = 'h';
    figure('Name','Duration exceedance','Units','normalized',...                
           'Resize','on','HandleVisibility','on','Tag','PlotFig'); 
    sledges = min(wetdur):(max(wetdur)-min(wetdur))/10:max(wetdur);
    histogram(wetdur,sledges,'Normalization', 'probability');
    title(sprintf('Duration of %.3g (mOD) threshold exceedance',z0))
    ylabel('Probability of occurrence (%)')
    xlabel('Duration (hours)')
end
%%
function thr_elev_freq(wl,z0)
    %Plot - 'Elevation frequency above threshold'
    wl(wl<z0) = NaN;
    wledges = floor(min(wl)):0.1:ceil(max(wl));
    fc = histcounts(wl,wledges,'Normalization', 'probability');
    wl = wledges(1)+0.05:0.1:wledges(end)-0.05;
    figure('Name','Elevation frequency','Units','normalized',...                
           'Resize','on','HandleVisibility','on','Tag','PlotFig'); 
    barh(wl,fc*100);
    title(sprintf('Elevation frequency above %.3g (mOD) threshold',z0))
    ylabel('Elevation (mOD)')
    xlabel('Probability of occurrence (%)')
end
%%
function movingtime_thr_duration(wl,t,z0)
    %Plot the rolling mean duration above the threshold
    [stid,edid] = zero_crossing(wl,z0);
    if isempty(stid)
        hw = warndlg('No exceedances found'); 
        waitfor(hw);
        return; 
    end
    %
    if stid(1)>edid(1)
        temp = stid;
        stid = edid;
        edid  = temp;
    end
    wetdur = t(edid)-t(stid);
    wetdur.Format = 'h';
    wetdur = hours(wetdur);
    tper = years(1);  %set to annual but movingtime allows this to be changed!
    [tm,vm] = movingtime(wetdur,t(stid),tper,tper,'mean');
    figure('Name','Duration exceedance','Units','normalized',...                
           'Resize','on','HandleVisibility','on','Tag','PlotFig'); 
    plot(tm,vm);
    title(sprintf('Rolling mean above %.3g (mOD) threshold',z0))
    ylabel('Mean annual duration of exceedances (hours)')
    xlabel('Time')
    reclen = t(end)-t(1);
    pcntexcdur = sum(wetdur)/hours(reclen)*100;
    numexc = mean(length(stid),length(edid));
    aveannumexc = numexc/years(reclen);
    msg1 = sprintf('Percentage time duration exceeded in %.3g years = %.3g%%',years(reclen),pcntexcdur);
    msg2 = sprintf('Average annual number of threshold exceedances = %.3g',aveannumexc);
    msgtxt = sprintf('%s\n%s',msg1,msg2);
    hm = msgbox(msgtxt,'Mean duration results');
    waitfor(hm)
end
