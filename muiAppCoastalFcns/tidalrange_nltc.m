function dst = tidalrange_nltc(var,dtime,issave,isplot,titletxt)
%var,tdt,varargin
%-------function help------------------------------------------------------
% NAME
%   tidalrange_nltc.m
% PURPOSE
%   fit trend and cycles to tidal range or HW/LW time series
% USAGE
%   [ftime,fvals] = tidalrange_nltc(tr,tdt,varargin)
% INPUTS
%   var    - tidal range or water level
%   dtime  - time in datetime format
%   issave - true if range is to be returned and saved (optional)
%   isplot - true if plotting options are to be called (optional)
%   titletxt - text to add to title (e.g. Case used) - optional
% OUTPUT (optional)
%   dst - dstable containing fit to var (m), 
% SEE ALSO
%  tidalrange.m
%
% Author: Ian Townend
% CoastalSEA (c)Dec 2023
%--------------------------------------------------------------------------
%
    if nargin<3          %catch call using tidalrange_nltc(wl,t)
        isplot = false;
        issave = false;
        titletxt = [];
    elseif nargin<4      %catch call using tidalrange_nltc(wl,t,issave)
        isplot = false;  
        titletxt = [];
    elseif nargin<5
        titletxt = [];
    end
    y2s = 31556952;
    CyclePeriod = [18.61, 180.0, 8.85];
    CyclePhase = [12.2,128.0, 0];
    om1 = 2*pi./CyclePeriod/y2s;
    om2 = om1.*CyclePhase*y2s;

    dt_zero = datetime('01-Jan-0000 00:00:00');
    jtime = seconds(dtime-dt_zero);
    dt = mean(diff(jtime));    

    %compute regression fits for the various models
    [TC,amp,rsq] = getModelFits(jtime,var,om1,om2,dt); 
    %time in datetime to match TC array size
    intime = interp1(dtime,dtime,dtime(1):seconds(dt/5):dtime(end));

    %setup output table
    outable = setOutputTable(CyclePeriod,CyclePhase,om1,amp,rsq,y2s);
    headtxt = sprintf('Tidal signal regression for %s',titletxt);
    ht = tablefigure('Tidal cycles analysis',headtxt,outable);
    ht.Position = [200,750,980,270];
    
    %plot the data and regression lines
    if isplot
        plotFits(dtime,var,intime,TC,titletxt); 
        plotDiffs(dtime,var,intime,TC,titletxt); 
        ok = 1;
        while ok>0
            [idx,ok] = listdlg('PromptString','Select a Case',...
                           'SelectionMode','single','ListSize',[300,200],...
                           'ListString',outable.Properties.RowNames);
            if ok<1, continue; end
            plotSelectedFit(dtime,var,intime,TC,titletxt,idx,outable);
        end
    end

    %option to save output
    if issave
        [idx,ok] = listdlg('PromptString','Select a Case',...
                           'SelectionMode','single','ListSize',[300,200],...
                           'ListString',outable.Properties.RowNames);
        if ok<1, dst = 'no output'; return; end

        dsp = setDSproperties();
        dst = dstable(TC(:,idx),'RowNames',intime','DSproperties',dsp);
        %dst.Source = sprintf('Function tidalrange_nltc for %s',titletxt);
        %Source and MetaData are set in muiUserModel. 
        %Put fit parameters in UserData
        dst.UserData = outable(idx,:);
    else
        dst = 'no output';
    end
end

%%
function [TC,amp,rsq] = getModelFits(jtime,var,om1,om2,dt)
    %compute the rgressions for various models with trend and cycles
    ptime = interp1(jtime,jtime,jtime(1):dt/5:jtime(end));    

    %linear trend
    [a,b,rsq(1,1),~,TC(:,1),~] = regression_model(jtime,var,'linear',length(ptime)-1,false);
    %rsq(1) = r_squared(var,a+b*jtime);  %check that both give the same
    amp{1} = [a,b];

    %trend + nodal cycle
    funtr2 = @(amp,t) amp(1)+amp(2)*t+real(amp(3).*exp(1i*(om1(1)*t+om2(1))));

    amp0 = [a, b, 0.2];
    [amp{2},~,~,~,~,~] = nlinfit(jtime,var,funtr2,amp0);
    TC(:,2)= funtr2(amp{2},ptime);
    rsq(1,2) = r_squared(var,funtr2(amp{2},jtime));

    %trend + two cycles
    funtr3 = @(amp,t) amp(1)+amp(2)*t+real(amp(3).*exp(1i*(om1(1)*t + ...
                            om2(1)))+amp(4).*exp(1i*(om1(2)*t+om2(2))));

    amp0 = [a, b, 0.2, 0.2];
    [amp{3},~,~,~,~,~] = nlinfit(jtime,var,funtr3,amp0);
    TC(:,3) = funtr3(amp{3},ptime);
    rsq(1,3) = r_squared(var,funtr3(amp{3},jtime));

    %two cycles
    funtr4 = @(amp,t) amp(1)+real(amp(2).*exp(1i*(om1(1)*t + ...
                            om2(1)))+amp(3).*exp(1i*(om1(2)*t+om2(2))));
    
    amp0 = [a, 0.2, 0.2];
    [amp{4},~,~,~,~,~] = nlinfit(jtime,var,funtr4,amp0);
    TC(:,4) = funtr4(amp{4},ptime);   
    rsq(1,4) = r_squared(var,funtr4(amp{4},jtime));

    %two cycles and LTC phase as a fit parameter
    funtr5 = @(amp,t) amp(1)+real(amp(2).*exp(1i*(om1(1)*t + ...
                            om2(1)))+amp(3).*exp(1i*(om1(2)*t+amp(4))));

    amp0 = [a, 0.2, 0.2, om2(2)];
    [amp{5},~,~,~,~,~] = nlinfit(jtime,var,funtr5,amp0);
    TC(:,5) = funtr5(amp{5},ptime);   
    rsq(1,5) = r_squared(var,funtr5(amp{5},jtime));

   %two cycles and LTC period and phase as a fit parameter
    funtr6 = @(amp,t) amp(1)+real(amp(2).*exp(1i*(om1(1)*t + ...
                            om2(1)))+amp(3).*exp(1i*(amp(4)*t+amp(5))));

    amp0 = [a, 0.2, 0.2, om1(2), om2(2)];
    [amp{6},~,~,~,~,~] = nlinfit(jtime,var,funtr6,amp0);
    TC(:,6) = funtr6(amp{6},ptime);   
    rsq(1,6) = r_squared(var,funtr6(amp{6},jtime));

    %trend, two cycles and NTC/LTC phases as fit parameters
    funtr7 = @(amp,t) amp(1)+amp(2)*t+real(amp(3).*exp(1i*(om1(1)*t + ...
                            amp(4)))+amp(5).*exp(1i*(om1(2)*t+amp(6))));

    amp0 = [a, 0, 0.2, om2(1), 0.2, om2(2)];
    [amp{7},~,~,~,~,~] = nlinfit(jtime,var,funtr7,amp0);
    TC(:,7) = funtr7(amp{7},ptime);   
    rsq(1,7) = r_squared(var,funtr7(amp{7},jtime));

    %three cycles 8.85, 18.61 and 180 year with phases as fit parameters
    funtr8 = @(amp,t) amp(1)+real(amp(2).*exp(1i*(om1(3)*t+amp(3))) + ...
                                  amp(4).*exp(1i*(om1(1)*t+amp(5))) + ...
                                  amp(6).*exp(1i*(om1(2)*t+amp(7))));
    
    amp0 = [a, 0.1, om2(3), 0.2, om2(1), 0.2, om2(2)];
    [amp{8},~,~,~,~,~] = nlinfit(jtime,var,funtr8,amp0);
    TC(:,8) = funtr8(amp{8},ptime);   
    rsq(1,8) = r_squared(var,funtr8(amp{8},jtime));
    disp(amp{8})
    disp(rsq(8))
end

%%
function plotFits(dtime,var,intime,TC,titletxt)
    %plot the data and the resulting fit lines
    hf = figure('Tag','FigPlot'); 
    ax = axes(hf);    

    s1 = subplot(2,1,1,ax);
    plot(s1,dtime,var,'x','DisplayName','Data','ButtonDownFcn',@godisplay);
    hold on
    plot(s1,intime,TC(:,1),'DisplayName','Linear','ButtonDownFcn',@godisplay);    
    plot(s1,intime,TC(:,2),'DisplayName','Linear+NTC(a)','ButtonDownFcn',@godisplay);
    plot(s1,intime,TC(:,3),'DisplayName','Linear+NTC(a)+LTC(a)','ButtonDownFcn',@godisplay);
    plot(s1,intime,TC(:,7),'DisplayName','Linear+NTC(a,phi)+LTC(a,phi)','ButtonDownFcn',@godisplay);
    hold off
    xlabel('Time')
    ylabel('Variable')
    title('Linear')
    legend

    s2 = subplot(2,1,2);
    plot(s2,dtime,var,'x','DisplayName','Data','ButtonDownFcn',@godisplay);
    hold on
    plot(s2,intime,TC(:,4),'DisplayName','NTC(a)+LTC(a)','ButtonDownFcn',@godisplay);
    plot(s2,intime,TC(:,5),'DisplayName','NTC(a)+LTC(a,phi)','ButtonDownFcn',@godisplay);
    plot(s2,intime,TC(:,6),'DisplayName','NTC(a)+LTC(a,per,phi)','ButtonDownFcn',@godisplay);
    plot(s2,intime,TC(:,8),'DisplayName','LPC(a,phi)+NTC(a,phi)+LTC(a,phi)','ButtonDownFcn',@godisplay);
    hold off
    xlabel('Time')
    ylabel('Variable')
    title('Cycles')
    legend
    
    sgtitle(sprintf('Tidal signal regression for %s',titletxt))
end
%%
function plotDiffs(dtime,var,intime,TC,titletxt)
    %plot the differences between estimates and data
    hf = figure('Tag','FigPlot'); 
    ax = axes(hf);    
    
    diffs = getDifferences(dtime,var,intime,TC);

    s1 = subplot(2,1,1,ax);
    plot(s1,dtime,diffs(:,1),'DisplayName','Linear','ButtonDownFcn',@godisplay);
    hold on    
    %bar(s1,dtime,diffs(:,[2,3,7]))
    stem(s1,dtime,diffs(:,2),'MarkerSize',4,'DisplayName','Linear+NTC(a)','ButtonDownFcn',@godisplay);
    stem(s1,dtime,diffs(:,3),'MarkerSize',4,'DisplayName','Linear+NTC(a)+LTC(a)','ButtonDownFcn',@godisplay);
    stem(s1,dtime,diffs(:,7),'MarkerSize',4,'DisplayName','Linear+NTC(a,phi)+LTC(a,phi)','ButtonDownFcn',@godisplay);
    hold off
    xlabel('Time')
    ylabel('Variable')
    title('Linear')
    legend

    s2 = subplot(2,1,2);
    plot(s2,dtime,diffs(:,1),'DisplayName','Linear','ButtonDownFcn',@godisplay);
    hold on
    stem(s2,dtime,diffs(:,4),'MarkerSize',4,'DisplayName','NTC(a)+LTC(a)','ButtonDownFcn',@godisplay);
    stem(s2,dtime,diffs(:,5),'MarkerSize',4,'DisplayName','NTC(a)+LTC(a,phi)','ButtonDownFcn',@godisplay);
    stem(s2,dtime,diffs(:,6),'MarkerSize',4,'DisplayName','NTC(a)+LTC(a,per,phi)','ButtonDownFcn',@godisplay);
    stem(s2,dtime,diffs(:,8),'MarkerSize',4,'DisplayName','LPC(a,phi)+NTC(a,phi)+LTC(a,phi)','ButtonDownFcn',@godisplay);
    hold off
    xlabel('Time')
    ylabel('Variable')
    title('Cycles')
    legend

    sgtitle(sprintf('Tidal signal regression for %s',titletxt))
end
%%
function diffs = getDifferences(dtime,var,intime,TC)
    %interpolate the model onto the observed time step and compute differences
    diffs = zeros(length(var),size(TC,2));
    for i=1:size(TC,2)
        tempdiff = interp1(intime,TC(:,i),dtime);
        diffs(:,i) = var-tempdiff;
    end
end
%%
function plotSelectedFit(dtime,var,intime,TC,titletxt,idx,outable)
    %pot the seleted fit against the data with lines to show the longer
    %term cycle
    hf = figure('Tag','FigPlot'); 
    ax = axes(hf);  
    plot(ax,dtime,var,'x','DisplayName','Data','ButtonDownFcn',@godisplay);
    hold on
    legtxt = outable.Properties.RowNames{idx};
    plot(ax,intime,TC(:,idx),'DisplayName',legtxt,'ButtonDownFcn',@godisplay);
    if abs(outable{idx,'ampLTC(m)'})>0
        amp = outable{idx,1:end-1}; 
        amp(2) = amp(2)/100; amp(4) = 0; amp(5) = 0;        
        om1 = 2*pi./amp(10);
        om2 = om1.*amp(11);
        funtr = @(amp,t) amp(1) +amp(2)*t +real(amp(3)+amp(6)+amp(9).*exp(1i*(om1*t+om2)));
        dt_zero = datetime('01-Jan-0000 00:00:00');
        ptime = years(intime-dt_zero);
        FC = funtr(amp,ptime);   
        plot(ax,intime,FC,'--g','DisplayName','LTC fit','ButtonDownFcn',@godisplay);
        amp(1) = amp(1)-2*(amp(3)+amp(6));
        FC = funtr(amp,ptime);  
        h1 = plot(ax,intime,FC,'--g','DisplayName','LTC fit','ButtonDownFcn',@godisplay);
        h1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    end
    hold off
    xlabel('Time')
    ylabel('Tidal signal')
    legend    
    title(sprintf('Tidal signal regression for %s',titletxt))
end
%%
function outable = setOutputTable(CyclePeriod,CyclePhase,om1,amp,rsq,y2s)
    %format the output table
    perNTC = CyclePeriod(1); perLTC = CyclePeriod(2); perLPC = CyclePeriod(3);
    phaNTC = CyclePhase(1);  phaLTC = CyclePhase(2);
    pha_LPC = @(amp) round(amp/om1(3)/y2s);  %years
    pha_NTC = @(amp) round(amp/om1(1)/y2s);  %years
    pha_LTC = @(amp) round(amp/om1(2)/y2s);  %years
    per_LTC = @(amp) round(2*pi./amp/y2s);   %years
    slp = @(amp) amp*y2s*100;                %m/100 years
    
    outable = table(amp{1}(1),slp(amp{1}(2)),0,0,0,0,0,0,0,0,0,rsq(1),...
        'VariableNames',{'intercept(m)','slope(m/100y)','ampLPC(m)','perLPC(y)','phaLPC(y)',...
        'ampNTC(m)','perNTC(y)','phaNTC(y)','ampLTC(m)','perLTC(y)','phaLTC(y)','R2'});
    
    addres = {amp{2}(1),slp(amp{2}(2)),0,0,0,amp{2}(3),perNTC,phaNTC,0,0,0,rsq(2);...
        amp{3}(1),slp(amp{3}(2)),0,0,0,amp{3}(3),perNTC,phaNTC,amp{3}(4),perLTC,phaLTC,rsq(3);...
        amp{4}(1),0,0,0,0,amp{4}(2),perNTC,phaNTC,amp{4}(3),perLTC,phaLTC,rsq(4);...
        amp{5}(1),0,0,0,0,amp{5}(2),perNTC,phaNTC,amp{5}(3),perLTC,pha_LTC(amp{5}(4)),rsq(5);...
        amp{6}(1),0,0,0,0,amp{6}(2),perNTC,phaNTC,amp{6}(3),per_LTC(amp{6}(4)),pha_LTC(amp{6}(4)),rsq(6);...
        amp{7}(1),slp(amp{7}(2)),0,0,0,amp{7}(3),perNTC,pha_NTC(amp{7}(4)),amp{7}(5),perLTC,pha_LTC(amp{7}(6)),rsq(7);...
        amp{8}(1),0,amp{8}(2),perLPC,pha_LPC(amp{8}(3)),amp{8}(4),perNTC,pha_NTC(amp{8}(5)),amp{8}(6),perLTC,pha_LTC(amp{8}(7)),rsq(8)};
    
    outable = [outable;addres];
    
    outable.Properties.RowNames = {'Linear','Linear+NTC(a)','Linear+NTC(a)+LTC(a)',...
        'NTC(a)+LTC(a)','NTC(a)+LTC(a,phi)',...
        'NTC(a)+LTC(a,per,phi)','Linear+NTC(a,phi)+LTC(a,phi)',...
        'LPC(a,phi)+NTC(a,phi)+LTC(a,phi)'};
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
        'Name',{'TRfit'},... 
        'Description',{'Tidal range'},...
        'Unit',{'m'},...
        'Label',{'Tidal Range (m)'},...
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