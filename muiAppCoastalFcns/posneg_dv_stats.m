function ok = posneg_dv_stats(v,t,varname,method)
%
%-------function help------------------------------------------------------
% NAME
%   posneg_dv_stats.m
% PURPOSE
%   Compute the rate of change of variable and plot histograms for positive
%   and negative components
% USAGE
%   ok = posneg_dv_stats(v,t,varname,method)
% INPUTS: 
%   v - variable
%   t - time
%   varname - a text description of the variable (optional)
%   method - defines the case: only 1 implemented = posnegVals (optional)
% OUTPUTS 
%   timeseries plot of change histogram plot of Poisson fits to data
% REQUIREMENTS
%   This function posneg_dv_stats uses fitdist and pdf from the
%   Statistics and Machine Learning toolbox
%
% Author: Ian Townend
% CoastalSEA (c) April 2017
%----------------------------------------------------------------------
%

    %dv and dt are padded so that they are the same length as the source
    %time series
    if nargin<3
        varname = 'Data';
        method = 1;
    elseif nargin<4 && ~isempty(varname)
        method = 1;        
    elseif isempty(varname)
        varname = 'Data';
    end

    dV = [NaN;diff(v)];              %volume change
    dt = [NaN;days(diff(t))];        %time interval
    dVdt = dV./dt;                   %rate of change of variable
    plot_dVdt(v,t,dV,dVdt,varname)   %plot input data and rates of change
    
    switch method
        case 1
            pdV = posnegVals(dV,dVdt,varname);
        case 2
    end
    %return message with results
    ok = sprintf('Exponential mu:\npostive dV %g\npostive dVdt %g\nnegative dV %g\nnegative dVdt %g',...
           pdV.mu_pdV.mu,pdV.mu_pdVdt.mu,pdV.mu_ndV.mu,pdV.mu_ndVdt.mu);
end
%%
function pnVr = posnegVals(dV,dVdt,varname)
    % -------------------------------------------------------------------------
    % Method 1
    % seperate positive and negative values and find pdf fit
    % this assumes that all profiles are sufficiently similar that the data can
    % be combined.
    % -------------------------------------------------------------------------
    pad = num2cell(zeros(1,12));
    pnV = table(pad{:},'VariableNames',{'pv','pdv','nv','ndv',...
                                         'pdV','pd_pdV',...
                                         'pdVdt','pd_pdVdt',...
                                         'ndV','pd_ndV',...
                                         'ndVdt','pd_ndVdt'});                                 
    nint = 20;
    pnV.pv = 0:max(dV)/nint:max(dV);
    pnV.pdv = 0:max(dVdt)/nint:max(dVdt);
    pnV.nv = -0:min(dV)/nint:min(dV);
    pnV.ndv = 0:min(dVdt)/nint:min(dVdt);   

    pdf_name = 'Exponential';
    pnV.pdV = dV(dV>=0)';
    pnVr.mu_pdV = fitdist(pnV.pdV',pdf_name);
    pnV.pd_pdV = pdf(pnVr.mu_pdV,pnV.pv);
    %
    pnV.pdVdt = dVdt(dVdt>=0)';
    pnVr.mu_pdVdt = fitdist(pnV.pdVdt',pdf_name);
    pnV.pd_pdVdt = pdf(pnVr.mu_pdVdt,pnV.pdv);
    %
    pnV.ndV = dV(dV<0)';
    pnVr.mu_ndV = fitdist(-pnV.ndV',pdf_name);
    pnV.pd_ndV = pdf(pnVr.mu_ndV,-pnV.nv);
    %
    pnV.ndVdt = dVdt(dVdt<0)';
    pnVr.mu_ndVdt = fitdist(-pnV.ndVdt',pdf_name);
    pnV.pd_ndVdt = pdf(pnVr.mu_ndVdt,-pnV.ndv);
    %plot results
    plotdVarStats(pnV,pnVr,varname);
end
%%
function plotdVarStats(pnV,pnVr,varname)
    % plot results as positive and negative pdfs
    %
    if iscell(varname), sgtxt = varname{2}; varname = varname{1}; end

    vartxt = sprintf('%s change',varname);
    hf = figure('Name','Rates of change','Units','normalized',...
                'Tag','PlotFig');
    nbins = 10;
    norm_type = 'pdf';
    %
    subplot(2,2,1)
    histogram(pnV.ndV,nbins,'Normalization',norm_type);
    title(sprintf('Negative values: mu=%0.2g',pnVr.mu_ndV.mu))
    xlabel(vartxt)
    ylabel('Probability')
    hold on
    plot(pnV.nv,pnV.pd_ndV,'--r')
    xlim([min(pnV.nv),0]);
    hold off
    %
    subplot(2,2,2)
    histogram(pnV.pdV,nbins,'Normalization',norm_type);
    title(sprintf('Positive values: mu=%0.2g',pnVr.mu_pdV.mu))
    xlabel(vartxt)
    hold on
    plot(pnV.pv,pnV.pd_pdV,'--r')
    xlim([0,max(pnV.pv)]);
    hold off
    %
    subplot(2,2,3)
    histogram(pnV.ndVdt,nbins,'Normalization',norm_type);
    title(sprintf('Negative rates: mu=%3.3g',pnVr.mu_ndVdt.mu))
    xlabel('Rate of change (per day)')
    ylabel('Probability')
    hold on
    plot(pnV.ndv,pnV.pd_ndVdt,'--r')
    xlim([min(pnV.ndv),0]);
    hold off
    %
    subplot(2,2,4)
    histogram(pnV.pdVdt,nbins,'Normalization',norm_type);
    title(sprintf('Positive rates: mu=%0.2g',pnVr.mu_pdVdt.mu))
    xlabel('Rate of change (per day)')
    hold on
    plot(pnV.pdv,pnV.pd_pdVdt,'--r')
    xlim([0,max(pnV.pdv),]);
    hold off
    
    sgtitle(sgtxt)
    out_mu = table(pnVr.mu_ndV.mu,pnVr.mu_pdV.mu,pnVr.mu_ndVdt.mu,pnVr.mu_pdVdt.mu,...
                   'VariableNames',{'mu_ndV','mu_pdV','mu_ndVdt','mu_pdVdt'});
    outable = [pnV,out_mu];
    add_copy_button(hf,outable);
end
%%
function plot_dVdt(v,t,dV,dVdt,varname)
    %generate plot volume timeseries and rates of change
    if iscell(varname), sgtxt = varname{2}; varname = varname{1}; end
    
    dt = [NaN;days(diff(t))];  %time interval between surveys in days    
    fprintf('Mean time interval = %0.3f days',mean(dt,'omitnan'))
    
    hf = figure('Name','Magnitude of change','Units','normalized',...
                'Tag','PlotFig');
    
    subplot(3,1,1)
    plot(t,v);
    ylabel(varname)
    title(sprintf('%s change',varname));

    s2 = subplot(3,1,2);
    coloroptions = s2.ColorOrder;      %yyaxis forces single ColorOrder option
    s2.ColorOrder = [0,0,0];           %set axis label color
    yyaxis(s2,'left')
    s2.ColorOrder = coloroptions(2,:);
    bar(t,dV,ceil(length(dV)/50));
    % plot(t,dV)
    ylabel('Change (dV)')
    %
    yyaxis(s2,'right') 
    s2.ColorOrder = coloroptions(4,:);
    stem(t,dt);
    % plot(t,dt);
    ylabel('Survey interval (d)');
    
    s3 = subplot(3,1,3);
    bar(t,dVdt,ceil(length(dV)/50),'FaceColor',s3.ColorOrder(3,:));
    % plot(t,dVdt,'Color',s3.ColorOrder(3,:));
    xlabel('Year')
    ylabel('Rate of change (per day)')
    
    sgtitle(sgtxt)
    outable = dstable(v,dV,dVdt,'RowNames',t,'VariableNames',{'V','dV','dVdt'});
    add_copy_button(hf,outable);
end




