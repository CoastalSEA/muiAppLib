function dst = littoraldriftstats(qs,tdt,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   littoraldriftstats.m
% PURPOSE
%   Estimate net drift over selected period from the time series of drift rates
% USAGE
%   dst = littoraldriftstats(qs,tdt,'month');
% INPUTS
%   qs     - alongshore drift rate (m3/s) as a time series
%   tdt    - time as a datetime array
%   varargin - period: optional variable to define annual or monthly output 
%              to be returned in var. Use 'year' or 'month'. Default is
%              'month'.
%              isprompt: logical optional variable set to false to suppress 
%              prompts for ouput. Default is true. NB must be logical input
%              OR
%              plotoption: numerical value - 0=default plot, 1-3 plot only
%              subplot as numbered
% OUTPUT
%   Plot of gaps and monthly/annual drift plot
%   Defintions: gap in record = difference between dt*ndt in year and nrecs in each year
%               gap in data = data values that are NaN
%               calms = |qs|<calmsthreshold
%   Option to return time series of results (e.g. monthly values) 
%   dst - dstable if defined period is chosen or a struct of dstables if 
%         annual and monthly are to be saved. Each table contains:
%         qtime, qdrift - depends on selection 
%   Otherwise returns the average annual values (return value required by muiUserModel) 
% NOTES
%	designed to be called from Derive Output in any muitoolbox model UI
% SEE ALSO
%   littoraldrift.m, xshore_bailard.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    maxgap = 12;  %to avoid drift being applied over large data gaps set 
                  %maximum gap that can be used to derive drift volume
    calmsthreshold = 100;  %"calms" are drift rates less than threshold
                           % 100m^3/yr ~= 3e-6 m^3/s; 
    answer = inputdlg('Calms threshold (Qs (m^3/yr):','Drift stats',...
                                        1,{num2str(calmsthreshold)});
    if isempty(answer), return; end %user cancelled
    calmsthreshold = str2double(answer{1})/31556952; %default y2s
    
    %check input
    isprompt = true;                 %default is to prompt
    plotoption = 0;                  %default is to plot all subplots
    if nargin<3
        period = 'month';        
    else
        period = varargin{1};
        if length(varargin)>1
            inp = varargin{2};
            if islogical(inp)
                isprompt = inp;      %logical input used to set silent mode
            else
                plotoption = inp;    %numerical input used to select plot option
            end
        end
    end

    %time variables
    dur = years(tdt(end)-tdt(1));  %duration of the selected period - Tian edited Dec 2018
    tinyears = year(tdt);          %year of each record in ts
    tinmonth = month(tdt);         %month of each record in ts (1-12)
    yearrecs = tinyears(1):1:tinyears(end);
    nyr = length(yearrecs);
    dt = hours(diff(tdt));         %individual times step in hours
    delt = [dt;dt(end)];           %pad to maintain record length
    delt(delt>maxgap) = maxgap;    %remove large gaps
    vs = qs.*delt*3600;            %volume transported (m^3)
    %statistics for gaps in record, in data (NaNs) and calms 
    %drift statsistics annual and monthly
    zm = zeros(nyr,1); gapsinrecord = zm; gapsindata = zm; calms = zm;
    andrift = zm; offshore = zm; 
    anposdrift = zm; annegdrift = zm;
    mdrift = zeros(nyr*12,1); mtime = NaT(nyr*12,1); pdrift = mdrift; ndrift = mdrift;
    ii=-11;
    for i=1:nyr
        idy = find(tinyears==yearrecs(i));  %indices for year(i)
        modedt = mode(delt(idy));           %typical time step in year
        if eomday(yearrecs(i),2)==29        %end of month day in Feb is 29
            anrecs = floor(366*24/modedt);  %number of records there should be in year
        else
            anrecs = floor(365*24/modedt);  %non-leap year
        end    
        %gaps in record as a percentage: ie difference between anrecs and recs in each year
        if isempty(idy)
            gapsinrecord(i) = 100;          %no records in year(i)
        else                                %gaps as % of anrecs
            gapsinrecord(i) = (anrecs-length(idy))/anrecs*100; 
        end

        %gaps in data as a %
        %gapsindata(i)   = length(find(isnan(qs(idy))))/anrecs*100;
        gapsindata(i)   = sum(isnan(qs(idy)))/anrecs*100;  %NaN values

        %calms in either direction when >0 and <threshold
        idcalm = find(abs(qs(idy))<calmsthreshold & abs(qs(idy))>0);
        calms(i) = length(idcalm)/anrecs*100;
        %qs=0 is assigned when waves are offshore
        offshore(i) = length(find(qs(idy)==0))/anrecs*100;
        %annual drift   
        vols = vs(idy);
        andrift(i) = sum(vols,'omitnan');        
        anposdrift(i) = sum(vols(vols>0),'omitnan');
        annegdrift(i) = sum(vols(vols<0),'omitnan');
        %now get monthly sums within year
        modrift = zeros(12,1); motime = NaT(12,1); 
        posdrift = modrift; negdrift = modrift; 
        anvs = vs(idy);
        for j=1:12                     
            idm = tinmonth(idy)==j; %index for month(j)
            modrift(j) = sum(anvs(idm),'omitnan');
            vols = anvs(idm);
            posdrift(j) = sum(vols(vols>0),'omitnan');
            negdrift(j) = sum(vols(vols<0),'omitnan');
            motime(j,1) = datetime(yearrecs(i),j,15);
        end
        %alt(i) = sum(modrift);  %check that this gives annual drift
        ii = ii+12; jj=i*12;
        mdrift(ii:jj,1) = modrift;
        mtime(ii:jj,1)  = motime;
        pdrift(ii:jj,1) = posdrift;
        ndrift(ii:jj,1) = negdrift;
    end

    %convert from datenum to datetime for plotting and output
    andtn  = datetime(yearrecs,7,1);

    %set output of drift and time
    if strcmp(period,'year') || strcmp(period,'annual')
        qdrift = andrift; qtime = andtn';
    else
        qdrift = mdrift; qtime = mtime;
    end
    
    %define the gaps for plotting
    qsgaps = NaN(size(qs));
    qsgaps(isnan(qs)) = 0;
    
    stdate = datetime(yearrecs(1),1,1);
    endate = datetime(yearrecs(end),12,31);
    
    hf = figure('Name','Drift statistics','Tag','PlotFig');
    if plotoption<2
        %plot gap statistics
        s(1) = subplot(3,1,1);
        if length(yearrecs)>1
            bar(yearrecs,[gapsinrecord,gapsindata,offshore,calms],'stacked')
        else
            bar([yearrecs;NaN],[gapsinrecord,gapsindata,offshore,calms;0,0,0,0],'stacked')
        end
        xlim([yearrecs(1)-0.5,yearrecs(end)+0.5])
        xlabel('Year');
        ylabel(sprintf('Percentage gaps,\ncalms and offshore conditions'));
        calmstxt = sprintf('Calms (qs<%g m^3/s)',calmsthreshold);
        legend({'Gaps in record','Void/offshore (NaN input)','Offshore (qs=0)',calmstxt},'Location','best');
    end

    if plotoption==0 || plotoption==2
        %plot annual and monthly drift volumes
        s(2) = subplot(3,1,2);
        yyaxis left
        plot(mtime,mdrift);                                      %monthly drift
        hold on
        hp = plot([mtime(1),mtime(end)],[0,0],'--k');            %zero axis
        excludefromlegend(hp);
        plot(tdt,qsgaps,'-r','LineWidth',3);                     %gap line on zero axis
        hp = plot([tdt(1),tdt(end)],[0,0],'dr','LineWidth',1.5); %start and end of record
        excludefromlegend(hp);
        hold off
        xlim([stdate,endate])
        xlabel('Time');
        ylabel('Monthly drift volumes (m^3)');
        yyaxis right
        mver = version('-release');   %plotting datetime changes in v2017a
        if str2double(mver(1:4))<2017
            andtn = datenum(andtn); %#ok<DATNM> 
        end
        bar(andtn,andrift,'FaceColor','none','BarWidth',1);     %annual drift
        ylabel('Annual drift volumes (m^3)');
        legend({'Monthly drift','Gaps','Annual drift'},'Location','best');
    end
    
    if plotoption==0 || plotoption==3
        %plot annual positive and negative drift volumes
        s(3) = subplot(3,1,3);
        bar(andtn,anposdrift,'BarWidth',1);
        hold on 
        bar(andtn,annegdrift,'BarWidth',1);
        plot(mtime,pdrift*2,'-y')  
        plot(mtime,ndrift*2,'-g') 
        hp = plot([tdt(1),tdt(end)],[0,0],'dr','LineWidth',1.5); %start and end of record
        excludefromlegend(hp);
        hold off
        xlim([stdate,endate])
        xlabel('Time');
        ylabel('Drift volumes (m^3)');
        legend({'+ve annual: left to right','-ve annual: right to left','+ve monthly x2','-ve monthly x2'},...
            'Location','best');
    end

    if plotoption==0
        sgtitle('Drift potential')
    else
        %NB this produces a plot with multiple coordinate systems and the
        %axes cannot be copied so cannot be combined using compile_tiled_figure
        ttxt = {'Gaps and Calms','Drift rates','Directional drift rates'};
        allAxes = findall(gcf, 'type', 'axes'); % Get all axes in the figure
        delete(setdiff(allAxes,s(plotoption))); % Delete all except selected
        set(s(plotoption),'Position',[0.13 0.11 0.775 0.815]); % Default full-axes position
        title(ttxt{plotoption})
    end
    
    if isprompt
        hqd = questdlg('Save results?','Drift','Defined','All','No','No');
        if strcmp(hqd,'No') %returns total drift over period of record
            reclen = dur;                   %duration of record in years
            qdrift = (sum(qdrift))/reclen;  %returns drift rate (m3/yr) - Tian edited
            pdrift = (sum(pdrift))/reclen; 
            ndrift = (sum(ndrift))/reclen; 
            %single valued answer is returned as char in first column:
            txt1 = sprintf('Average annual drift rate %.1f m^3/y',qdrift);
            txt2 = sprintf('Average annual positive drift %.1f m^3/y',pdrift);
            txt3 = sprintf('Average annual negative drift %.1f m^3/y',ndrift);
            dst = sprintf('%s\n%s\n%s',txt1,txt2,txt3);
        elseif strcmp(hqd,'Defined')
            dsp = modelDSproperties(period);
            qtime.Format = dsp.Row.Format;   %force format to defined format
            dst = dstable(qdrift,'RowNames',qtime,'DSproperties',dsp);
        else
            adsp = modelDSproperties('year');
            andtn.Format = adsp.Row.Format;  %force format to defined format
            dst.Year = dstable(andrift,'RowNames', andtn,'DSproperties',adsp);
            mdsp = modelDSproperties('month');
            mtime.Format = mdsp.Row.Format;  %force format to defined format
            dst.Month = dstable(mdrift,'RowNames', mtime,'DSproperties',mdsp);
        end
    end
end

%%
function excludefromlegend(hp)
    set(get(get(hp,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Exclude line from legend
end
%%
function dsp = modelDSproperties(period) 
    %define a dsproperties struct and add the model metadata
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique
    
    %struct entries are cell arrays and can be column or row vectors
    switch period
        case 'year'             %Drift
            dsp.Variables = struct(...                       
                'Name',{'aQs'},...
                'Description',{'Annual alongshore drift potential'},...
                'Unit',{'m^3'},...
                'Label',{'Transport (m^3)'},...
                'QCflag',{'model'});
        case 'month'
            dsp.Variables = struct(...                       
                'Name',{'mQs'},...
                'Description',{'Monthly alongshore drift potential'},...
                'Unit',{'m^3'},...
                'Label',{'Transport (m^3)'},...
                'QCflag',{'model'});
    end
    %
    dsp.Row = struct(...
        'Name',{'Time'},...
        'Description',{'Time'},...
        'Unit',{'h'},...
        'Label',{'Time'},...
        'Format',{'dd-MM-yyyy'});        
    dsp.Dimensions = struct(...    
        'Name',{''},...
        'Description',{''},...
        'Unit',{''},...
        'Label',{''},...
        'Format',{''});    
end