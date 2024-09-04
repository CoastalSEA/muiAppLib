function res = wave_scatter(dst)
%
%-------function help------------------------------------------------------
% NAME
%   wave_scatter.m
% PURPOSE
%   plot H-T scatter diagram that uses the wave celerity based on water depth 
%   to add contours of wave steepness, plot a 3D contoured histogram or
%   various combinations of height, period, depth and steepness.
%   the  
% USAGE
%   res = wave_scatter(dst)
% INPUTS: 
%   dst - struct array of dstables containing the following:
%         dst(1) - wave height (m)
%         dst(2) - wave period (s)
%         dst(3) - water depth (m)
% OUTPUTS 
%   varioius plots of H-T scatter format contoured histograms
%
% Author: Ian Townend
% CoastalSEA (c) March 2023
%----------------------------------------------------------------------
%
    res = 'no output'; %cell ouput required by call from DataManip.createVar 

    ok = 1;
    while ok>0
        %allow user to generate various plots
        plotlist = {'Scatter plot','Scatter histogram','Steepness plots'};
        [idx,ok] = listdlg('Name','Plot options', ...
            'PromptString','Select a plot:', ...
            'SelectionMode','single', ...
            'ListSize',[250,150],...
            'ListString',plotlist);
        if ok<1, return; end

        switch plotlist{idx}
            case 'Scatter plot'
                scatter_plot(dst);
            case 'Scatter histogram'
                scatter_histogram(dst);
            case 'Steepness plots'
                steepness_plot(dst);
        end
    end
end
%%
function scatter_plot(dst)
    %create a scatter plot with steepness contours from wave Hs and Tp or
    %Tz and water depth
    H = dst(1).data.DataTable{:,1};
    T = dst(2).data.DataTable{:,1};
    d = dst(3).data.DataTable{:,1}; 

    Tmnmx = minmax(T,'omitnan');
    Tr = Tmnmx(1):0.1:Tmnmx(2);
    Steep = @(t,dep,N) celerity(t,dep).*t./N;
    %wave heights with 1/10, 1/15, 1/20 steepness at minimum and maximum depth
    [Hmn,Hmx] = wave_steepness_heights(Steep,Tr,H,d,[10,15,20]);

    %create plot
    hf = figure('Name','Rates of change','Units','normalized',...
                'Tag','PlotFig');
    ax = axes(hf);
    scatter(ax,T,H,[],d,'fill')
    xlabel(dst(2).data.VariableDescriptions{1})
    ylabel(dst(1).data.VariableDescriptions{1})
    cb = colorbar;
    cb.Label.String = dst(3).data.VariableDescriptions{1};

    mnlines = {'-k','--k','-.k'};
    mxlines = {'-r','--r','-.r'};
    mntxt = {'S=1/10 min depth','S=1/15 min depth','S=1/20 min depth'};
    mxtxt = {'S=1/10 max depth','S=1/15 max depth','S=1/20 max depth'};
    hold on
    for i=1:3
        plot(ax,Tr,Hmn(i,:),mnlines{i},'DisplayName',mntxt{i})
        plot(ax,Tr,Hmx(i,:),mxlines{i},'DisplayName',mxtxt{i})
    end
    hold off
    legend('Location','northwest')
    title(dst(1).data.Description)
end
%%
function scatter_histogram(dst)
    %create a contoured 3D histogram from Hs and Tp or Tz and water depth
    H = dst(1).data.DataTable{:,1};
    T = dst(2).data.DataTable{:,1};
    d = dst(3).data.DataTable{:,1};  

    answer = questdlg('X-Y combination?','Steepness','T-H','d-H','d-S','T-H');
    switch answer
        case 'T-H'
            v1 = T; v2 = H;
            xlabtxt = dst(2).data.VariableDescriptions{1};
            ylabtxt = dst(1).data.VariableDescriptions{1};
        case 'd-H'
            v1 = d; v2 = H;
            xlabtxt = dst(3).data.VariableDescriptions{1};
            ylabtxt = dst(1).data.VariableDescriptions{1};
        case 'd-S'
            v1 = d;
            [~,v2] = wave_steepness(H,T,d,[]);
            xlabtxt = dst(3).data.VariableDescriptions{1};
            ylabtxt = 'Wave Steepness';
    end


    %create plot
    hf = figure('Name','Rates of change','Units','normalized',...
                'Tag','PlotFig');
    ax = axes(hf);
    h = histogram2(ax,v1,v2,'Normalization','probability','FaceColor','flat',...
                                'EdgeColor','none','DisplayStyle','tile');
    cb = colorbar;
    cb.Label.String = 'Probability of occurrence';
    x = h.XBinEdges(1:end-1)+h.BinWidth(1)/2;
    y = h.YBinEdges(1:end-1)+h.BinWidth(2)/2;
    hold on
    contour(ax,x,y,h.Values','k')
    hold off
    xlabel(xlabtxt)
    ylabel(ylabtxt)
    title(dst(1).data.Description)
end
%%
function steepness_plot(dst)
    %plot wave steepness against another property
    H = dst(1).data.DataTable{:,1};
    T = dst(2).data.DataTable{:,1};
    d = dst(3).data.DataTable{:,1};
    t = dst(1).data.RowNames;

    [~,S] = wave_steepness(H,T,d,[]);

    answer = questdlg('Plot against?','Steepness','Time','S-H (d)',...
                                                      'T-H (S)','T-H (S)');

    %create plot
    hf = figure('Name','Rates of change','Units','normalized',...
                'Tag','PlotFig');
    ax = axes(hf);
    switch answer
        case 'Time'
            plot(ax,t,S)
            xlabel('Time')
        case 'S-H (d)'
            scatter(ax,S,H,10,d,'fill','DisplayName','data')
            xlabel('Wave Steepness')
            ylabel(dst(1).data.VariableDescriptions{1})
            cb = colorbar;
            cb.Label.String = dst(3).data.VariableDescriptions{1};
        case 'T-H (S)'
            Tmnmx = minmax(T,'omitnan');
            Tr = Tmnmx(1):0.1:Tmnmx(2);
            Steep = @(t,dep,N) celerity(t,dep).*t./N;
            %wave heights with 1/10, 1/15, 1/20 steepness at minimum and maximum depth
            [Hmn,Hmx] = wave_steepness_heights(Steep,Tr,H,d,[10,15,20]);

            scatter(ax,T,H,10,S,'fill')
            xlabel(dst(2).data.VariableDescriptions{1})
            ylabel(dst(1).data.VariableDescriptions{1})
            cb = colorbar;
            cb.Label.String = 'Wave Steepness';

            mnlines = {'-k','--k','-.k'};
            mxlines = {'-r','--r','-.r'};
            mntxt = {'S=1/10 min depth','S=1/15 min depth','S=1/20 min depth'};
            mxtxt = {'S=1/10 max depth','S=1/15 max depth','S=1/20 max depth'};
            hold on
            for i=1:3
                plot(ax,Tr,Hmn(i,:),mnlines{i},'DisplayName',mntxt{i})
                plot(ax,Tr,Hmx(i,:),mxlines{i},'DisplayName',mxtxt{i})
            end
            hold off
            legend('Location','northwest')
    end 
    title(dst(1).data.Description)
end
%%
function [HSmn,HSmx] = wave_steepness_heights(Steep,Tr,H,d,N)
    %wave heights for given steepness
    Hmnmx = minmax(H,'omitnan');
    dmnmx = minmax(d,'omitnan');
    for i=1:length(N)
        H4Smn = Steep(Tr,max([1,dmnmx(1)]),N(i));   %wave heigths with 1/16 steepness at minimum depth
        H4Smx = Steep(Tr,dmnmx(2),N(i));   %wave heigths with 1/16 steepness at maximum depth
        H4Smn(H4Smn>Hmnmx(2)+0.5) = NaN; %limit to range of input waves
        H4Smx(H4Smx>Hmnmx(2)+0.5) = NaN;
        HSmn(i,:) = H4Smn;
        HSmx(i,:) = H4Smx;
    end
end