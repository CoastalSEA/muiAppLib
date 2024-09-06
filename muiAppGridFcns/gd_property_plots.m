function gd_property_plots(obj,irec,src)
%
%-------function help------------------------------------------------------
% NAME
%   gd_property_plots.m
% PURPOSE
%   plots displayed on Proprety tab in ChannelForm model
% USAGE
%   gd_model_tabs(obj,irec,src)
% INPUTS
%   obj - instance of any form model that uses the GDinterface abstract class
%   irec - row/timestep in dstable to use to access grid data
%   src - handle to calling popup menu, create Figure button, or figure handle
% OUTPUT
%   generate properties table and drop down plot selection on tab or figure
% SEE ALSO
%   used to handle different types of plot accessed from the Properties tab
%   Used in ChannelForm (cf_model_tabs.m) and ModelSkill (getInletTools.m) Apps
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    
    %trap user making selection on tab after properties have been deleted
    if ~isfield(obj.Data,'Hypsometry') || isempty(obj.Data.Hypsometry)
        getdialog('No Form Properties available for selected grid')
        return;
    end

    tabsrc = src.Parent;           %handle to Form-Props tab
    
    %zprop = getDSTable(obj.Data.Grid,irec,[]);        %dstable of elevations
    pprop = getDSTable(obj.Data.Plan,irec,[]);        %dstable of half-form widths
    wprop = getDSTable(obj.Data.WaterLevels,irec,[]); %dstable of along-channel water levels
    hprop = getDSTable(obj.Data.Hypsometry,irec,[]);  %dstable of hypsometry data
    sprop = getDSTable(obj.Data.SectionProps,irec,[]);%dstable of section properties
    
    grid = getGrid(obj,irec);
    grid.desc = sprintf('%s (%s)',grid.desc,char(grid.t));
    
    %xi = zprop.Dimensions.X;               %full grid x-axis
    xi = grid.x;                           %full grid x-axis
    zcentre = hprop.Dimensions.Z;          %z-axis centre of intervals
    [indices,ixM] = gd_basin_indices(grid);
    xj = xi(indices);                      %x-axis distances (mouth-head)

%   %CF_HydroData only holds the current values and so cannot be used with 
%   %a timeseries of records. Soulution would be to save cst results in
%   %Data.WaterLevels table
%     if ~isempty(obj.RunParam) && isfield(obj.RunParam,'CF_HydroData')
%         hydobj = obj.RunParam.CF_HydroData;
%     else
%         hydobj = [];
%     end
    
    %water levels at mouth
    wl_0 = zeros(1,3);
    wl_0(1) = wprop.zhw(ixM);                %NEEDS TESTING***********************
    wl_0(2) = wprop.zmt(ixM);    
    wl_0(3) = wprop.zlw(ixM); 
 
    if strcmp(src.Tag,'FigButton')
        hfig = figure('Tag','PlotFig');
        ax = axes('Parent',hfig,'Tag','PlotFig','Units','normalized','NextPlot','add');
        src = findobj(tabsrc,'Tag','PlotList'); 
    else
        ht = findobj(tabsrc,'Type','axes'); %remove any existing plot
        delete(ht); 
        ax = axes('Parent',tabsrc,'Tag','PropertyPlot','NextPlot','add',...
           'Units','normalized','Position', [0.08,0.11,0.86,0.58]);
    end

    switch src.String{src.Value}
        case 'Hypsommetry'            
            zsurf = hprop.SurfaceArea;
            zvol = hprop.Volume;            
            hypsommetryplot(ax,zcentre,zsurf,zvol,wl_0,grid.desc)
            
        case 'Cross-sections'    
            crossectionplot(ax,grid,wl_0)   
            
        case 'Thalweg + Plan width'            
            yz = reshape(pprop.DataTable{1,:},length(grid.x),3);
            centrelineplot(ax,grid,yz);
            
        case 'Along-channel properties'
            varname = selectChannelProperty();
            [var,txt] = getChannelVariables(varname,sprop,indices,grid.Rv);
            alongchannelplot(ax,xj,var,grid,txt)

        case 'Elevation-Area histogram'
            zcentre = hprop.Dimensions.Z;
            zhist = hprop.SAfreq; %SAfreq = histogram bin count x grid area
            elevationareaplot(ax,zcentre,zhist,wl_0,grid)
            
        case 'Hydraulic depth'
            zcentre = hprop.Dimensions.Z;
            zsurf = hprop.SurfaceArea;
            zvol = hprop.Volume;          
            hydraulicdepthplot(ax,zcentre,zsurf,zvol,wl_0,grid)
            
        case 'Prism-Area ratio'
            varnames = {'CSAhw','CSAmt','CSAlw'};
            txt.names = {'P-A_h_w','P-A_m_t','P-A_l_w'};
            txt.label = 'Prism-Area ratio';
            pr = sprop.PrV(indices);
            var = zeros(length(xj),3);  
            for i=1:3
                var(:,i) = pr./sprop.(varnames{i})(indices); 
                txt.desc{i} = sprintf('%s',txt.names{i});
            end
            alongchannelplot(ax,xj,var,grid,txt)
            
        case 'a/h and Vs/Vc'
            Ahw = sprop.CSAhw(indices);
            Amt = sprop.CSAmt(indices);
            Alw = sprop.CSAlw(indices);   
            Whw = sprop.Whw(indices);
            Wmt = sprop.Wmt(indices);
            Wlw = sprop.Wlw(indices);
            amp = sprop.amp(indices); 
            % xi = wprop.Dimensions.X';
            % amp = interp1(xi,amp,xj);
            % 
            ah(:,1) = amp./(Ahw./Whw); %section specific not cumulative upstream
            ah(:,2) = amp./(Amt./Wmt);
            ah(:,3) = amp./(Alw./Wlw);
            ah(ah>1.2) = NaN; ah(ah<0) = NaN;
            VsVc = (Ahw-Alw-2*amp.*Wlw)./(Alw+amp.*Wlw);
            VsVc(VsVc>3) = NaN;
            %VsVc = movmean(VsVc,5);
            asymmratiosplot(ax,xj,ah,VsVc,grid);  
            
        case 'Hydraulics'
%             if ~isempty(hydobj) && ~isempty(hydobj.cstres)
%                 watobj = hydobj;
%             else
%                 watobj = wprop;
%             end
            hydraulicsplot(ax,wprop,xi,grid); 
            
        case 'Transgression'
            if ~isa(obj,'CF_TransModel')
                warndlg('Trangression model required for this plot');
                return;
            end
            transgressionplot(ax,obj.Data.Transgression);            
    end   
end
%%
function hypsommetryplot(ax,zcentre,zsurf,zvol,wl,casedesc)
    %plot volume and surface area as function of elevation
    ax.XDir = 'normal';
    plot(ax,zvol,zcentre,'-r');
    hold on
    plot(ax,zsurf,zcentre,'-.b');
    
    %add water levels at mouth
    plot(ax,xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(ax,xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(ax,xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    
    xlabel('Volume (m^3) and Area (m^2)'); 
    ylabel('Elevation (mAD)');
    legend('Volume','Surface area','Location','SouthEast');
    title(casedesc,'FontWeight','normal','FontSize',10);
    hold off
end
%%
function crossectionplot(ax,grid,wl)
    %plot a series of cross-sections along length of channel
    ax.XDir = 'reverse'; %to be consistent with the Q-Plot projection
    gd_dir = gd_ax_dir(grid);
    if gd_dir.x==1 || gd_dir.x==4                 
        %orientation of x-axis, x=0 is nearest the mouth if ishead=false
        zgrd = flipud(grid.z);
        ix0 = find(grid.x>=grid.xM-eps,1,'first');
    else
        %NB this finds the sections from >0 so geogrids with negative
        %coordinates may not display correctly
        zgrd = grid.z;
        ix0 = find(grid.x>=grid.xM-eps,1,'first'); 
    end
    
    noxi=length(grid.x);
    noxi = noxi-ix0;
    nx1=ix0; nx2=ix0+ceil(0.1*noxi); nx3=ix0+ceil(0.2*noxi);
    nx4=ix0+ceil(0.4*noxi); nx5=ix0+ceil(0.6*noxi); nx6=ix0+ceil(0.8*noxi);
    
    green = mcolor('green');
    orange = mcolor('orange');
    purple = mcolor('purple');    
    
    plot(ax,grid.y,zgrd(nx1,:),'-r','LineWidth',0.6);
    hold on
    plot(ax,grid.y,zgrd(nx2,:),'-.b','LineWidth',0.59);
    plot(ax,grid.y,zgrd(nx3,:),'--','Color',green,'LineWidth',0.58);
    plot(ax,grid.y,zgrd(nx4,:),'-.','Color',orange,'LineWidth',0.56);
    plot(ax,grid.y,zgrd(nx5,:),'--','Color',purple,'LineWidth',0.54);
    plot(ax,grid.y,zgrd(nx6,:),'-.m','LineWidth',0.52);
    plot(ax,grid.y,zgrd(noxi-1,:),':k');
    
    %add water levels at mouth
    plot(xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    
    xlabel('Width (m)'); 
    ylabel('Elevation (mAD)');
    hL=legend('0','0.1L','0.2L','0.4L','0.6L','0.8L','R','Location','SouthEast');
    set(hL, 'Color', 'none');
    title(grid.desc,'FontWeight','normal','FontSize',10);
    hold off
end
%%
function centrelineplot(ax,grid,yz)
    %plot centrel-line depth and width variation along channel
    ax.Position(3) = 0.8;  
    ax = set_ax_dir(ax,grid); %check orientation of x axis in grid    
    zi = grid.z(:,ceil(size(grid.z,2)/2):end); 
    xsgn = sign(grid.x(2)-grid.x(1));  %direction of axes (-ve descending)            
    %if obj.Grid.ishead  %orientation of x-axis, x=0 is nearest the mouth if ishead=false
    if xsgn<0
        zi = flipud(zi);
    end
    yyaxis left
    plot(ax,grid.x,zi(:,1),'Color','k'); 
    ylabel('Elevation (mAD)');
    yyaxis right
    plot(ax,grid.x,yz(:,1),'Color','g','LineStyle','--','LineWidth',0.6);
    hold on
    plot(ax,grid.x,yz(:,2),'Color','r','LineStyle',':','LineWidth',0.75);
    plot(ax,grid.x,yz(:,3),'Color','b','LineStyle','-.','LineWidth',0.55);
    hold off
    ylabel('Width (m)');
    hL=legend('Grid centre line','High water', 'Mean tide level','Low water',...
                                                    'Location','East');
    set(hL, 'Color', 'none');
    title(grid.desc,'FontWeight','normal','FontSize',10);    
end
%%
function alongchannelplot(ax,xj,var,grid,txt)
    %plot alnong-channel variation of selected property
    ax.Position(1) = 0.11;
    ax = set_ax_dir(ax,grid); %check orientation of channel in grid 
    plot(ax,xj,var(:,1),'DisplayName',txt.desc{1})
    hold on
    for i=2:size(var,2)
        plot(ax,xj,var(:,i),'DisplayName',txt.desc{i})
    end
    hold off
    ylabel(txt.label);
    hL = legend;                                                                                       
    set(hL, 'Color', 'none');
    title(grid.desc,'FontWeight','normal','FontSize',10);
end
%%
function elevationareaplot(ax,zcentre,zhist,wl,grid)
    %histogram of amount of surface area at each elevation
    ax.XDir = 'normal';
    barh(ax,zcentre, zhist, 'histc'); 
    hold on   
    %add water levels at mouth
    plot(xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    hold off
    xlabel('Surface area (m^2)');
    ylabel('Elevation (mAD)');
    title(grid.desc,'FontWeight','normal','FontSize',10);
end
%%
function hydraulicdepthplot(ax,zcentre,zsurf,zvol,wl,grid)
    %plot variation of hydraulic depth with elevation
    ax.XDir = 'normal';
    plot(ax,zvol./zsurf,zcentre);
    hold on   
    %add water levels at mouth
    plot(xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    hold off
    xlabel('Hydraulic depth (m)'); 
    ylabel('Elevation (mAD)');
    title(grid.desc,'FontWeight','normal','FontSize',10);
end
%%
function asymmratiosplot(ax,xj,ah,VsVc,grid)
    %plot a/h and Vs/Vc along channel
    ax.Position = [0.08,0.11,0.83,0.58];
    ax = set_ax_dir(ax,grid); %check orientation of x axis in grid 
    yyaxis left
    plot(ax,xj,ah(:,1),'Color','g','LineStyle','-');
    hold on
    plot(ax,xj,ah(:,2),'Color','r','LineStyle','-');
    plot(ax,xj,ah(:,3),'Color','b','LineStyle','-');
    hold off
    ylabel('a/h');
    yyaxis right
    plot(ax,xj,VsVc);
    ylabel('Vs/Vc');
    hL=legend('a/h High water', 'a/h Mean tide level','a/h Low water',...
                                               'Vs/Vc','Location','best');
    set(hL, 'Color', 'none');
    title(grid.desc,'FontWeight','normal','FontSize',10);
    %subtitle('Estimates are for each x-interval (not cumulative upstream)')
end
%%
function hydraulicsplot(ax,hydobj,x,grid)
    %plot the results from the CST hydraulic model or linear decay
    ax = set_ax_dir(ax,grid); %check orientation of x axis in grid 
    green = mcolor('green');
    orange = mcolor('orange');

    if isprop(hydobj,'cstres')
        zmt = hydobj.cstres.z;     %along channel mean level
        zhw = zmt+hydobj.cstres.a; %high water level
        zlw = zmt-hydobj.cstres.a; %low water level
        U = hydobj.cstres.U;       %tidal velocity amplitude
        v = hydobj.cstres.v;       %river velocity 
        d = hydobj.cstres.d;       %hydraulic depth
        incvelocity = true;  %used to plot water levels without velocities
    else
        zhw = hydobj.zhw;    %high water level
        zmt = hydobj.zmt;    %mean tide level
        zlw = hydobj.zlw;    %low water level
        incvelocity = false; %used to plot water levels without velocities
        if length(x)~=length(zhw)
            x = [min(x),max(x)];
        end
    end
    yyaxis(ax,'left')
    cla                                    %clear any existing plot lines
    plot(ax,x,zhw,'-.b','LineWidth',0.8)   %plot high water level    
    hold on
    plot(ax,x,zmt,'-r','LineWidth',1.0);   %plot time v elevation
    plot(ax,x,zlw,'-.b','LineWidth',0.8)   %plot low water level  
    ax.YLim(1) = min(zlw)-0.1;
    ax.YLim(2) = max(zhw)+0.1;
    ylabel('Elevation (mOD)'); 
    
    if incvelocity
        plot(ax,x,zmt-d,'-k','LineWidth',0.6); %hydraulic depth below mean tide level
        ax.YLim(1) = min(zmt-d);
        yyaxis(ax,'right')
        cla                                %clear any existing plot lines
        plot(ax,x,U,'--','Color',orange,'LineWidth',0.6)%plot tidal velocity
        plot(ax,x,v,'--','Color',green,'LineWidth',0.6) %plot river velocity
        ylabel('Velocity (m/s)'); 
            legend('HWL','MTL','LWL','Hydraulic depth',...
                'Tidal velocity','River velocity','Location','west');
    else
       legend('HWL','MTL','LWL','Location','west');
    end
    hold off
    ax.XLimMode = 'manual';
    ax.XLim = [min(x),max(x)];
    ax.Position(3) = ax.Position(3)*0.95;
    title(grid.desc,'FontWeight','normal','FontSize',10);
end
%%
function transgressionplot(ax,trans)
    %plot transgression time step results
    t = trans.RowNames;
    vardesc = trans.VariableDescriptions;
    ok=1; idlast = [];
    hold on
    while ok>0
        [idx,ok] = listdlg('Name','Properties','PromptString','Select variable',...                         
                           'ListSize',[220,200],'SelectionMode','single',...
                           'ListString',vardesc);
        if ok==0, continue; end
        var = trans.VariableNames{idx};
        plot(ax,t,trans.(var),'DisplayName',vardesc{idx});
        idlast = idx;
    end
    hold off
    if ~isempty(idlast)
        xlabel('Time (years)')
        ylabel(trans.VariableLabels{idlast})
        ax.Position(1) = 0.1;
        ytickformat('%3.2g')
        legend
        title(trans.Description)
    end
end
%%
function ax = set_ax_dir(ax,grid)
    %check orientation of x-axis relative to head of channel to set ax.XDir
    %this is needed because of the fixed annotation of Mouth/Head on x-axis  
    ax = gd_ax_dir(ax,grid.x);  %check whether axis direction is reversed
    %
    gd_dir = gd_ax_dir(grid);   %find orientation of channel in grid
    if gd_dir.x==1 || gd_dir.x==4                          
        ax.XLabel.String = '<= head         Distance along channel (m)        mouth =>';                                 
    else 
        ax.XLabel.String = '<= mouth         Distance along channel (m)        head =>'; 
    end
end
%%
function varname = selectChannelProperty()
    %select the property to be be used for the plot
    varnames = {'Width','Cross-sectional Area','Deepest Depth','Surface Area','Volume',...
                'Tidal prism','Gamma','Hydraulic Volumes',...
                'Hydraulic Depth'};
    answer = listdlg('Name','Properties','PromptString','Select variable',...                         
                      'ListSize',[160,140],'SelectionMode','single',...
                      'ListString',varnames);
    if isempty(answer)
        varname = [];
    else
        varname = varnames{answer};
    end
end
%%
function [var,txt] = getChannelVariables(varname,sprop,indices,Rv)
    %get the variables required for specific along-channel property plot
    xj = sprop.Dimensions.X(indices);
    switch varname
        case 'Width'
            txt.names = {'Whw','Wmt','Wlw'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar(xj,sprop,txt,indices,Rv.Wr);
        case 'Cross-sectional Area'
            txt.names = {'CSAhw','CSAmt','CSAlw'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar(xj,sprop,txt,indices,Rv.Ar);
        case 'Deepest Depth'
            txt.names = {'Dhw','Dmt','Dlw'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar_no_convergence(xj,sprop,txt,indices);
        case 'Surface Area'
            txt.names = {'Shw','Smt','Slw'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar(xj,sprop,txt,indices,0);
        case 'Volume'
            txt.names = {'Vhw','Vmt','Vlw'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar(xj,sprop,txt,indices,0);
        case 'Tidal prism'            
            txt.names = {'PrA','PrV'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar(xj,sprop,txt,indices,0);
        case 'Gamma'
            txt.names = {'Gamma'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar_no_convergence(xj,sprop,txt,indices);
        case 'Hydraulic Volumes'
            txt.names = {'Vs','Vc'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};            
            [var,txt] = getvar(xj,sprop,txt,indices,0);
        case 'Hydraulic Depth'
            txt.names = {'hyd'};
            idv = strcmp(txt.names{1},sprop.VariableNames);
            txt.label = sprop.VariableLabels{idv};
            [var,txt] = getvar_no_convergence(xj,sprop,txt,indices);
    end
    %nested function to assign variable
    function [var,txt] = getvar(xj,sprop,txt,idx,Rv)
        nvar = length(txt.names);
        var = zeros(length(xj),nvar);  
        for i=1:nvar
            var(:,i) = sprop.(txt.names{i})(idx);  %selects variable
            if isempty(Rv)                         %no river offset
                evar = var(:,i);
            else
                evar = var(:,i)-Rv;  evar(evar<0) = 0; %remove any offset at the head
            end
            Lv = -getconvergencelength(xj,evar);
            txt.desc{i} = sprintf('%s (L_v = %.0f)',txt.names{i},Lv);
        end
    end
    %
    function [var,txt] = getvar_no_convergence(xj,sprop,txt,idx)
        nvar = length(txt.names);
        var = zeros(length(xj),nvar);  
        for i=1:nvar
            var(:,i) = sprop.(txt.names{i})(idx); %selects variable
            txt.desc{i} = sprintf('%s',txt.names{i});
        end
    end
end