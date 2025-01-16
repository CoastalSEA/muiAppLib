function res = wave_scatter_3d(H,T,D)
%
%-------function help------------------------------------------------------
% NAME
%   wave_scatter_3d.m
% PURPOSE
%   plot H-T-D scatter diagram  
% USAGE
%   res = wave_scatter_3d(dst)
% INPUTS: 
%   H - wave height (m)
%   T - wave period (s)
%   D - wave direction (dTN)
% OUTPUTS 
%   varioius plots of H-T-D scatter format contoured histograms
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%----------------------------------------------------------------------
%
    res = 'no output'; %cell ouput required by call from DataManip.createVar 
    ok = 1;
    while ok>0
        %allow user to generate various plots
        plotlist = {'Scatter plot (H,T,D)'};
        [idx,ok] = listdlg('Name','Plot options', ...
            'PromptString','Select a plot:', ...
            'SelectionMode','single', ...
            'ListSize',[250,150],...
            'ListString',plotlist);
        if ok<1, return; end

         switch plotlist{idx}
            case 'Scatter plot (H,T,D)'
                scatter_plot(H,T,D);                
        end
    end
end

%%
function scatter_plot(H,T,D)
    %create a scatter plot for height, period and direction
    
    %create plot
    hf = figure('Name','Rates of change','Units','normalized',...
                'Tag','PlotFig');
    ax = axes(hf);

    scatter3(ax,T,D,H,[],H,'.');
    xlabel('Wave period (s)')
    ylabel('Wave direction (dTN)')
    zlabel('Wave height (m)')
    cmap = cmap_selection;
    colormap(cmap);
end