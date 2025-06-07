function res = binned_wave_climate(H,T,D)
%
%-------function help------------------------------------------------------
% NAME
%   binned_wave_climate.m
% PURPOSE
%    compute equal energy flux bins for the wave height-direction scatter
% USAGE
%   res =  binned_wave_climate(H,T,D)
%   res =  binned_wave_climate(dst)
% INPUTS: 
%   H - wave height (m) H,T and D are vector data of same length
%   T - wave period (s)
%   D - wave direction (deg)
%   Or:
%   dst - struct array of dstables containing the following:
%         dst(1) - wave height (m)
%         dst(2) - wave period (s)
%         dst(3) - wave direction (deg)
% OUTPUTS 
%   representative wave scatter bins using energy flux method
% SEE ALSO
%   Walstra D J R, Hoekstra R, Tonnon P K and Ruessink B G, 2013, Input reduction for long-term morphodynamic simulations in wave-dominated coastal settings. Coastal Engineering, 77, pp. 57-70.
%   Benedet L, Dobrochinski J P F, Walstra D J R, Klein A H F and Ranasinghe R, 2016, A morphological modeling study to compare different methods of wave climate schematization and evaluate strategies to reduce erosion losses from a beach nourishment project. Coastal Engineering, 112, pp. 69-86, 10.1016/j.coastaleng.2016.02.005.
%   Boechat Albernaz M, Brückner M Z M, van Maanen B, van der Spek A J F and Kleinhans M G, 2023, Vegetation Reconfigures Barrier Coasts and Affects Tidal Basin Infilling Under Sea Level Rise. Journal of Geophysical Research: Earth Surface, 128 (4), 10.1029/2022jf006703.
%
% Author: Ian Townend
% CoastalSEA (c) March 2023
%----------------------------------------------------------------------
%
    res = 'no output'; %null ouput required for exit in muiUserModel.setEqnData
    %extract data
    if nargin>1
        meta{1} = 'Wave height (m)';
        meta{2} = 'Wave period (s)';
        meta{3} = 'Wave direction (degTN)';        
    else
        dst = H;
        H = dst(1).data.DataTable{:,1};
        T = dst(2).data.DataTable{:,1};
        D = dst(3).data.DataTable{:,1}; 
        %create a scatter plot for wave height and direction
        meta{1} = dst(3).data.VariableLabels{1};
        meta{2} = dst(1).data.VariableLabels{1};
        meta{3} = dst(2).data.VariableLabels{1};
    end
    % scatter_plot(D,H,T,meta)

    promptxt = {'Number of direction bins','Number of Height bins'};
    defaults = {'5','4'};
    answer = inputdlg(promptxt,'Bins',1,defaults);
    if isempty(answer), return; end               %user cancelled
    Ndir = str2double(answer{1});
    Nhgt = str2double(answer{2});

    Ef = wave_energyflux(H,T,100,1025,9.81);
    meta{3} = 'Energy flux (J/ms)';
    hf = scatter_plot(D,H,Ef,meta);
    
    %get the cumulate energy by directional increments
    TotalE = sum(Ef,'omitnan');
    DirE = TotalE/Ndir;
    dirint = 1;                                   %hard coded interval
    diri = min(D):dirint:max(D);
    E_CumDir = zeros(size(diri));
    for i=1:length(diri)
        E_CumDir(i) = sum(Ef(D<=diri(i)),'omitnan');        
    end
    %bins that are empty result in duplicate values.
    [E_CumDir,idd] = unique(E_CumDir(:).');       %remove duplicates
    diri = diri(idd);
    E_dirbins =(1:Ndir)*DirE;
    dirs = interp1(E_CumDir,diri,E_dirbins);
    dirs(Ndir) = max(D);
    direction_bins(hf,dirs,max(H))

    %get the cumulative energy in each direction bin in height increments
    dirsplus = [min(D),dirs];
    heightint = 0.1;                              %hard coded interval
    hghts = zeros(Ndir,Nhgt);
    HghtE = DirE/Nhgt;
    for i=1:Ndir
        idx_dir = D>dirsplus(i) & D<=dirsplus(i+1);
        heights = 0:heightint:max(H(idx_dir));          
        E_CumHght = zeros(size(heights));
        for j=1:length(heights)
            idx_bin = idx_dir & H<=heights(j);
            E_CumHght(j) =  sum(Ef(idx_bin),'omitnan');              
        end        
        %bins that are empty result in duplicate values.
        [E_CumHght,idh] = unique(E_CumHght(:).'); %remove duplicates
        heights = heights(idh);
        E_hghtbins =(1:Nhgt)*HghtE;
        % check_plot(E_CumHght,heights);
        hghts(i,:) = interp1(E_CumHght,heights,E_hghtbins);
        hghts(i,end) = heights(end);
        clear idx_dir heights
    end
    height_bins(hf,dirsplus,hghts)

    %compute the average conditions in each bin
    hghts = [zeros(Ndir,1),hghts];
    RepH = zeros(Ndir,Nhgt); RepT = RepH; RepD = RepH;
    for i=1:Ndir
        idx_dir = D>dirsplus(i) & D<=dirsplus(i+1);
        for j=1:Nhgt
            idx_hgt = H>hghts(i,j) & H<=hghts(i,j+1);
            index = idx_dir & idx_hgt;
            RepH(i,j) = mean(H(index),'omitnan');  
            RepT(i,j) = mean(T(index),'omitnan');  
            RepD(i,j) = mean(D(index),'omitnan');  
        end
    end

    %format output table to show bins ordered to look like plot
    Height = flipud(RepH');
    Period = flipud(RepT');
    Dir = flipud(RepD');
    HeightBins = cellstr(num2str(flipud((1:Nhgt)')));
    RepTable = table(Height,Period,Dir,'RowNames',HeightBins);
    RepTable = splitvars(RepTable);
    tablefigure('Climate','Representatgive wave climate',RepTable);
end
%%
function hf = scatter_plot(D,H,V,meta)
    %create plot
    hf = figure('Name','Scatter','Units','normalized','Tag','PlotFig');
    ax = axes(hf);
    scatter(ax,D,H,[],V,'.')
    xlabel(meta{1})
    ylabel(meta{2})
    cb = colorbar;
    cb.Label.String = meta{3};
end
%%
function direction_bins(hf,dirs,maxH)
    %add the direction bins to the plot
    ax = findobj(hf,'type','axes');

    hold on
    for i=1:length(dirs)
        plot(ax,[dirs(i),dirs(i)],[0,maxH],'k','Tag','DirBins')
    end
    hold off
end
%%
function height_bins(hf,dirs,hghts)
    %add the height bins to the plot
    ax = findobj(hf,'type','axes');
    hv = findobj(ax,'Tag','DirBins');
    delete(hv)
    hold on
    for i=1:length(dirs)-1
        % dirhght = max(hghts(i,end),hghts(i+1,end));
        plot(ax,[dirs(i+1),dirs(i+1)],[0,hghts(i,end)],'k')
        for j = 1:size(hghts,2)
            plot(ax,[dirs(i),dirs(i+1)],[hghts(i,j),hghts(i,j)],'k')
        end
    end 
    hold off
end
%%
function check_plot(x,y) %#ok<DEFNU> 
    hf = figure('Name','Scatter','Units','normalized','Tag','PlotFig');
    ax = axes(hf);
    plot(ax,x,y)
end
    