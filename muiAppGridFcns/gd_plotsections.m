function [xlines,hf] = gd_plotsections(grid,cplines,inp)
%
%-------function help------------------------------------------------------
% NAME
%   gd_plotsections.m
% PURPOSE
%   display grid and allow user to interactively define start and
%   end points of a section line to be plotted in a figure.   
% USAGE
%   gd_plotsections(grid,cplines,inp)
% INPUTS
%   grid - struct of x,y,z values that define grid 
%   cplines - cell array of section lines
%   inp - struct with fields for zmax, sint and method, optional if empty 
%         user is prompted to enter the values in an input dialog 
% OUTPUT
%   xlines - interpolated distances and elevations along section lines
%   hf - handle to plot of sections
% NOTES
%
% SEE ALSO
%   used in PL_Sections and PL_PlotSections
%   
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%
    %clean grid for interpolation - check for NaNs and high points
    if nargin<3 || isempty(inp)
        promptxt = {'Maximum elevation (mAD)','sampling interval (m)','Interpolatoin method'};
        defaults = {'5','10','makima'};
        input = inputdlg(promptxt,'Plot sections',1,defaults);
        if isempty(input), input = defaults; end    
        inp.zmax = str2double(input{1});
        inp.sint = str2double(input{2});
        inp.method = input{3};
    end

    grid.z(isnan(grid.z)) = inp.zmax; 
    grid.z(grid.z>inp.zmax) = inp.zmax;
    [X,Y] = meshgrid(grid.x,grid.y);

    %extract the elevations for the defined section
    nlines = length(cplines);            
    for i=1:nlines
        points = cplines{1,i};
        xlen = diff([points(1:end-1).x]);
        ylen = diff([points(1:end-1).y]);
        slen = hypot(xlen,ylen);            %length of section        
        spnts = 0:inp.sint:slen;                %points along section
        xq = points(1).x+spnts/slen*xlen;
        yq = points(1).y+spnts/slen*ylen;  

        zline = interp2(X,Y,grid.z',xq,yq,inp.method);
        if zline(1)<inp.zmax          %adjust start point if below zmax
            zline = [inp.zmax,zline]; %#ok<AGROW> 
            spnts = [0,spnts+inp.sint];
        end
        if zline(end)<inp.zmax        %adjust end point if below zmax
            zline = [zline,inp.zmax]; %#ok<AGROW> 
            spnts = [0,spnts+inp.sint];
        end

        zplines(1,i) = struct('x',[spnts,NaN],'y',[zline,NaN]); %#ok<AGROW> 
    end            
    xlines = gd_points2lines(zplines,2);
    hf = plotSections(xlines,grid.desc);
end
%%
function hf = plotSections(xlines,casedesc)
    %plot along-channel sections
    hf = figure('Name','Sections','Units','normalized',...
                     'Tag','PlotFig','Visible','on');
    ax = axes(hf);
    glines = {'-','--',':','-.'};
    hold on
    cplines = gd_plines2cplines(gd_lines2points(xlines));
    nlines = length(cplines);
    for i=1:nlines
        aline = cplines{1,i};
        spnts = [aline(:).x];
        zpnts =[ aline(:).y];
        nline = length(findobj(ax,'Tag','asection'));
        lname = sprintf('Section %d',nline+1);        
        plot(ax,spnts,zpnts,'LineStyle',glines{rem(nline,4)+1},...
              'LineWidth',1,'Tag','asection','DisplayName',lname,...
              'ButtonDownFcn',@godisplay)
    end
    hold off
    if nlines<20
        legend
    end
    title(sprintf('Along-channel sections for %s',casedesc))  
end