function clines = gd_boundary(grid,paneltxt,outype,isdel)
%
%-------function help------------------------------------------------------
% NAME
%   gd_boundary.m
% PURPOSE
%   Accept figure to interactively generate a contour boundary
% USAGE
%   points = gd_boundary(grid,figtitle,promptxt,outype,npts,isdel);
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   paneltxt- character string used for title of figure
%   outype - format of output - see Outputs for details
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   points - if lines are input, format is the same as the input, otherwise
%            outype=0: array of structs with x, y and z fields defining selected points,
%            outype=1: Nx2 or Nx3 array.
%            outype=2: struct with x, y (and z) vector fields
%            outype=3: table with x, y (and z) vector fields
%            points = [] if user closes figure, or no points defined
% NOTES
% 
% SEE ALSO
%   called in GD_Sections and uses gd_getcontour
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
    if nargin<4
        isdel = false; 
    end
    figtitle = 'Define contour boundary';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Smooth','Resample','Reset','Use'};
    tooltips = {'Smooth the current contour lines',...
                'Resample the current contour lines at a specified interval',...
                'Reset - choose a different contour level (or NaN mask)',...
                'Use digitised points and exit. Close figure window to Quit without saving'};
    position = [0.3,0.4,0.35,0.5];
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position,tooltips);
    ax = gd_plotgrid(h_plt,grid);
    axis equal  %assume geographical projection or grid of similar dimensions
    axis tight
    zlevel = setLevel();
    clines =  gd_getcontour(grid,zlevel,false);
    ax = plotLines(ax,clines);
%     if isempty(clines)
%          delete(h_plt.Parent); return
%     else
%         ax = plotLines(ax,clines);
%     end

    ok = 0;
    while ok<1
        waitfor(h_but,'Tag');
        if ~ishandle(h_but) %this handles the user deleting figure window
            clines = []; return;
        elseif strcmp(h_but.Tag,'Smooth')
            %smooth the shoreline
            hlines = findobj(ax,'Tag','slines');    %remove any existing smoothed lines
            delete(hlines)  
            %get the user to define the upper limit to use for the hypsomety
            promptxt = {'Method (0-moving av, 1-smooth)','Window size (1 or 2 elements)',...
                        'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
            defaults = {'0','10','4','50'};
            inp = inputdlg(promptxt,'Boundary',1,defaults);
            if isempty(inp), inp=defaults; end
            idm= logical(str2double(inp{1}));
            win = str2num(inp{2}); %#ok<ST2NM> allow input of vector [1x2]
            deg = str2double(inp{3});  
            npnts = str2double(inp{4});
            if idm==0, method = 'movmean'; else, method = 'sgolay'; end
            shore = smoothLines(clines,method,win,deg,npnts);   
            hold on
            plot(ax,clines.x,shore,'-.g','LineWidth',1,'Tag','slines')
            hold off
            clines.y = shore;
        elseif strcmp(h_but.Tag,'Resample')
            %resample as selected interval
            cint = setInterval();
            clines = resampleLines(clines,cint);
            ax = plotLines(ax,clines);
        elseif strcmp(h_but.Tag,'Reset')
            zlevel = setLevel();
            clines =  gd_getcontour(grid,zlevel,false);            
            ax = plotLines(ax,clines);
        else
            ok = 1;
        end
        h_but = resetbutton(ax,h_but);
    end

    %convert format of output if required
    clines = gd_vec2pnt(clines);
    clines = gd_pnt2vec(clines,outype);
    
    %delete figure if isdel has been set by call.
    if isdel
        delete(h_plt.Parent)
    end
end

%%
function h_but = resetbutton(ax,h_but)
    %reset button and panel title string to get new input    
    if ishandle(h_but)   %check that figure has not been deleted
        ax.Title.String = 'Select a button';
        h_but.Tag = 'reset';
    end
end

%%
function ax = plotLines(ax,inlines)
    %plot any points or lines thar are imported   
    hlines = findobj(ax,'Tag','clines');    
    delete(hlines)
    hlines = findobj(ax,'Tag','slines');    
    delete(hlines)    
    % idN = find(isnan(inlines.x));
    % idN = [1,idN,length(inlines.x)];
    % i=2;
    % idl = idN(i):idN(i+1);
    if ~isempty(inlines)
        idl = 1:length(inlines.x);
        hold on
            plot(ax,inlines.x(idl),inlines.y(idl),'r','LineWidth', 1,'Tag','clines');
        hold off
    end
end

%%
function zlevel = setLevel()
    %prompt user to set the level for the contour to be extracted
    promptxt = sprintf('Elevation of shoreline contour\n(or NaN for NaN mask):');
    inp = inputdlg({promptxt},'Boundary',1,{'NaN'});
    if isempty(inp), return; end  %user cancelled
    zlevel = str2double(inp{1});
end

%%
function cint = setInterval()
    %prompt user to set point spacing interval for the contour
    promptxt = sprintf('Sampling interval (m)');
    inp = inputdlg({promptxt},'Boundary',1,{'100'});
    if isempty(inp), return; end  %user cancelled
    cint = str2double(inp{1});
end

%%
function clines = resampleLines(inlines,cint)
    %resample the contour lines at intervals of cint
    nrec = length(inlines.x);
    idN = find(isnan(inlines.x));
    idN = [1,idN,nrec]; 
    %inlines = gd_pnt2vec(inlines,1);  %convert to matrix
    [points,~] = gd_vec2pnt(inlines);  %convert to points
    nlines = [];
    for i=1:length(idN)-1
        cline = getIndex(points,idN,i);
        clength = sum(vecnorm(diff(cline),2,2));  %cline is a column vector [Nx2]
        cpoints = round(clength/cint);            %number of points in new line
        newcline = curvspace(cline,cpoints);
        nlines = [nlines;newcline;[NaN,NaN]]; %#ok<AGROW>  
    end  
    clines.x = nlines(:,1);
    clines.y = nlines(:,2);
end

%%
function smoothline = smoothLines(inlines,method,win,deg,npnts)
    %smooth each individaul line segment
    nrec = length(inlines.x);
    idN = find(isnan(inlines.x));
    idN = [1,idN,nrec]; 
    smoothline = [];
    %find each line and smooth
    for i=1:length(idN)-1
        aline = inlines.y(idN(i):idN(i+1));
        %first line has no leading nan
        if i>1, aline = aline(2:end); end
        %last line has not trailing nan
        if idN(i)==nrec, lend = nrec; else, lend = length(aline)-1; end

        if length(aline)>npnts
            if strcmp(method,'sgolay')
                aline = smoothdata(aline,method,win,'Degree',deg);
            else %use movmean
                aline = smoothdata(aline,method,win,'omitnan');
            end
        end
        smoothline = [smoothline,aline(1:lend),NaN]; %#ok<AGROW>
    end
end

%%
function [aline,idL] = getIndex(points,idN,line)
    %get the indices of the line to be extracted
    % idL - indices of line including trailing NaNs
    newpnts.x = NaN; newpnts.y = NaN;        %line termination
    idL =idN(line):idN(line+1);
    aline = points(idL);
    if isnan(aline(1).x)
        idL = idN(line)+1:idN(line+1);
        aline(1) = [];
    end
    %
    if ~isnan(aline(end).x)
        aline = [aline,newpnts];
    end
    aline = gd_pnt2vec(aline(1:end-1),1);
end