function blines = gd_setboundary(grid,paneltxt,inlines,isdel)
%
%-------function help------------------------------------------------------
% NAME
%   gd_setboundary.m
% PURPOSE
%   Accept figure to interactively generate a contour boundary
% USAGE
%   blines = gd_setboundary(grid,paneltxt,outype,isdel);
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   paneltxt- character string used for title of figure
%   inlines - struct or table of x,y vectors to be edited or the format 
%             of output if no lines are being input (see Outputs for details)
%   outype - format of output - see Outputs for details
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   blines - if lines are input, format is the same as the input, otherwise
%            outype=0: array of structs with x, y and z fields defining selected points,
%            outype=1: Nx2 or Nx3 array.
%            outype=2: struct with x, y (and z) vector fields
%            outype=3: table with x, y (and z) vector fields
%            blines = [] if user closes figure, or no points defined
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
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position,0.8,tooltips);
    ax = gd_plotgrid(h_plt,grid);

    %handle input of existing boundary lines
    if isempty(inlines) 
        outype = 2;  blines = [];   zlevel = NaN;
    elseif isnumeric(inlines) && isscalar(inlines)   
        %handle call to function with no lines
        outype = inlines;  blines = [];  zlevel = NaN;   
    else                            %plot imported lines
        blines = inlines; 
        clear inlines
        [~,outype] = gd_lines2points(blines); 
        %check that lines are terminated with a NaN
        if ~isnan(blines.x(end))
            blines = [blines;[Nan,NaN]];
        end
        if length(blines)>5000
            getdialog(sprintf('Large number of points (N=%d)\nLoading linework my take some time',length(blines)));
        end
        ax = plotLines(ax,blines);
        zlevel = [];
    end

    %if no boundary imported extract boundary using contour
    if ~isempty(zlevel)
        zlevel = setLevel();
        blines =  gd_getcontour(grid,zlevel,false);
        ax = plotLines(ax,blines);
    end

    %now work in plines rather than lines
    b_plines = gd_lines2points(blines);

    ok = 0;
    while ok<1
        waitfor(h_but,'Tag');
        if ~ishandle(h_but) %this handles the user deleting figure window
            blines = []; return;

        elseif strcmp(h_but.Tag,'Smooth')
            %smooth the shoreline
            b_plines = smoothLines(ax,b_plines);
%             clearLines(ax,{'slines'})    %remove any existing smoothed lines 
%             %get the user to define the upper limit to use for the hypsomety
%             promptxt = {'Method (0-moving av, 1-smooth)','Window size',...
%                         'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
%             defaults = {'0','10','4','10'};
%             inp = inputdlg(promptxt,'Boundary',1,defaults);
%             if ~isempty(inp)
%                 idm= logical(str2double(inp{1}));
%                 win = str2num(inp{2}); %#ok<ST2NM> allow input of vector [1x2]
%                 deg = str2double(inp{3});  
%                 npnts = str2double(inp{4});
%                 if idm==0, method = 'movmean'; else, method = 'sgolay'; end
%                 blines = gd_smoothlines(blines,method,win,deg,npnts);   
%                 hold on
%                 plot(ax,blines.x,blines.y,'-.g','LineWidth',1,'Tag','slines')
%                 hold off
%             end

        elseif strcmp(h_but.Tag,'Resample')
            %resample as selected interval
            cint = setInterval();
            if ~isempty(cint)
                b_plines = resampleLines(b_plines,cint);
                ax = plotLines(ax,blines);
            end

        elseif strcmp(h_but.Tag,'Reset')
            zlevel = setLevel();
            if ~isempty(zlevel)
                blines =  gd_getcontour(grid,zlevel,false); 

                ax = plotLines(ax,blines);
            end

        else
            ok = 1; 
            delete(h_but);   %keep figure but delete buttons
            title(ax,'')
        end
        h_but = resetbutton(ax,h_but);
    end

    %convert format of output if required
%     blines = gd_lines2points(blines);
    blines = gd_points2lines(b_plines,outype);
    
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
    %plot any points or lines that are imported 
    clearLines(ax,{'clines','slines'});  
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
    if isempty(inp), zlevel = []; return; end  %user cancelled
    zlevel = str2double(inp{1});
end

%%
function c_plines = smoothLines(ax,c_plines)
     %smooth the centre-line with option to use or revert to original
    clearLines(ax,{'smoothlines'})   %remove any existing smoothed lines
    %get the user to define the upper limit to use for the hypsomety
    promptxt = {'Method (0-moving av, 1-smooth)','Window size',...
        'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
    defaults = {'0','10','4','10'};
    inp = inputdlg(promptxt,'Boundary',1,defaults);
    if isempty(inp), return; end          %return line unchanged
    idm= logical(str2double(inp{1}));
    win = str2num(inp{2}); %#ok<ST2NM> allow input of vector [1x2]
    deg = str2double(inp{3});  
    npnts = str2double(inp{4});
    if idm==0, method = 'movmean'; else, method = 'sgolay'; end
    c_lines = gd_points2lines(c_plines,2);
    c_lines = gd_smoothlines(c_lines,method,win,deg,npnts);   
    hold on
    plot(ax,c_lines.x,c_lines.y,'-.g','LineWidth',1,'Tag','smoothlines')
    hold off
    answer = questdlg('Accept new line or retain previous version?',...
                      'Sections','Accept','Reject','Reject');

    if strcmp(answer,'Accept')
        c_plines = gd_lines2points(c_lines);
    end
    clearLines(ax,{'smoothlines','clines'}) %remove any existing centre-lines
    plotCLines(ax,c_plines);
end

%%
function cint = setInterval()
    %prompt user to set point spacing interval for the contour
    promptxt = sprintf('Sampling interval (m)');
    inp = inputdlg({promptxt},'Boundary',1,{'100'});
    if isempty(inp), cint = []; return; end  %user cancelled
    cint = str2double(inp{1});
end

%%
function clines = resampleLines(inlines,cint)
    %resample the contour lines at intervals of cint
    idN = [0;find(isnan(inlines.x))];
    [points,~] = gd_lines2points(inlines);            %convert to points
    nlines = [];
    for i=1:length(idN)-1
        cline = points(idN(i)+1:idN(i+1));            %extract line        
        cline = gd_points2lines(cline(1:end-1),1);    %convert to matrix omit trailing NaN
        if size(cline,1)>1                            %trap single point lines
            clength = sum(vecnorm(diff(cline),2,2));  %cline is a column vector [Nx2]
            cpoints = round(clength/cint);            %number of points in new line
            newcline = curvspace(cline,cpoints);
        else
            newcline = cline;
        end
        nlines = [nlines;newcline;[NaN,NaN]]; %#ok<AGROW>  
    end  
    clines.x = nlines(:,1);    %return struct of column vectors
    clines.y = nlines(:,2);
end

%%
function clearLines(ax,types)
    % delete one or line types based on Tag names
    for i=1:length(types)
        h_lines = findobj(ax,'Tag',types{i});
        delete(h_lines)
    end
end