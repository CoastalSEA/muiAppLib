function [ax,H] = gd_plotpoints(ax,points,tagname,type,isxyz)
%
%-------function help------------------------------------------------------
% NAME
%   gd_plotpoints.m
% PURPOSE
%   plot a set of points or plines so that they can be highlighted and selected
% USAGE
%   ax = gd_plotpoints(ax,points,tagname,type)
% INPUTS
%   ax - figure axes to use to interactively select point
%   points - struct with x, y fields defining points
%   tagname - character vector of text to be used as tag for plotted points
%   type - 1=points; 2=lines; 3=point labels; 4=end points; 5=centre-line
%   isxyz - logical flag true to input z values - optional, default is false
% OUTPUTS
%   ax = updated axes
%   H - handle to graphical objects - points, lines or text
% SEE ALSO
%   called in GD_Sections and follows convention used in gd_setpoint and
%   gd_getpoint for points; and gd_setpline and gd_getpline for plines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    if nargin<5, isxyz = false; end
    if isempty(points), H = []; return; end

    switch type
        case 1      %points
            [ax,H] = setPoints(ax,points,tagname,isxyz);
        case 2      %lines
            [ax,H] = setLines(ax,points,tagname,isxyz);
        case 3      %labels
            [ax,H] = setLabels(ax,points,tagname);
        case 4      %end points
            [ax,H] = setEndPoints(ax,points,tagname);
        case 5      %centre-line markers
            [ax,H] = setClines(ax,points,tagname);
    end
end

%%
function [ax,H] = setPoints(ax,points,tagname,isxyz)
    %plot as a set of points with callback function
    points = fliplr(points);
    H = gobjects(size(points));
    hold on
    for i=1:length(points)
        H(i) = plot(ax,points(i).x,points(i).y,'ok','MarkerSize',4,...
            'MarkerFaceColor','w','Tag',tagname);
        H(i).ButtonDownFcn = {@pointSelected, H(i)};
        H(i).UserData = int32(0);
        if isxyz, ax = addZlabel(ax,points(i)); end
    end
    hold off
    %nested function-------------------------------------------------------
    function pointSelected(src, evt, H)
            if evt.Button==1
                H(H==src).Color = 'r';
            elseif evt.Button==3
                H(H==src).Color = 'k';
            end
            H(H==src).UserData = evt.Button;
        end
    end

%%
function [ax,H] = setLines(ax,plines,tagname,isxyz)
    %plot set of lines with callback function
    cplines = gd_plines2cplines(plines);
    H = gobjects(size(cplines));
    hold on
    for i=1:size(cplines,2)
        H(i) = plot(ax,[cplines{1,i}.x],[cplines{1,i}.y],'-r',...
                                             'LineWidth',1,'Tag',tagname);
        H(i).ButtonDownFcn = {@lineSelected, H(i)};
        H(i).UserData = int32(0);
        if isxyz, ax = addZlabel(ax,cplines{1,i}); end
    end
    hold off
    %nested function-------------------------------------------------------
        function lineSelected(src, evt, H)
            if evt.Button==1
                H(H==src).Color = 'g';
            elseif evt.Button==3
                H(H==src).Color = 'r';        
            end
            H(H==src).UserData = evt.Button;
        end
end

%%
function [ax,H] = setLabels(ax,pline,nlinetxt)
    %plot labels as non-pickable text
    %NB 'tagname' input variable is used to input text to display or is
    %empty to number lines in sequence
    if isempty(nlinetxt)
        nlinetxt = num2str(length(findobj(ax,'Tag','mylines'))+1); 
    end
    hpts = findobj(ax,'Tag','mypoints');    
    hold on
    if ~isempty(hpts)
        hpts(1).MarkerSize = 8;                %make start point marker bigger
        Hp = [];                               %handle not returned
    else
        Hp = plot(ax,pline(1).x,pline(1).y,'ok','MarkerSize',8,...
                                    'MarkerFaceColor','w','Tag','mytext');
    end
    Ht = text(ax,pline(1).x,pline(1).y,sprintf('%s',nlinetxt),'Clipping', 'on',...
        'HorizontalAlignment','center','PickableParts','none','FontSize',6,'Tag','mytext');
    hold off
    H = [Hp,Ht];
end

%%
function [ax,H] = setEndPoints(ax,plines,tagname)
    %plot the end points of a line as points
    idN = [0,find(isnan([plines(:).x]))];
    idS = idN(1:end-1)+1;
    idE = idN(2:end)-1;
    idp = sort([idS,idE]);
    points = plines(idp);
    [ax,H] = setPoints(ax,points,tagname,false);
end

%%
function [ax,H] = setClines(ax,plines,tagname)
    %plot the centreline as non-pickable linework
    %NB only returns handle to the end point highlight marker
    hold on
    plot(ax,[plines(:).x],[plines(:).y],'+k','MarkerSize',4,...
                              'PickableParts','none','Tag',tagname);
    plot(ax,[plines(:).x],[plines(:).y],'ok','MarkerSize',3,...
                              'PickableParts','none','Tag',tagname);
    plot(ax,[plines(:).x],[plines(:).y],'.w','MarkerSize',2,...
                              'Tag',tagname);
    H = plot(ax,[plines(1).x],[plines(1).y],'ok','MarkerSize',3,...
               'PickableParts','none','MarkerFaceColor','w','Tag',tagname);
    hold off
end

%%
function ax = addZlabel(ax,points)
    %add the values of Z as a text point label if defined
    if isfield(points,'z') && ~isempty(points.z)
        for i=1:length(points)
        text(ax,points(i).x,points(i).y,sprintf('  %.1f',points(i).z),...
                              'Color','white','FontSize',6,'Tag','ztext');
        end
    end
end