function ax = gd_plotpoints(ax,points,tagname,type)
%
%-------function help------------------------------------------------------
% NAME
%   gd_plotpoints.m
% PURPOSE
%   plot a set of points so that they can be highlighted and selected, or 
%   a set of lines
% USAGE
%   ax = gd_plotpoints(ax,points,tagname,type)
% INPUTS
%   ax - figure axes to use to interactively select point
%   points - struct with x, y fields defining points
%   tagname - character vector of text to be used as tag for plotted points
%   type - 1=points; 2=lines; 3=
% OUTPUTS
%   ax = updated axes
% SEE ALSO
%   called in GD_Sections and follows convention used in gd_setpoint and
%   gd_getpoint (gd_set
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    switch type
        case 1
            ax = setPoints(ax,points,tagname);
        case 2
            ax = setLines(ax,points,tagname);
        case 3
            ax = setLabels(ax,points,tagname);
    end
end

%%
function ax = setPoints(ax,points,tagname)
        %plot as a set of points
        points = fliplr(points);  
    hold on
    for i=1:length(points)
        H = plot(ax,points(i).x,points(i).y,'ok','MarkerSize',4,...
                                   'MarkerFaceColor','w','Tag',tagname);
        H.ButtonDownFcn = {@LineSelected, H};
        H.UserData = int32(0);
    end
    hold off
    %nested function-------------------------------------------------------
        function LineSelected(src, evt, H)
            if evt.Button==1
                H(H==src).Color = 'r';
            elseif evt.Button==3
                H(H==src).Color = 'k';        
            end
            H(H==src).UserData = evt.Button;
        end
end

%%
function ax = setLines(ax,plines,tagname)
    %plot set of lines
    cplines = gd_plines2cplines(plines);
    hold on
    for i=1:size(cplines,2)
        H = plot(ax,[cplines{1,i}.x],[cplines{1,i}.y],'-r',...
                                             'LineWidth',1,'Tag',tagname);
        H.ButtonDownFcn = {@lineSelected, H};
        H.UserData = int32(0);
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
function ax = setLabels(ax,line,nlinetxt)
    %plot a label
    if isempty(nlinetxt)
        nlinetxt = num2str(length(findobj(ax,'Tag','mylines'))+1); 
    end
    hpts = findobj(ax,'Tag','mypoints');    
    hold on
    hpts(1).MarkerSize = 8;                %make start point marker bigger
    text(ax,line(1).x,line(1).y,sprintf('%s',nlinetxt), 'Clipping', 'on',...
        'HorizontalAlignment','center','PickableParts','none','FontSize',6,'Tag','mytext');
    hold off
end
