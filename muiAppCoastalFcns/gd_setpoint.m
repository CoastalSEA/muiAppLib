function [point,H] = gd_setpoint(ax,promptxt,tagname,isxyz)
%
%-------function help------------------------------------------------------
% NAME
%   gd_setpoint.m
% PURPOSE
%   interactively create a single point on a plot and return the point
%   coordinates. Includes an option to enter an additional value at the
%   selected point (e.g. elevation).
% USAGE
%   point = gd_setpoint(ax,promptxt,tagname,isxyz)
% INPUTS
%   ax - figure axes to use to interactively select point
%   promptxt - prompt to be used for point being defined
%   tagname - character vector of text to be used as tag for plotted points
%   isxyz - logical flag true to input z values - optional, default is false
% OUTPUTS
%   point - struct with x, y fields defining added point and z if included 
%   H - handle to graphical point object
% NOTES
%   Use button 3 (right click for right handed mouse), or press Return to
%   quit and return point = [].
%   When a point is defined the graphical callback stores mouse event in
%   the UserData property (left or right mouse click).
% SEE ALSO
%   called in gd_digitisepoints and gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%  
    point = []; H = [];    
    hold on
    title(ax,promptxt)
    but=0;  
    pause(0.1)
    while (but ~= 1) %Repeat until the Left button is not clicked
        try
            [xval,yval,but] = ginput(1);  %use mouse to select points
        catch        %user closes figure while selection is in progress 
            return;
        end
        if but==3, return; end %user right clicks mouse
    end
    %
    if isempty(xval) 
        %user presses Return during ginput selection
    else
        point.x = xval;  point.y = yval; 
        if isxyz, point = setZvalue(ax,point); end
    end

    if ~isempty(point)
        [~,H] = gd_plotpoints(ax,point,tagname,1,isxyz);
%         H = plot(ax,xval,yval,'ok','MarkerSize',4,'MarkerFaceColor','w','Tag','mypoints');
%         H.ButtonDownFcn = {@pointSelected, H};
%         H.UserData = int32(0);
%         if isxyz
%             text(ax,xval,yval,sprintf('  %.1f',point.z),'Color','white',...
%                                               'FontSize',6,'Tag','ztext');
%         end
    end
    hold off
end

%%
function point = setZvalue(ax,point)
    %prompt user to define z value
    z = find_zvalue(ax,point);
    promptxt = 'Define Z value for point';
    answer = inputdlg(promptxt,'Digitize Z',1,{num2str(z)});
    if isempty(answer) 
        point.z = NaN;
    else        
        point.z = str2double(answer{1});
    end
end

%%
function z = find_zvalue(ax,point)
    %find the z-value at the selected point if there is a surface
    hsurf = findobj(ax,'Tag','PlotGrid');
    [X,Y] = meshgrid(hsurf.XData,hsurf.YData);
    z = interp2(X,Y,hsurf.CData,point.x,point.y);
end

% %%
% function pointSelected(src, evt, H)
%     if evt.Button==1
%         H(H==src).Color = 'r';
%     elseif evt.Button==3
%         H(H==src).Color = 'k';        
%     end
%     H(H==src).UserData = evt.Button;
% end