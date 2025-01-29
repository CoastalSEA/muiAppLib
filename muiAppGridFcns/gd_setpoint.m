function point = gd_setpoint(ax,promptxt,isxyz)
%
%-------function help------------------------------------------------------
% NAME
%   gd_setpoint.m
% PURPOSE
%   interactively create a sinle point on a plot and return the point
%   coordinates. Includes an option to enter an additional value at the
%   selected point (e.g. elevation).
% USAGE
%   point = gd_setpoint(ax,promptxt,isxyz)
% INPUTS
%   ax - figure axes to use to interactivvely select point
%   promptxt - prompt to be used for point being defined
%   isxyz - logical flag true to input z values - optional, default is false
% OUTPUTS
%   point - struct with x, y fields defining added point and z if included 
% NOTES
%   Use button 3 or press Return to quit and return point = []
% SEE ALSO
%   called in gd_digitisepoints and gd_selectpoints
%
% Author: Ian Townend
% CoastalSEA (c) Aug 2022
%--------------------------------------------------------------------------
%  
    hold on
    title(ax,promptxt)
    but=0;  
    pause(0.1)
    while (but ~= 1) %Repeat until the Left button is not clicked
        try
            [xval,yval,but] = ginput(1);  %use mouse to select points
        catch        %user closes figure while selection is in progress
            point = [];
            return;
        end
        if but==3, point = []; return; end %user right clicks mouse
    end
    %
    if isempty(xval) %user presses Return during ginput selection
        point = [];
    else
        point.x = xval;  point.y = yval; 
        if isxyz, point = setZvalue(ax,point); end
    end

    if ~isnan(xval)
        H = plot(ax,xval,yval,'ok','MarkerSize',4,'MarkerFaceColor','w','Tag','mypoints');
        H.ButtonDownFcn = {@LineSelected, H};
        H.UserData = int32(0);
        if isxyz
            text(ax,xval,yval,sprintf('  %.1f',point.z),'Color','white',...
                                              'FontSize',6,'Tag','ztext');
        end
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

%%
function LineSelected(src, evt, H)
    if evt.Button==1
        H(H==src).Color = 'r';
    elseif evt.Button==3
        H(H==src).Color = 'k';        
    end
    H(H==src).UserData = evt.Button;
end