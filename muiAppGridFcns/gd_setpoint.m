function point = gd_setpoint(ax,promptxt,isxyz)
%
%-------function help------------------------------------------------------
% NAME
%   gd_setpoint.m
% PURPOSE
%   interactively select a point on a plot and return the point
%   coordinates. Includes an option to enter an additional value at the
%   selected point (e.g. elevation).
% USAGE
%   point = gd_setpoint(ax,promptxt,isxyz)
% INPUTS
%   ax - figure axes to use to interactivvely select point
%   promptxt - prompt to be used for point being defined
%   isxyz - logical flag true to input z values - optional, default is false
% OUTPUTS
%   point - struct with x and y fields defining selected point
% NOTES
%   
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
    %pause(1)
    while (but ~= 1) %Repeat until the Left button is not clicked
        try
            [xval,yval,but]=ginput(1);  %use mouse to select points
        catch        %user closes figure while selection is in progress
            point.x = NaN;   point.y = NaN;
            return;
        end
    end
    %
    if isempty(xval) %user presses Return during ginput selection
        point.x = NaN;   point.y = NaN;
    else
        point.x = xval;  point.y = yval; 
        if isxyz, point = setZvalue(point); end
    end

    if ~isnan(xval)
        plot(ax,xval,yval,'ok','MarkerSize',4,'MarkerFaceColor','w','Tag','mypoints');
        if isxyz
            text(ax,xval,yval,sprintf('  %.1f',point.z),'Color','white','FontSize',6);
        end
    end
    hold off
end
%%
function point = setZvalue(point)
    %prompt user to define z value
    promptxt = 'Define Z value for point';
    answer = inputdlg(promptxt,'Digitize Z',1,{'0'});
    if isempty(answer) 
        point.z = NaN;
    else        
        point.z = str2double(answer{1});
    end
end