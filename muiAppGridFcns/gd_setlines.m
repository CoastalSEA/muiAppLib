
function lines = gd_setlines(ax,promptxt,linepoints)
%%-------function help------------------------------------------------------
% NAME
%   gd_setlines.m
% PURPOSE
%   converts a set of points to a set of lines and based on NaN separators,
%   of digistise a set of points and return as a line
% USAGE
%   points = gd_setliness(ax,promptxt,linepoints);
% INPUTS
%   ax - figure axes to use to interactivvely select point
%   promptxt - prompt to be used for point being defined
%   linepoints - input points to be set as lines x,y struct with each line
%                define by a NaN separator (optional)
% OUTPUTS
%   lines - struct with x, y fields defining selectable lines
% NOTES
%   uses gd_setpoints to create a line of points
%   if linepoints is empty the user can digitise a new line of points
% SEE ALSO
%   called in gd_sectionlines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    [rows, cols] = size(linepoints);
    if nargin<3 || isempty(linepoints)
        %no lines defined so get user to create line points
        linepoints = gd_setpoints(ax,promptxt,false);  %creates a single line
        nanpnts.x = NaN; nanpnts.y = NaN;              %line termination
        linepoints = [linepoints,nanpnts];  
    elseif isstruct(linepoints) && length(linepoints)<2
        %not enough points for a line may be a single struct of lines
        try
            linepoints = gd_vec2pnt(linepoints);
        catch
            lines = []; return;
        end
    elseif ismatrix(linepoints) && ((rows==2) || (cols==2))
        %is a matrix of xy rows 
        if cols==2, linepoints = linepoints'; end
        linepoints = gd_vec2pnt(linepoints);
    end
    %
    if isempty(linepoints), return; end

    idN = [0,find(isnan([linepoints(:).x]))];  
    hold on
    for i=1:length(idN)-1
        lines{1,i} = getLine(linepoints,idN,i); %#ok<AGROW> %cell so that number of points can vary
        H = plot(ax,[lines{1,i}.x],[lines{1,i}.y],'-r','LineWidth',1,'Tag','mylines');
        H.ButtonDownFcn = {@lineSelected, H};
        H.UserData = int32(0);
    end
    hold off
end

%%
function lineSelected(src, evt, H)
    if evt.Button==1
        H(H==src).Color = 'g';
    elseif evt.Button==3
        H(H==src).Color = 'r';        
    end
    H(H==src).UserData = evt.Button;
end

%%
function [aline,idL] = getLine(points,idN,line)
    %get a line of points and the indices of the line to be extracted
    % idL - indices of line excluding trailing NaN
    % where idN = [0,find(isnan([points(:).x]))]; 
    idL =idN(line)+1:idN(line+1);
    aline = points(idL);
end