function cplines = gd_setcplines(ax,promptxt,linepoints)
%%-------function help------------------------------------------------------
% NAME
%   gd_setcplines.m
% PURPOSE
%   converts a set of points to a set of lines based on NaN separators,
%   plots a graphical line, or digistise a set of points and return as a pline
% USAGE
%   cplines = gd_setcplines(ax,promptxt,linepoints);
% INPUTS
%   ax - figure axes to use to interactively select point
%   promptxt - prompt to be used for point being defined
%   linepoints - set of x,y points struct array with each line defined by a 
%                NaN separator (optional)
% OUTPUTS
%   cplines - a cell array of plines, or a new digitised pline
% NOTES
%   If linepoints is empty the user can digitise a new line of points.
%   Uses gd_setpoints to create a line of points.
%   When a graphical line is defined, the callback stores mouse event in
%   the UserData property (left or right mouse click).
% SEE ALSO
%   called in gd_sectionlines and see gd_getline
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    if nargin<3, linepoints = []; end
    [rows, cols] = size(linepoints);
    
    if isempty(linepoints)
        %no lines defined so get user to create line points
        linepoints = gd_setpoints(ax,promptxt,false);  %creates points
        nanpnts.x = NaN; nanpnts.y = NaN;              %line termination
        linepoints = [linepoints,nanpnts];             %single line
    elseif isstruct(linepoints) && length(linepoints)<2
        %not enough points for a line may be a single struct of lines
        try
            linepoints = gd_lines2points(linepoints);
        catch
            cplines = []; return;
        end
    elseif ismatrix(linepoints) && ((rows==2) || (cols==2))
        %is a matrix of xy rows 
        if cols==2, linepoints = linepoints'; end
        linepoints = gd_lines2points(linepoints);
    end
    %
    if isempty(linepoints), return; end

    idN = [0,find(isnan([linepoints(:).x]))];  
    hold on
    for i=1:length(idN)-1
        cplines{1,i} = getLine(linepoints,idN,i); %#ok<AGROW> %cell so that number of points can vary
        H = plot(ax,[cplines{1,i}.x],[cplines{1,i}.y],'-r','LineWidth',1,'Tag','mylines');
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