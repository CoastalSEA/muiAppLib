function output = gd_ax_dir(ax,x,y,ischeck)
%
%-------function help------------------------------------------------------
% NAME
%   gd_ax_dir.m
% PURPOSE
%   check direction of grid axes and reverse if descending, OR
%   find grid orientation using ishead and direction of x-axis, OR
%   check grid axis direction by prompting user
% USAGE
%    output = gd_ax_dir(ax,x,y,ischeck)
%    e.g. option 1: ax = gd_ax-dir(ax,x,y)  %axes directions match x and y
%                2: gd_dir = gd_ax_dir(grid);  %find orientation of axes in grid
%                3: output = gd_ax_dir([],x,y,[1,0]) %returns output struct
%                   with fields x (if included), y (if included) and isrev,
%                   where isrev holds a [1x2] logical array for x and y axis. 
%                   True if axis coordinates are reversed.
% INPUTS
%   ax - handle to graphic axes, a grid struct,or empty to assign axes data
%   x - column vector of x-axis data, or empty if not used
%   y - column vector of y-axis data, or empty if not used
%   ischeck - used to control options when checking grid axis direction
%             > empty or omitted - prompts user to select axis direction
%             > [1x2] logical array for x and y axis. True if axis
%               coordinates are to be reversed.
% OUTPUT
%   output - updated axes directions, or struct of x & y orientation flags, 
%            or struct of x & y data
% NOTES
%   if ax is an axes handle - checks direction of axes to be plotted
%   if ax is a grid struct - determine orientation of grid using ishead and x-direction
%   else uses ischeck to set grid direction or prompt user to select.
%
%   If only x or y passed with ax=[], only a partial struct is returned 
%   e.g. if only x input then only output.x is assigned
%
%   When used to test grid orientation: 
%   For .x 1 = x-increasing, head nearest min x, mouth nearest max x
%          2 = x-increasing, head nearest max x, mouth nearest min x
%          3 = x-decreasing, head nearest min x, mouth nearest max x
%          4 = x-decreasing, head nearest max x, mouth nearest min x
%   For .y 1 = y-increasing
%          3 = y-decreasing
% SEE ALSO
%   used in GDinterface.m and GD_ImportData.m
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    if nargin<4, ischeck = []; end
    if nargin<3, y = []; ischeck = []; end
    if nargin<2, x = []; y = []; ischeck = []; end
    
    if isgraphics(ax,'Axes') %check direction of axes to be plotted
        output = ax;
        if ~isempty(x) && issorted(x,'strictdescend')  %only change if descending
            output.XDir = 'reverse'; 
        else
            output.XDir = 'normal'; 
        end
        %
        if ~isempty(y) && issorted(y,'strictdescend')  %only change if descending
            output.YDir = 'reverse';
        else
            output.YDir = 'normal';     
        end
    elseif isstruct(ax) %determine orientation of grid using ishead and x-direction
        grid = ax;
        output.x = findaxesdirection(grid,'x');
        output.y = findaxesdirection(grid,'y');
    else                %prompt user to confirm direction of each axis
        if ~isempty(x)
            [output.x,output.isrev(1)] = checkaxesdirection(x,'X',ischeck);
        end
        %
        if ~isempty(y)
            [output.y,output.isrev(2)] = checkaxesdirection(y,'Y',ischeck);
        end
    end
end
%%
function [var,isreversed] = checkaxesdirection(input,ptxt,ischeck)
    %nested function to check the direction of the axes data
    isreversed = false;
    var = unique(input,'stable'); %stable so that orientatioin of 
                                  %grid is preserved to input direction
    dirtxt = 'ascending';
    if issorted(var,'strictdescend'), dirtxt = 'descending'; end
    if isempty(ischeck)
        promptxt = sprintf('%s-data is %s',ptxt,dirtxt);
        %prompt user to confirm direction to be used
        answer = questdlg(promptxt,'Grid-orientation','Keep','Reverse','Keep');
        if strcmp(answer,'Reverse')
            var = flipud(var);    %assumes a column vector
            isreversed = true;
        end
    elseif (ischeck(1) && strcmp(ptxt,'X')) ||...
                                      (ischeck(2) && strcmp(ptxt,'Y'))
        var = flipud(var);        %assumes a column vector
        isreversed = true;
    end
end
%%
function gd_dir = findaxesdirection(grid,anax)
    %determine orientation of grid using ishead and axis-direction of data
    %in ax2test (NB cchannel alignment currently restricted to x-axis).
    if strcmp(anax,'x')
        if grid.ishead && issorted(grid.x,'ascend') 
            gd_dir = 1; %x-increasing, head at min x, mouth at max x
        elseif ~grid.ishead && issorted(grid.x,'ascend')   
            gd_dir = 2; %x-increasing, head at max x, mouth at min x
        elseif grid.ishead && issorted(grid.x,'descend')
            gd_dir = 3; %x-decreasing, head at min x, mouth at max x
        else
            gd_dir = 4; %x-decreasing, head at max x, mouth at min x
        end  
    else
        if issorted(grid.y,'ascend')
            gd_dir = 1;
        else
            gd_dir = 3;
        end
    end
end