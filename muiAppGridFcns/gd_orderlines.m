function plines = gd_orderlines(ax,plines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_orderplines.m
% PURPOSE
%   amend the order of the lines in an array of plines
% USAGE
%   plines = gd_orderlines(ax,plines);;
% INPUTS
%   ax - figure axes being used to edit points
%   plines - struct of x,y vectors defining one or more lines
% OUTPUTS
%   plines - re-ordered version of input plines
% NOTES
%   
% SEE ALSO
%   called in gd_editlines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    [ax,cplines] = plotLines(ax,plines);    %plot the input lines
    
    nlines = length(cplines);
    promptxt = sprintf('Lines are numbered 1:%d. List in required order',nlines);
    defaults = {num2str(1:nlines)};
    ok = 0;
    while ok<1
        inp = inputdlg(promptxt,'Order',1,defaults);
        if isempty(inp), clearLines(ax); return; end
        order = str2num(inp{1}); %#ok<ST2NM> handles vectors
        if length(order)~=nlines, continue; end
        cplines = cplines(order);
        plines = gd_cplines2plines(cplines);
        [ax,cplines] = plotLines(ax,plines);
        answer = questdlg('Accept the revised order?','Order','Yes','No','Quit','Yes');
        if strcmp(answer,'Yes')
            ok = 1;
        elseif strcmp(answer,'Quit')
            clearLines(ax); return; 
        end
    end
    %clean up order linework
    clearLines(ax);
end
%%
function [ax,cplines] = plotLines(ax,plines)
    %plot the input lines as numbered lines
    clearGraphics(obj,{'mylines','mytext'})
    ax = clearLines(ax);
    gd_plotpoints(obj.Axes,plines,'mylines',2); 
    gd_plotpoints(obj.Axes,plines,'mypoints',1); 
    cplines = gd_plines2cplines(plines);
    hold on
    for i=1:length(cplines) 
        gd_plotpoints(obj.Axes,cplines{i},num2str(i),3); %set labels
    end
    hold off
end
%%
function ax = clearLines(ax)
    %clean up order linework
    h_pnts = findobj(ax,'Tag','mypoints');        delete(h_pnts)
    h_lines = findobj(ax,'Tag','mylines');        delete(h_lines)
    h_pnts = findobj(ax,'Tag','mylinepoints');    delete(h_pnts)
    h_txt = findobj(ax,'Tag','mytext');           delete(h_txt)
end

