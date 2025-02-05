function [subdomain,sublimitxt] = gd_subdomain(grid,isdel)
%
%-------function help------------------------------------------------------
% NAME
%   gd_subdomain.m
% PURPOSE
%   Accept figure to interactively select a subdomain of a grid
% USAGE
%   [subdomain,sublimitxt] = gd_subdomain(grid)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   isdel - logical flag is true to to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   subdomain - [min(x),max(x),min(y),max(y)] of selected domain
%               returns empty if user cancels
%   sublimitxt - text description of the subdomain
% SEE ALSO
%   called in GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jun 2022
%--------------------------------------------------------------------------
%
    if nargin<2, isdel = false; end
    x = grid.x; y = grid.y;
    gdims = gd_dimensions(grid);
     
    subdomain0 = [min(x),max(x),min(y),max(y)];    %domain limits        
    subdomain1 = [min(x)+gdims.delx,max(x)-gdims.delx, ... %initial subdomain
                  min(y)+gdims.dely,max(y)-gdims.dely];
    figtitle = sprintf('Subgrid selection');
    paneltxt = 'Subgrid definition';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Accept','Edit','Quit'};
    position = [0.2,0.4,0.4,0.4];
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position);
    ax = gd_plotgrid(h_plt,grid);
    [~,ix0,iy0] = gd_subgrid(grid,subdomain1);%muifunction
    hold on
    plot(ax,x(ix0),y(iy0),'--r','LineWidth',0.8)
    hold off
    promptxt = {'Min X','Max X','Min Y','Max Y'};
    title = 'Define subgrid';
    subdomain = subdomain1; 
    ok = 0; 
    while ok<1
        waitfor(h_but,'Tag');
        if ~ishandle(h_but) %this handles the user deleting figure window
            subdomain = []; sublimitxt = '';
            return;
        elseif strcmp(h_but.Tag,'Edit')
            %Get user to redfine subgrid                   
            sublimitxt = string(subdomain');
            asub = inputdlg(promptxt,title,1,sublimitxt);
            if isempty(asub), subdomain = []; return; end
            subdomain = cellfun(@str2double,asub)';
            
            %check for out of range and wrong order
            check(1) = any(subdomain(1,[1,3])<subdomain0(1,[1,3]));
            check(2) = any(subdomain(1,[2,4])>subdomain0(1,[2,4]));
            if gdims.xsgn<0  %x-axis is descending
                check(3) = subdomain(1,1)>subdomain(1,2);  %minX>maxX
                check(4) = subdomain(1,3)>subdomain(1,4);  %minY>maxY
            end
            
            if  any(check)  
                if any(check(1:2))
                    warndlg('Selection out of bounds. Please make a new selection');
                else
                    warndlg('Min and Max values in wrong order. Please make a new selection');
                end
                figure(h_plt.Parent); %make figure current again
                subdomain = subdomain1;
            end
            [~,ixo,iyo] = gd_subgrid(grid,subdomain);
            h_sd = findobj(ax,'Type','line');  %remove existing subdomain rectangle
            delete(h_sd);
            hold on
            plot(ax,x(ixo),y(iyo),'--r','LineWidth',0.8) %updated subdomain bounding rectangle
            hold off
            h_but.Tag = '';
        elseif strcmp(h_but.Tag,'Quit') 
            subdomain = []; sublimitxt = '';
            delete(h_plt.Parent)  %delete figure
            return;
        else
            ok = 1;                    
            delete(h_plt.Parent)  %delete figure
        end
    end
    limitxt = string(subdomain');
    sublimitxt = 'subdomain: ';
    for i=1:4
        sublimitxt = sprintf('%s%s = %s, ',sublimitxt,promptxt{i},limitxt(i));
    end
    sublimitxt = sublimitxt(1:end-2);  %remove final comma   
    
    %delete figure if isdel has been set by call.
    if isdel
        delete(h_plt.Parent)
    end
end




