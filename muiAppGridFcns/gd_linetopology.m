function topo = gd_linetopology(grid,plines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_linetopology.m
% PURPOSE
%   Accept figure to interactively define line connectivity
% USAGE
%   topo = gd_linetopology(plines)
% INPUTS
%   plines - array of structs with x, y and z fields defining points in 
%            lines and separated by NaN points,
% OUTPUTS
%   topo - 
% NOTES
% 
% SEE ALSO
%   called in GD_Sections and when constructin centre-lines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    hf = figure('Name','Thalwegs','Units','normalized','Tag','PlotFig');  
    hf.Position = [0,0,1,1];
    ax = gd_plotgrid(hf,grid);
    colormap(ax,'gray');

    if isstruct(plines) && size(plines,2)==1
        plines = gd_lines2points(plines);   %ensure input is pline type
    end
    
    cplines = gd_plines2cplines(plines); 
    nlines = length(cplines);

    %no need to do anything further if only a single centre-line
    %topo = [plines(1),plines(end)];%single line topology: start-end points    
    if nlines==1  
        topo = [1,length(plines)]; %single line topology: start-end points
        return; 
    end      
    %plot the centre-lines as pickable points and lines
    ax = plotCLine(ax,plines);
    for j=1:nlines                              %call one at a time
        apline = (cplines{1,j});                   %to order numbering

        %ax = plotPoints(ax,aline,'mypoints');     %call one at a time
        %sline = gd_setcplines(ax,'',aline);       %to order numbering
        %plotLine(ax,sline{1},num2str(j));   

        ax = gd_plotpoints(ax,apline,'mypoints',1); %set points
        ax = gd_plotpoints(ax,apline,'mylines',2);  %set line  
        gd_plotpoints(ax,apline,num2str(j),3);      %set labels  
    end    
    
    firstnode = cplines{1,1}(1);
    start = 1;
    join = [];
    %select each line
    for j=1:nlines-1
        apline = cplines{1,j};
        points(j,:) = getNode(ax,j);
        joinpnts(j,:) = [firstnode,points(1)];
        firstnode = points(j,2);
        links(j,1) = find([apline(:).x]==points(1).x & [apline(:).y]==points(1).y); 
        links(j,2) = find([apline(:).x]==points(2).x & [apline(:).y]==points(2).y); 
    end
    

    %topo = subgraph(digraph(s,n),unique(s,n));
end

%%
function points = getNode(ax,nline)
    %define start and end of section line and plot it   
    prompt1 = sprintf('Select first node on line %d',nline);
    points(1) = gd_getpoint(ax,prompt1);
    prompt2 = sprintf('Select select connecting point on line %d',nline);
    points(2) = gd_getpoint(ax,prompt2);
  



    %gd_plotpoints(ax,points,[],3);      %set labels 
end

%% ------------------------------------------------------------------------
% Plotting functions
%%-------------------------------------------------------------------------
function ax = plotCLine(ax,plines)
    %plot the clentreline as non-pickable linework
    hold on
    plot(ax,[plines(:).x],[plines(:).y],'-r','LineWidth',2,...
                                   'PickableParts','none','Tag','clines');
    hold off
end