function [cumlen,G,hf,hg] = gd_linetopology(grid,plines)
%
%-------function help------------------------------------------------------
% NAME
%   gd_linetopology.m
% PURPOSE
%   figure to interactively define line connectivity
% USAGE
%   [cumlen,G,hf,hg] = gd_linetopology(grid,plines)
% INPUTS
%   plines - array of structs with x, y and z fields defining points in 
%            lines and separated by NaN points,
% OUTPUTS
%   cumlen - cumulative lengths along lines from first point in line set
%   G - directed graph of the channl network
%   hf - handle to figure used to define links;
%   hg - handle to figure of the resultant network
% NOTES
%   the Nodes are labelled with the point indices for the centre-line and
%   the Edge weights are the line number of each connecting reach. 
%   Connecting links between lines have a weight of 0.
% SEE ALSO
%   called in GD_Sections and when constructin centre-lines
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    hf = figure('Name','Topology','Units','normalized','Tag','PlotFig');  
    hf.Position = [0,0.03,1,0.93];
    ax = plotBackGround(hf,grid);

    if isstruct(plines) && size(plines,2)==1
        plines = gd_lines2points(plines);   %ensure input is pline type
    end
    
    cplines = gd_plines2cplines(plines); 
    nlines = length(cplines);

    %no need to do anything further if only a single centre-line  
    if nlines==1  
        topo = [1,length(plines)-1]; %single line topology: start-end points
        G = setChannelGraph(topo,1);
        return; 
    end   

    %plot the centre-lines as pickable points and lines
    for j=1:nlines                                  %call one at a time
        apline = (cplines{1,j});                    %to order numbering  
        ax = gd_plotpoints(ax,apline,'mylines',2);  %set line  
        % ax = gd_plotpoints(ax,apline,'mypoints',1); %set points
        ax = gd_plotpoints(ax,apline,num2str(j),3);      %set labels  
    end    

    %%get user to define links and set up graph matrix
    idN = find(isnan([plines(:).x]));
    ids = [1,idN(1:end-1)+1];
    ide = idN-1;

    for j=1:nlines-1
        [topo(j,:),lineids(j,:)] = getNode(ax,plines,[1,idN]);
        resetpoints(ax);
        resetlines(ax);
    end
    getdialog('Connections defined')
    
    %now sort the connections into a directed graph and plot the result
    [topo,idp] = sortrows(fliplr(topo));         %order graph matrix
    lineids = lineids(idp,:);
    edges = zeros(1,size(topo,1));
    for i=1:length(ide)
        ends = [1,ide];
        idx = topo(:,1)>ends(i) & topo(:,1)<ends(i+1);
        if any(idx)
            addpnts = [ids(i);topo(idx,1);ide(i)];
            addpnts = [addpnts(1:end-1),addpnts(2:end)];
            topo = [topo;addpnts]; %#ok<*AGROW> 
            np = size(addpnts,1);
            edges = [edges,ones(1,np)*i];
            lineids = [lineids;ones(np,2)*i];
        else
            topo = [topo;[ids(i),ide(i)]];
            edges = [edges,i]; 
            lineids = [lineids;[i,i]];
        end
    end
    [cumlen,G,hg] = setChannelGraph(topo,edges,lineids,cplines);
end

%%
function [cumlen,G,hg] = setChannelGraph(topo,edges,lineids,cplines)
    %use the start and end nodes to create a directed graph of the network
    [topo,idd] = sortrows(topo);
    lineids = lineids(idd,:);
    names = string(unique(topo));
    G = subgraph(digraph(topo(:,1),topo(:,2)),unique(topo));
    G.Nodes.Names = names;
    G.Edges.Weight = edges(idd)';
    %rebuild digraph with addition of line no. data
    edgeTable = G.Edges;
    edgeTable = addvars(edgeTable,topo(:,1),topo(:,2),lineids(:,1),lineids(:,2),...
                     'NewVariableNames',{'Node1','Node2','Line1','Line2'});
    G = digraph(edgeTable,G.Nodes);

    %add the cumulative lenghts from mouth to the nodes
    cumlen = alongLineLengths(cplines,G);
    nodeids = str2double(G.Nodes.Names);
    nodedist = string(round(cumlen(nodeids)',0));
    %plot the resultant network
    nodetable = G.Nodes;
    nodetable = addvars(nodetable,nodedist,'NewVariableNames',{'Distance'});
    G = digraph(edgeTable,nodetable);

    %plot the resultant network
    hg = figure('Name','Topo','Units','normalized','Tag','PlotFig');
    ax = axes(hg);
    nlabel = strcat(cellstr(G.Nodes.Names),{' ('},cellstr(G.Nodes.Distance),{'m)'});
    plot(ax,G,'EdgeLabel',G.Edges.Weight,'NodeLabel',nlabel)
end

%%
function [idp,idl] = getNode(ax,plines,idN)
    %define start and end of section line and plot it
    % idp - [1x2] indices of line end point and join point
    % idl - index of line being joined
    orange = [0.8500 0.3250 0.0980];
    prompt1 = 'Select a branch line';
    ok = 0;    
    while ok<1        
        pline = gd_getpline(ax,prompt1,'mylines');
        if ~isempty(pline), ok = 1; end
    end
    epoint = pline(1);
    %assign first index point
    idp(1) = find([plines(:).x]==epoint.x & [plines(:).y]==epoint.y); 
    %find the line number
    nline = find(idp(1)>idN,1,'last');
    [ax,Hpts] = gd_plotpoints(ax,plines,'mypoints',1); %set points


    prompt2 = sprintf('Select point on a line that branch %d connects to',nline);
    ok = 0;
    while ok<1
        jpoint = gd_getpoint(ax,prompt2,'mypoints');
        if ~isempty(jpoint)
            %highlight branch line label
            h_txt = findobj(ax,'Type','Text','Tag','mytext');
            pnts = reshape([h_txt(:).Position],3,[]);
            idx = pnts(1,:)==epoint.x & pnts(2,:)==epoint.y;
            h_txt(idx).Color = orange;
            delete(Hpts)
            %plot join point and add label
            hold on
            H = plot(ax,jpoint.x,jpoint.y,'ok','MarkerSize',8,...
                                    'MarkerFaceColor','w','Tag','mytext');
            % H.MarkerSize = 10;    %highlight selected point
            answer = questdlg('Accept point','Topo','Yes','No','Yes');
            if strcmp(answer,'Yes')
                text(ax,jpoint.x,jpoint.y,num2str(nline),'Clipping', 'on',...
                         'HorizontalAlignment','center','Color',orange,...
                                             'FontSize',6,'Tag','mytext');
                ok = 1; 
            else
                [ax,Hpts] = gd_plotpoints(ax,plines,'mypoints',1); %set points
                delete(H)
            end
            hold off
        end
    end
    %assign send index point
    idp(2) = find([plines(:).x]==jpoint.x & [plines(:).y]==jpoint.y); 
    %indices of lines being joined
    idl =  [find(idp(2)>idN,1,'last'),nline];
end

%%
function cumlen = alongLineLengths(cplines,G)
    [~,~,reachlen] = gd_curvelineprops(cplines,1);
    %now correct lengths of reaches to be distance from mouth
    plines = gd_cplines2plines(cplines);        
    idN = [0,find(isnan([plines(:).x]))];
    etable = G.Edges;
    clength = reachlen;
    nlines = length(cplines);
    for i=1:nlines            
        idl = etable.Weight==i;
        idn1 = etable.Node1(idl);
        idn2 = etable.Node2(idl);
        idx = idn1(ismember(idn1,idn2));
        idp1 = idx-idN(i);
        linelen = reachlen{i}(idp1);
        for j=1:length(linelen)
            idj = etable.Weight==0 & etable.Node1==idx(j);                
            idp2 = etable.Node2(idj)-idN(i);
            joinlength = getLength(plines,idp1(j),idp2);
            joinlength = linelen(j)+joinlength;
            iline = etable.Line2(idj);
            clength{iline} = clength{iline}+joinlength;
        end
    end            

    cumlen = [];
    for i=1:nlines
        cumlen = [cumlen,clength{i}]; 
    end
end

%%
function joinlength = getLength(plines,id1,id2)
    joinpts = plines([id1,id2]);
    dx = joinpts(1).x-joinpts(2).x;
    dy = joinpts(1).y-joinpts(2).y;
    joinlength = hypot(dx,dy);               %length between points
end

%%
function newpnts = resetpoints(ax)
    %gd_getpoint sets point UserData and color when clicked on. Reset in 
    %case user clicked on points without making an action selection
    newpnts = [];    %reset newpnts in case user Quits when adding, etc
    h_pnts = findobj(ax,'Tag','mypoints');
    if isempty(h_pnts), return; end
    idx = [h_pnts.UserData]>0;
    if any(idx)
        [h_pnts(idx).UserData] = deal(int32(0));
        [h_pnts(idx).Color] = deal([0,0,0]);   
    end    
end

%%
function resetlines(ax)
    %gd_getpoint sets point UserData and color when clicked on. Reset in 
    %case user clicked on points without making an action selection
    h_lines = findobj(ax,'Tag','mylines');
    if isempty(h_lines), return; end
    idx = [h_lines.UserData]>0;
    if any(idx)
        [h_lines(idx).UserData] = deal(int32(0));
        [h_lines(idx).Color] = deal([1,0,0]);
    end    
end

%% ------------------------------------------------------------------------
% Plotting functions
%%-------------------------------------------------------------------------
function ax = plotCLine(ax,plines)
    %plot the centreline as non-pickable linework
    hold on
    plot(ax,[plines(:).x],[plines(:).y],'-r','LineWidth',2,...
                                   'PickableParts','none','Tag','clines');
    hold off
end

%%
function ax = plotBackGround(hf,grid)
    %determin if input is a grid or an image and plot as background
    if isgraphics(grid,'axes') %if graphics already exists
        h_grid = findobj(grid.Children,'Tag','PlotGrid');
        clear grid
        ax = axes(hf);
        copyobj(h_grid,ax);
    elseif isstruct(grid)                       %if called with a grid then plot the grid
        ax = gd_plotgrid(hf,grid);
        hplt = findobj(ax,'Tag','PlotGrid');
        hplt.Annotation.LegendInformation.IconDisplayStyle = 'off';      
    elseif isa(grid,'dstable')
        ax = axes(hf);
        %image is passed to class as the 'grid' input variable
        img = grid.geoimage;     %image object
        h_im = imagesc(ax,'XData',img.XData,'YData',img.YData,'CData',img.CData);
        set(h_im, 'AlphaData', 1-isnan(img.CData)); %set Nan values to be transparent
    else
        warndlg('Unknown grid/image type'); return;
    end
    axis equal tight
    colormap(ax,'gray');
end