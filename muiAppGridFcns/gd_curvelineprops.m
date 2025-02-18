function [clinedir,ncplines,cumlen] = gd_curvelineprops(cplines,idL)
    %for each point from idL to the end use the centre-line coordinates and
    %direction to find the lengths and directions along the centre-line
    
    % cumlen is local reach specific cumulative lengths from start point
    if nargin<3,  end
    nlines = length(cplines);
    cumlen = cell(1,nlines); clinedir = cumlen; 
    nrec = length(cplines{1,1});               %length of first line
    %nrec = 0;
    j = 1;                                     %count of lines included
    for i=1:nlines
        lp = cplines{1,i};
        nl = length(lp);
        if idL>=nrec                           %start point not in line
            nrec = nrec+length(cplines{1,i+1});
            continue;
        elseif ~exist('dx','var')              %start point in line
            idl = idL-(nrec-nl);               %index of start point in line
            dx = diff([lp(idl:end-1).x]);      %omit trailing NaN
            dy = diff([lp(idl:end-1).y]);   
            ncplines{1,j} = lp(idl:end);       %#ok<AGROW> %crop line to start point
        else                                   %subsequent lines
            dx = diff([lp(1:end-1).x]);        %omit trailing NaN
            dy = diff([lp(1:end-1).y]);  
            ncplines{1,j} = lp;                %#ok<AGROW> %add subsequent lines
        end
        
        if ~isempty(dx)                        %trap single point at end of line
            %pad to make same length as lines
            dx = [dx(1),dx,dx(end)];           %#ok<AGROW> 
            dy = [dy(1),dy,dy(end)];           %#ok<AGROW> 
    
            slen = hypot(dx,dy);               %length between points
            cumlength = cumsum(slen);          %cumulative length
            cumlen{j} = [0,cumlength(2:end-1),NaN];
            theta = atan2(dy,dx);              %direction between points
            clinedir{j}(1) = theta(1);
            for k=2:length(theta)              %mean direction at point
                clinedir{j}(k) = (theta(k-1)+theta(k))/2;
            end  
            j = j+1; 
        end
    end  
end    

% FOLLOWING MOVED TO gd_linetopology so that it is added to digraph
%     if ~isempty(props)
%         %now correct lengths of reaches to be distance from mouth
%         plines = gd_cplines2plines(cplines);        
%         idN = [0,find(isnan([plines(:).x]))];
%         etable = props.topo.Edges;
%         clength = cumlen;
%         for i=1:nlines            
%             idl = etable.Weight==i;
%             idn1 = etable.Node1(idl);
%             idn2 = etable.Node2(idl);
%             idx = idn1(ismember(idn1,idn2));
%             idp1 = idx-idN(i);
%             linelen = cumlen{i}(idp1);
%             for j=1:length(linelen)
%                 idj = etable.Weight==0 & etable.Node1==idx(j);                
%                 idp2 = etable.Node2(idj)-idN(i);
%                 joinlength = getLength(plines,idp1(j),idp2);
%                 joinlength = linelen(j)+joinlength;
%                 iline = etable.Line2(idj);
%                 clength{iline} = clength{iline}+joinlength;
%             end
%         end            
%         
%         %check plot
%         cumlen = [];
%         for i=1:nlines
%             cumlen = [cumlen,clength{i}]; %#ok<AGROW> 
%         end
%         nodeids = str2double(props.topo.Nodes.Names);
%         nodedist = string(round(cumlen(nodeids)',0));        
%         %plot the resultant network
%         hg = figure('Name','Topo','Units','normalized','Tag','PlotFig');
%         ax = axes(hg);        
%         nodetable = props.topo.Nodes;
%         nodetable = addvars(nodetable,nodedist,'NewVariableNames',{'Distance'});
%         G = digraph(etable,nodetable);
%         plot(ax,G,'EdgeLabel',G.Edges.Weight,'NodeLabel',G.Nodes.Distance)
%     end
% end
% 
% %%
% function joinlength = getLength(plines,id1,id2)
%     joinpts = plines([id1,id2]);
%     dx = joinpts(1).x-joinpts(2).x;
%     dy = joinpts(1).y-joinpts(2).y;
%     joinlength = hypot(dx,dy);               %length between points
% end