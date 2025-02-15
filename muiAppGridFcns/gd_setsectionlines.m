function [s_lines,c_lines,cumlen] = gd_setsectionlines(obj,cobj,paneltxt,isdel)
%
%-------function help------------------------------------------------------
% NAME
%   gd_setsectionlines.m
% PURPOSE
%   extract the section lines that are normal to the channel centre line
%   and extend to the bounding shoreline
% USAGE
%   [s_lines,c_lines] = gd_setsectionlines(obj,cobj,paneltxt,isdel);
% INPUTS
%   obj - instance of GD_Sections with Boundary and ChannelLine defined
%   cobj - instance of EDBimport class with grid or geoimages
%   paneltxt- character string used for title of figure
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   s_lines - struct of x,y vectors defining the section lines
%   c_lines = struct of x,y vectors defining the channel centre-lines
% SEE ALSO
%   called in GD_Sections
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    if nargin<4, isdel = false; end

    figtitle = 'Define contour boundary';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Set sections','Clip sections','Adjust centre-line',...
                'Add section','Edit section','Delete section','Clear/View','Reset','Use'};
    tooltips = {'Reset the sections using current centre-line',...
                'Clip sections to the boundary line',...
                'Smooth centre-line, or reset points onto sections',...
                'Add a cross-section',...
                'Edit a point from a cross-section',...
                'Delete a cross-section',...         
                'Toggle display of section lines on and off',...
                'Reset centre-line and sections to intial state',...
                'Use digitised points and exit. Close figure window to Quit without saving'};
    % position = [0.3,0.4,0.35,0.5];
    position = [0.01,0.01,0.98,0.94];
    [h_plt,h_but] = acceptfigure(figtitle,paneltxt,tag,butnames,position,0.8,tooltips);
    ax = plotGrid(cobj,h_plt);
    msg = 'Use ''Add section'' to initialise a section';
    b_lines = obj.Boundary;          %Boundary coordinates
    hold on
    plot(ax,b_lines.x,b_lines.y,'k');
    hold off;

    c_lines = obj.ChannelLine;     %Centre line coordinates
    [c_plines,outype] = gd_lines2points(c_lines);
    ax = plotCLines(ax,c_plines);
    [s_plines,c_plines] = setInitialSections(ax,c_plines,2); %plot without labels
    h_but = resetbutton(ax,h_but);

    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window
            s_lines = []; cumlen = []; return;
        elseif strcmp(h_but.Tag,'Set sections')
            [s_plines,c_plines] = setInitialSections(ax,c_plines,1);   

        elseif strcmp(h_but.Tag,'Clip sections')
            s_plines = clipSections(ax,s_plines,b_lines,c_plines);

        elseif strcmp(h_but.Tag,'Adjust centre-line')
            %smooth the shoreline or get final centreline points that 
            %lie on the sections
            c_plines = adjustCentreLine(ax,c_plines,s_plines);

        elseif strcmp(h_but.Tag,'Add section')
            pline = addLinePoints(ax);
            if ~isempty(pline)
                gd_setcplines(ax,'',pline);
                s_plines = [s_plines,pline]; %#ok<AGROW>
            end

        elseif strcmp(h_but.Tag,'Edit section')
            if isempty(s_plines), warndlg(msg); continue; end %no sections defined
            promptxt = 'Select point to edit';
            edpoint = gd_getpoint(ax,promptxt);
            if ~isempty(edpoint)
                promptxt = 'Left click to create points, right click on any point to quit';
                newpnt = gd_setpoint(ax,promptxt,false);                
                if ~isempty(newpnt)
                   s_plines = editpoint(ax,s_plines,edpoint,newpnt);
                end                
            end

        elseif strcmp(h_but.Tag,'Delete section')
            promptxt = 'Select section to Delete, right click on any section to quit';
            deline = gd_getpline(ax,promptxt,true);
            %delpnt = gd_getpoint(ax,promptxt);   %get point to delete
            if ~isempty(deline)
                s_plines = deleteline(ax,s_plines,deline); %delete the section
            end

        elseif strcmp(h_but.Tag,'Clear/View')
            ax = toggle_view(ax,s_plines);

        elseif strcmp(h_but.Tag,'Reset')
            c_lines = obj.ChannelLine;   %Centre line coordiantes
            [c_plines,~] = gd_lines2points(c_lines);
            ax = plotCLines(ax,c_plines);
            [s_plines,c_plines] = setInitialSections(ax,c_plines,2);    
            %to plot labels on section lines and change 2 to 1 above

        else   %user accepted    
            %get cumulative length of centre-line
            cplines = gd_plines2cplines(c_plines);
            [~,~,cumlen] = clineProperties(cplines,1);
            ok = 1; 
            delete(h_but);   %keep figure but delete buttons
            title(ax,'')
        end
        h_but = resetbutton(ax,h_but);
        resetpoints(ax);
        resetlines(ax);
    end

    %convert format of output if required
    s_lines = gd_points2lines(s_plines,outype);
    c_lines = gd_points2lines(c_plines,outype);
    
    %delete figure if isdel has been set by call.
    if isdel
        delete(h_plt.Parent)
    end
end

%%
function h_but = resetbutton(ax,h_but)
    %reset button and panel title string to get new input
     if ishandle(h_but)   %check that figure has not been deleted
        ax.Title.String = 'Select a button';
        h_but.Tag = 'reset';        
     end     
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

%%
function pline = addLinePoints(ax)
    %define start and end of section line and plot it   
    prompt1 = 'Select start point';
    points(1) = gd_setpoint(ax,prompt1,false);
    prompt2 = 'Select end point';
    points(2) = gd_setpoint(ax,prompt2,false);
    gd_plotpoints(ax,points,[],3);           %set labels
    nanpnts.x = NaN; nanpnts.y = NaN;        %line termination
    pline = [points,nanpnts];   
end

%%
function plines = editpoint(ax,plines,edpoint,newpnt)
    %edit point to new value as defined in newpnt   
    h_pnts = findobj(ax,'Tag','mypoints');
    idxp = [h_pnts(:).XData]==edpoint.x & [h_pnts(:).YData]==edpoint.y; 
    
    if ~any(idxp)
       warndlg('Not current section'); return;
    else                            %call to edit point (newpnt is new point)
        h_pnts(idxp).XData = newpnt.x;
        h_pnts(idxp).YData = newpnt.y;             
        idd = [h_pnts(:).XData]==newpnt.x & [h_pnts(:).YData]==newpnt.y;
        delete([h_pnts(idd & ~idxp)]);         %remove new point       
    end 

    %find which line the point to be edited is on
    cplines = gd_plines2cplines(plines);
    [idl,idp] = findPointinLines(cplines,edpoint);
    if idl<1, return;  end          %line not found
    cplines{1,idl}(idp).x = newpnt.x; 
    cplines{1,idl}(idp).y = newpnt.y;
    plines = gd_cplines2plines(cplines);
    
    %find the indices of the line point being edited in graphics and lines
    h_lns = findobj(ax,'Tag','mylines');
    idxl = findPointinLines(h_lns,edpoint);
    idxlp = [h_lns(idxl).XData]==edpoint.x & [h_lns(idxl).YData]==edpoint.y;
    %update the point to the new location in graphics and lines
    h_lns(idxl).XData(idxlp) = newpnt.x;
    h_lns(idxl).YData(idxlp) = newpnt.y; 
    %update text marker if first point in section moved
    h_txt = (findobj(ax,'Tag','mytext'));
    if idp(1) && ~isempty(h_txt)        
        pnts = reshape([h_txt(:).Position],3,[]);
        idtl = pnts(1,:)==edpoint.x & pnts(2,:)==edpoint.y;
        nline = h_txt(idtl).String;
        delete(h_txt(idtl));
        gd_plotpoints(ax,cplines{1,idl},nline,3);  %set labels
    end
end

%%
function plines = deleteline(ax,plines,deline)
    %delete lines and points for section that includes delpnt
    delpoint = deline(1);         %first point in line
    %find which line the point to be edited is on
    cplines = gd_plines2cplines(plines);
    [idl,~] = findPointinLines(cplines,delpoint);
    if idl<1, return;  end          %line not found
    deline = cplines{1,idl};
    cplines(idl) = [];                %delete from lines cell array
    plines = gd_cplines2plines(cplines);

    %find the indices of the line point being edited in graphics and lines
    h_lns = findobj(ax,'Tag','mylines');
    idxl = findPointinLines(h_lns,delpoint);
    delete(h_lns(idxl));
    
    %update the text label
    h_txt = findobj(ax,'Tag','mytext');
    delete(h_txt(idxl));
    h_txt = flipud(findobj(ax,'Tag','mytext'));
    for i=1:length(h_txt)
        h_txt(i).String = num2str(i);
    end
    %update end points 
    h_pnts = findobj(ax,'Tag','mypoints');
    idx1 = [h_pnts(:).XData]==deline(1).x & [h_pnts(:).YData]==deline(1).y; 
    idx2 = [h_pnts(:).XData]==deline(2).x & [h_pnts(:).YData]==deline(2).y; 
    delete(h_pnts(idx1 | idx2));
end

%%
function [idl,idp] = findPointinLines(cplines,qrypoint)
    %find which line and position of point on a set of lines (can be a
    %cplines cell array containing a struct array of points, or
    %a graphical array)
    nlines = length(cplines);
    if  isgraphics(cplines,'line')   %convert graphical array to cplines
        glines = cplines; cplines = [];
        for i=1:nlines
            lp.x = glines(i).XData; 
            lp.y = glines(i).YData;
            cplines{1,i} = gd_lines2points(lp); %#ok<AGROW> 
        end   
    end

    idl = -1; idp = [];
    for i=1:nlines
        isline = gd_findline(cplines{1,i},qrypoint);
        if isline>0
            idl = i;
            idp = [cplines{1,i}.x]==qrypoint.x & [cplines{1,i}.y]==qrypoint.y;
            break
        end    
    end
end

%%
function   [s_plines,c_plines] = setInitialSections(ax,c_plines,islabel)
    %use centreline and boundary to initialise the sections lines
    % islabel determines whether sections are plotted and numbered
    clearLines(ax,{'mylines','mypoints','mytext'});
    c_cplines = gd_plines2cplines(c_plines);
    %get user to define the mouth by defining a point near to a point on the
    %channel centre-line
    promptxt = 'Left click to set mouth point, right click to quit';
    ok = 0;
    tol = (abs(diff(xlim))+abs(diff(ylim)))/2/100; %(1e4+1e5)/2/100~O[300m]
    while ok<1
        [mpnt,hp] = gd_setpoint(ax,promptxt,false);
        delete(hp)
        if isempty(mpnt)
            s_plines = []; return;   %figure has been deleted
        else
           [isNear,idL] = isPointNearLine([c_plines(:).x],[c_plines(:).y],mpnt.x,mpnt.y,tol); 
           if isNear, ok = 1; end
        end
    end
    
    %for each point from idL to the end use the centreline coordinates and
    %direction to define a section at right angles to the centreline
    [clinedir,c_cplines,~] = clineProperties(c_cplines,idL);

    %update plot of centre-line to extend from the selected mouth point
    clearLines(ax,{'clines'})
    c_plines = gd_cplines2plines(c_cplines);  %cpoints to match selected centre-line****
    ax = plotCLines(ax,c_plines);

    %generate the section lines for clinedir +pi/2 and -pi/2
    maxlen = 1000;
    s_plines = setSectionLines(c_cplines,clinedir,maxlen);

    if islabel==1                                      %labelled sections
        s_cplines = gd_plines2cplines(s_plines); 
        for j=1:length(s_cplines)                      %call one at a time
            aline = (s_cplines{1,j});                  %to order numbering
            ax = gd_plotpoints(ax,aline,'mypoints',1); %set points
            ax = gd_plotpoints(ax,aline,'mylines',2);  %set line  
            gd_plotpoints(ax,aline,num2str(j),3);      %set labels  
        end        
    elseif islabel==2                                  %unlabled sections
        ax = gd_plotpoints(ax,s_plines,'mypoints',1);  %plot as points
        gd_plotpoints(ax,s_plines,'mylines',2);        %plot as lines        
    else
        %dont plot
    end
end
%%
function [clinedir,ncplines,cumlen] = clineProperties(cplines,idL)
    %for each point from idL to the end use the centre-line coordinates and
    %direction to find the lengths anad directions along the centre-line
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
            cumlen{j} = cumsum(slen);          %cumulative length
            theta = atan2(dy,dx);              %direction between points
            clinedir{j}(1) = theta(1);
            for k=2:length(theta)              %mean direction at point
                clinedir{j}(k) = (theta(k-1)+theta(k))/2;
            end  
            j = j+1; 
        end
    end  
end

%%
function plines = setSectionLines(cplines,clinedir,maxlen)
    %create lines at right angles to the centre-line to define initial
    %sections
    plines = [];
    for i=1:length(cplines)
        sldirpos = [clinedir{i}]+pi()/2;
        sldirneg = [clinedir{i}]-pi()/2;
        lpnts = cplines{1,i};
        pospnts = sectionEndPoints(lpnts,sldirpos,maxlen);
        %gd_plotpoints(ax,pospnts,'endpoints',1);
        negpnts = sectionEndPoints(lpnts,sldirneg,maxlen);
        %gd_plotpoints(ax,negpnts,'endpoints',1);      
        %sections along ith length of centre-line
        nanpnts.x = NaN; nanpnts.y = NaN;        %line termination
        for j=1:length(pospnts)
            plines = [plines,pospnts(j),negpnts(j),nanpnts]; %#ok<AGROW> 
        end
    end
        
    %nested function-------------------------------------------------------
    function endpnts  = sectionEndPoints(lpnts,sldir,maxlen)
        npnts = length(sldir);
        endpnts = []; 
        for k=1:npnts-1
            p.x = lpnts(k).x+maxlen*cos(sldir(k));
            p.y = lpnts(k).y+maxlen*sin(sldir(k));
            endpnts = [endpnts,p]; %#ok<AGROW> 
        end
    end %------------------------------------------------------------------
end

%%
function s_plines = clipSections(ax,s_plines,b_lines,c_plines)
    %clip the sections to boundary line

    s_cplines = gd_plines2cplines(s_plines);        %cell array of lines
    nlines = length(s_cplines);   
    c_plines(isnan([c_plines(:).x])) = [];          %remove Nans so only points
    ncentres = length(c_plines);
    if nlines~=ncentres
        warndlg('Number of sections do not match number of points on centre-lines')
        return;
    else
        clines = gd_points2lines(c_plines,1);
    end
    slines = [];  spos = []; sneg = [];   %initialise loop variables
    for i=1:nlines
        sline = gd_points2lines(s_cplines{1,i},1);
        P = InterX(sline',[b_lines.x';b_lines.y']); %intersections for 1 line
        if ~isempty(P)
            cpnt = clines(i,:)';
            %find first intersections between centre-line and initial 
            %section end points
            k = 1; l = 1;
            for j=1:size(P,2)  
                %for each intersection point check whether point if on
                %positive or negative side of line
                posline = [cpnt,sline(1,:)'];   
                ison = ispointonline(posline,P(:,j),true); %true=check ends
                if ison
                    spos(:,k) = P(:,j); %#ok<AGROW> 
                    k = k+1;
                end

                negline = [cpnt,sline(2,:)'];
                ison = ispointonline(negline,P(:,j),true); %true=check ends
                if ison
                    sneg(:,l) = P(:,j); %#ok<AGROW> 
                    l = l+1;
                end
            end

            spos = findNearest(spos,cpnt);
            sneg = findNearest(sneg,cpnt);
            if ~isempty(spos) &&  ~isempty(sneg)
                slines = [slines;spos';sneg';[NaN,NaN]];     %#ok<AGROW> 
            end
            spos = []; sneg = [];   %clear for next iteration
            % hold on
            % pp = [spos,sneg];
            % plot(ax,pp(1,:),pp(2,:),'ob')
            % hold off     
        end
    end
    s_plines = gd_lines2points(slines);
    clearLines(ax,{'mylines','mypoints','mytext'});
	ax = plotPoints(ax,s_plines,'mypoints');
    %ax = gd_plotpoints(ax,s_plines,'mypoints',1);
    gd_setcplines(ax,'',s_plines); 

    %nested function-------------------------------------------------------
    function point = findNearest(points,cpnt)
        if size(points,2)>1
            p = abs(points-cpnt);
            D = hypot(p(1,:),p(2,:));
            [~,idx] = min(D);
            point = points(:,idx);
        else
            point = points;
        end
    end %------------------------------------------------------------------
end

%%
function [isNear,idL] = isPointNearLine(xData,yData,xPoint,yPoint,tol)
    % Calculate the distance from the point to the line
    distances = sqrt((xData - xPoint).^2 + (yData - yPoint).^2);
    [~,idL] = min(distances,[],'omitnan');
    isNear = any(distances < tol); % Threshold for proximity
end

%%
function c_plines= adjustCentreLine(ax,c_plines,s_plines)
    %options to smooth centre-line or reset the points onto section lines
    %when sections have been edited or deleted
    answer = questdlg('Smooth centre-line, Add points, or Reset ponts on sections?',...
                      'CentreLine','Smooth','Add point','Reset','Smooth');
    if strcmp(answer,'Smooth')            
        c_plines = smoothLines(ax,c_plines);
    elseif strcmp(answer,'Reset') && ~isempty(s_plines)
        c_plines = resetCentreLine(ax,c_plines,s_plines);
    elseif strcmp(answer,'Add point')
        promptxt = 'Left click to create point, right click on any point to quit';
        [newpnt,H] = gd_setpoint(ax,promptxt,false);  %single point
        if ~isempty(newpnt)
            c_plines = addpoints(c_plines,newpnt);  
            %clear point and replot centre-line points
            delete(H)
            clearLines(ax,{'clines'})
            plotCLines(ax,c_plines);
        end
    end
end

%%
function c_plines = smoothLines(ax,c_plines)
     %smooth the centre-line with option to use or revert to original
    clearLines(ax,{'smoothlines'})   %remove any existing smoothed lines
    %get the user to define the upper limit to use for the hypsomety
    promptxt = {'Method (0-moving av, 1-smooth)','Window size',...
        'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
    defaults = {'0','10','4','10'};
    inp = inputdlg(promptxt,'Boundary',1,defaults);
    if isempty(inp), return; end          %return line unchanged
    idm= logical(str2double(inp{1}));
    win = str2num(inp{2}); %#ok<ST2NM> allow input of vector [1x2]
    deg = str2double(inp{3});  
    npnts = str2double(inp{4});
    if idm==0, method = 'movmean'; else, method = 'sgolay'; end
    c_lines = gd_points2lines(c_plines,2);
    c_lines = gd_smoothlines(c_lines,method,win,deg,npnts);   
    hold on
    plot(ax,c_lines.x,c_lines.y,'-.g','LineWidth',1,'Tag','smoothlines')
    hold off
    answer = questdlg('Accept new line or retain previous version?',...
                      'Sections','Accept','Reject','Reject');

    if strcmp(answer,'Accept')
        c_plines = gd_lines2points(c_lines);
    end
    clearLines(ax,{'smoothlines','clines'}) %remove any existing centre-lines
    plotCLines(ax,c_plines);
end

%%
function nc_plines = resetCentreLine(ax,c_plines,s_plines)
    %update the centre-line so that the points lie on the sections and
    %return updated points and distances from mouth    
    c_cplines = gd_plines2cplines(c_plines);     %centre-line lines
    nlines = size(c_cplines,2);
    s_cplines = gd_plines2cplines(s_plines);     %cross-section lines
    nsections = size(s_cplines,2);    
    % figure;  %checkplot of sections and centrelines
    % clines = gd_points2lines(c_plines,1);   
    % hold on, plot(clines(:,1),clines(:,2)); hold off

    for j=1:nlines
        plines = [];
        cline = gd_points2lines(c_cplines{1,j},1);
        for i=1:nsections
            pline = s_cplines{1,i};              %for each cross-section
            sline = gd_points2lines(pline,1);
            % hold on, plot(sline(:,1),sline(:,2)); hold off

            P = InterX(sline',cline'); %intersections for 1 line
            %line intersection misses points at the start and end of each
            %centre-line line. Check if the point is a least on the section.
            if isempty(P)
                ison = ispointonline(sline(1:end-1,:)',cline(1,:)',1,1e3);
                if ison
                    cp = gd_lines2points(cline(1,:));
                    plines = [plines,cp];             %#ok<AGROW> 
                else                    
                    ison = ispointonline(sline(1:end-1,:)',cline(end-1,:)',1,1e3);
                    if ison
                        cp = gd_lines2points(cline(end-1,:));
                        plines = [plines,cp];         %#ok<AGROW>
                    end
                end
            else
                cp = struct('x',P(1,1),'y',P(2,1));
                plines = [plines,cp];                 %#ok<AGROW> 
            end

        end
        cplines{1,j} = plines;                        %#ok<AGROW> 
    end
    %reinstate NaN line terminators
    nanpnts.x = NaN; nanpnts.y = NaN;                 %line termination
    nc_plines = [];
    for j=1:nlines
        nc_plines = [nc_plines,cplines{1,j},nanpnts]; %#ok<AGROW> 
    end

    clearLines(ax,{'clines'}) %remove any existing centre-lines
    plotCLines(ax,nc_plines);
end

%%
function plines = addpoints(plines,point)
    %find nearest point on existing centre-line to insert new point
    dist = sqrt(([plines(:).y] - point.y.').^2 + ([plines(:).x] - point.x.').^2);
    %index for the closest points, i.e. those with min distance.
    [~,idp] = min(dist,[],2); 
    
    %find the relative position on the line
    isfirst = false; isstart = false;
    if idp==1
        isfirst = true;              %special case first point in lines
    else
        isstart = isnan(plines(idp-1).x); %point is at start of a line
    end    
    isend = isnan(plines(idp+1).x);   %end point is at end of a line
    
    %insert the new point (cases ensure NaNs are maintained)
    if isfirst
        plines = [point,plines];     %in front of all lines
    elseif isstart
        plines = [plines(1:idp-1),point,plines(idp:end)];    
    elseif isend
        plines = [plines(1:idp),point,plines(idp+1:end)];
    else
        %find the next nearest point
        xplines = plines;  
        xplines(idp) = [];
        dist = sqrt(([xplines(:).y] - point.y.').^2 + ([xplines(:).x] - point.x.').^2);
        %index for the closest points, i.e. those with min distance.
        [~,idx] = min(dist,[],2); 
        idp = min([idp,idx]);
        plines = [plines(1:idp),point,plines(idp+1:end)];
    end
end

%% ------------------------------------------------------------------------
% Plotting functions
%%-------------------------------------------------------------------------
 function ax = plotPoints(ax,points,tagname)
    %plot the imported lines
    points = fliplr(points);  
    hold on
    for i=1:length(points)
        H = plot(ax,points(i).x,points(i).y,'ok','MarkerSize',4,...
                                   'MarkerFaceColor','w','Tag',tagname);
        H.ButtonDownFcn = {@LineSelected, H};
        H.UserData = int32(0);
    end
    hold off
    %nested function
        function LineSelected(src, evt, H)
            if evt.Button==1
                H(H==src).Color = 'r';
            elseif evt.Button==3
                H(H==src).Color = 'k';        
            end
            H(H==src).UserData = evt.Button;
        end
 end

%%
function ax = plotGrid(cobj,src)
    %plot either the bathymetry grid or an image of it as a backdrop
    isgrid = false; isimage = false;
    if isfield(cobj.Data,'Grid')
        %dst = cobj.Data.Grid;
        grid = getGrid(cobj,1);                 %grid selected
        ax = gd_plotgrid(src,grid);
        hplt = findobj(ax,'Tag','PlotGrid');
        hplt.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        isgrid = true;
    elseif isfield(cobj.Data,'GeoImage')pnts
        dst = cobj.Data.GeoImage;
        im = dst.geoimage;                      %image object
        ax = axes(src);
        h_im = imagesc(ax,'XData',im.XData,'YData',im.YData,'CData',im.CData);
        set(h_im, 'AlphaData', 1-isnan(im.CData)); %set Nan values to be transparent              
        isimage = true;
    end
    %
    if isgrid || isimage 
        axis equal tight
        cb = colorbar;
        cb.Label.String = 'Elevation (mAD)';    
    end
end

%%
function ax = plotCLines(ax,points)
    %plot the clentreline as non-pickable linework
    hold on
    plot(ax,[points(:).x],[points(:).y],'+k','MarkerSize',4,...
                              'Tag','clines');
    plot(ax,[points(:).x],[points(:).y],'ok','MarkerSize',3,...
                              'Tag','clines');
    hold off
end

%%
function ax = toggle_view(ax,cplines)
    %switch a line of the lines and points defined with the start point
    %of each line emphasised with a circle marker.
    hline = findobj(ax,'Tag','mylines');
    if ~iscell(cplines)
        cplines = gd_plines2cplines(cplines);
    end
    if isempty(hline)           %toggle line and points on
        for j=1:length(cplines)                         %call one at a time
            aline = (cplines{1,j});                     %to order numbering
            ax = gd_plotpoints(ax,aline,'mypoints',1);  %set points
            ax = gd_plotpoints(ax,aline,'mylines',2);   %set line 
            ax = gd_plotpoints(ax,aline,num2str(j),3);  %set labels
        end 
    else                        %toggle line and points off
        clearLines(ax,{'mylines','mypoints','mytext'});
    end
end

%%
function clearLines(ax,types)
    % delete one or line types based on Tag names
    for i=1:length(types)
        h_lines = findobj(ax,'Tag',types{i});
        delete(h_lines)
    end
end



