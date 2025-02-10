function [slines,clines] = gd_sectionlines(obj,cobj,paneltxt,isdel)
%
%-------function help------------------------------------------------------
% NAME
%   gd_sectionlines.m
% PURPOSE
%   extract the section lines that are normal to the channel centre line
%   and extend to the bounding shoreline
% USAGE
%   [slines,clines] = gd_sectionlines(obj,cobj,paneltxt,isdel);
% INPUTS
%   obj - instance of GD_Sections with Boundary and ChannelLine defined
%   cobj - instance of EDBimport class with grid or geoimages
%   paneltxt- character string used for title of figure
%   isdel - logical flag true to delete figure on completion - optional, 
%           default is false
% OUTPUTS
%   slines - struct of x,y vectors defining the section lines
%   clines = struct of x,y vectors defining the channel centre-lines
% SEE ALSO
%   called in GD_Sections
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%
    if nargin<4, isdel = false; end
   % nanpnts.x = NaN; nanpnts.y = NaN;        %line termination
    figtitle = 'Define contour boundary';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Set sections','Clip sections','Smooth centre-line',...
                'Add section','Edit point','Delete section','Clear/View','Reset','Use'};
    tooltips = {'Reset the sections using current centre-line',...
                'Clip sections to the boundary line',...
                'Smooth the centre-line',...
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
    bline = obj.Boundary;          %Boundary coordinates
    hold on
    plot(ax,bline.x,bline.y,'k');
    hold off;

    clines = obj.ChannelLine;     %Centre line coordinates
    [cpoints,outype] = gd_vec2pnt(clines);
    ax = plotCLine(ax,cpoints);
    [slines,clines] = setInitialSections(ax,cpoints,2); %plot without labels
    h_but = resetbutton(ax,h_but);

    ok = 0;
    while ok<1
        waitfor(h_but,'Tag')
        if ~ishandle(h_but) %this handles the user deleting figure window
            slines = []; return;
        elseif strcmp(h_but.Tag,'Set sections')
            [slines,clines] = setInitialSections(ax,cpoints,1);   

        elseif strcmp(h_but.Tag,'Clip sections')


        elseif strcmp(h_but.Tag,'Smooth centre-line')
            %smooth the shoreline
            clearLines(ax,{'clines'})   %remove any existing smoothed lines
            %get the user to define the upper limit to use for the hypsomety
            promptxt = {'Method (0-moving av, 1-smooth)','Window size',...
                'Savitzky-Golay degree (<window)','Mininum points in line to smooth'};
            defaults = {'0','10','4','10'};
            inp = inputdlg(promptxt,'Boundary',1,defaults);
            if ~isempty(inp)
                idm= logical(str2double(inp{1}));
                win = str2num(inp{2}); %#ok<ST2NM> allow input of vector [1x2]
                deg = str2double(inp{3});  
                npnts = str2double(inp{4});
                if idm==0, method = 'movmean'; else, method = 'sgolay'; end
%                 lines = gd_pnt2vec(clines,2);
                lines = gd_smoothlines(clines,method,win,deg,npnts);   
                hold on
                plot(ax,clines.x,clines.y,'-.g','LineWidth',1,'Tag','slines')
                hold off
            end

        elseif strcmp(h_but.Tag,'Add section')
            newpoints = addLinePoints(ax);
            newline = gd_setlines(ax,'',newpoints);
            slines = [slines,newline]; %#ok<AGROW>

        elseif strcmp(h_but.Tag,'Edit point')
            if isempty(slines), warndlg(msg); continue; end %no sections defined
            promptxt = 'Select point to edit';
            edpoint = gd_getpoint(ax,promptxt);
            if ~isempty(edpoint)
                promptxt = 'Left click to create points, right click on any point to quit';
                newpnt = gd_setpoint(ax,promptxt,false);
                slines = editpoint(ax,slines,edpoint,newpnt);
            end

        elseif strcmp(h_but.Tag,'Clear/View')
            ax = toggle_view(ax,slines);

        elseif strcmp(h_but.Tag,'Delete section')
            promptxt = 'Select section to Delete, right click on any section to quit';
            deline = gd_getline(ax,promptxt,true);
            %delpnt = gd_getpoint(ax,promptxt);   %get point to delete
            if ~isempty(deline)
                slines = deleteline(ax,slines,deline); %delete the section
            end

        elseif strcmp(h_but.Tag,'Reset')
            clines = obj.ChannelLine;   %Centre line coordiantes
            [cpoints,outype] = gd_vec2pnt(clines);
            ax = plotCLine(ax,cpoints);
            [slines,clines] = setInitialSections(ax,cpoints,2);    
            %to plot labels on section lines and change 2 to 1 above

        else   %user accepted
            %get final centreline points that lie on the sections and
            %update ChannelLine and SectionLines Properties

            ok = 1;     delete(h_but);   %keep figure but delete buttons

        end
        h_but = resetbutton(ax,h_but);
        resetpoints(ax);
        resetlines(ax);
    end

    %convert format of output if required
    slines = gd_pnt2vec(slines,outype);
    clines = gd_pnt2vec(clines,outype);
    
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
        [h_pnts(idx).UserData] = repmat(int32(0),sum(idx),1);
        cellobj = {[h_pnts(idx)]};
        [cellobj{:}.Color] = repmat(zeros(1,3),sum(idx),1);
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
        [h_lines(idx).UserData] = repmat(int32(0),sum(idx),1);
        cellobj = {[h_lines(idx)]};
        [cellobj{:}.Color] = repmat([1,0,0],sum(idx),1);
    end    
end

%%
function points = addLinePoints(ax)
    %define start and end of section line and plot it   
    prompt1 = 'Select start point';
    points(1) = gd_setpoint(ax,prompt1,false);
    prompt2 = 'Select end point';
    points(2) = gd_setpoint(ax,prompt2,false);
    plotLine(ax,points);
    nanpnts.x = NaN; nanpnts.y = NaN;        %line termination
    points = [points,nanpnts];   
end

%%
function lines = editpoint(ax,lines,edpoint,newpnt)
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
    [idl,idp] = findPointinLines(lines,edpoint);
    if idl<1, return;  end          %line not found
    lines{1,idl}(idp).x = newpnt.x; 
    lines{1,idl}(idp).y = newpnt.y;
    
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
        plotLine(ax,lines{1,idl},nline);
    end
end

%%
function lines = deleteline(ax,lines,deline)
    %delete lines and points for section that includes delpnt
    delpoint = deline(1);         %first point in line
    %find which line the point to be edited is on
    [idl,~] = findPointinLines(lines,delpoint);
    if idl<1, return;  end          %line not found
    deline = lines{1,idl};
    lines(idl) = [];                %delete from lines cell array

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
function [idl,idp] = findPointinLines(lines,qrypoint)
    %find which line and position of point on a set of lines (can be a
    %graphical array, or a cell array containing a struct array of points)
    nlines = length(lines);
    if  isgraphics(lines,'line') 
        glines = lines; lines = [];
        for i=1:nlines
            lp.x = glines(i).XData; 
            lp.y = glines(i).YData;
            lines{1,i} = gd_vec2pnt(lp); %#ok<AGROW> 
        end   
    end

    idl = -1; idp = [];
    for i=1:nlines
        isline = gd_findline(lines{1,i},qrypoint);
        if isline>0
            idl = i;
            idp = [lines{1,i}.x]==qrypoint.x & [lines{1,i}.y]==qrypoint.y;
            break
        end    
    end
end

%%
function   [slines,cpoints] = setInitialSections(ax,cpoints,islabel)
    %use centreline and boundary to initialise the sections lines
    clearLines(ax,{'mylines','mypoints','mytext'});
    clines = points2lines(cpoints);
    %get user to define the mouth by defining a point near to a point on the
    %channel centre-line
    promptxt = 'Left click to set mouth point, right click to quit';
    ok = 0;
    tol = (abs(diff(xlim))+abs(diff(ylim)))/200;
    while ok<1
        [mpnt,hp] = gd_setpoint(ax,promptxt,false);
        delete(hp)
        if isempty(mpnt)
            slines = []; return;   %figure has been deleted
        else
           [isNear,idL] = isPointNearLine([cpoints(:).x],[cpoints(:).y],mpnt.x,mpnt.y,tol); 
           if isNear, ok = 1; end
        end
        
    end
    
    %for each point from idL to the end use the centreline coordinates and
    %direction to define a section at right angles to the centreline
    [clinedir,clines,~] = clineProperties(clines,idL);

    %update plot of centre-line to extend from the selected mouth point
    clearLines(ax,{'clines'})
    cpoints = lines2points(clines);  %cpoints to match selected centre-line
    ax = plotCLine(ax,cpoints);

    %generate the section lines for clinedir +pi/2 and -pi/2
    maxlen = 1000;
    spoints = setSectionLines(clines,clinedir,maxlen);

    if islabel==1                                     %labelled sections
        slines = points2lines(spoints); 
        for j=1:length(slines)
            aline = (slines{1,j});
            ax = plotPoints(ax,aline,'mypoints');     %call one at a time
            sline = gd_setlines(ax,'',aline);         %to order numbering
            plotLine(ax,sline{1},num2str(j));     
        end        
    elseif islabel==2                                 %unlabled sections
        ax = plotPoints(ax,spoints,'mypoints');
        slines = gd_setlines(ax,'',spoints); 
    else
        slines = points2lines(spoints);               %dont plot
    end



%     hw = waitbutton(ax);
%     waitfor(hw,'Tag')
% 
% 
% 
%     slines = gd_pnt2vec(spoints,1);
%     P = InterX(slines',[bline.x';bline.y']);
% 
%     
%     %sections = points2lines(spoints);
%     newsections = cell(size(sections));
%     for i=1:length(sections)
%         for j=1:size(P,2)
%             ison = ispointonline(sections{i}(1:2),P(:,j));
%             if ison
%                 if size(newsections{i},2)==1
%                     newsections{i}(:,1) = P(:,j);
%                 else
%                     newsections{i}(:,2) = P(:,j);
%                     newsections{i}(:,3) = [NaN;NaN];
%                 end
% 
%             end
%         end
%     end
% 
%     hold on
%     plot(ax,P(1,:),P(2,:),'or')
%     hold off
%     spoints = [];
%     %spoints.x = []; spoints.y = [];
%     nlines = size(P,2);
%     for i=1:2:nlines-1
% %         spoints.x = [spoints.x,P(1,i),P(1,i+1),NaN];
% %         spoints.y = [spoints.y,P(2,i),P(2,i+1),NaN];
% 
%     end
% 
%     sections = points2lines(spoints);       
%     for j=1:length(sections)
%         ax = plotPoints(ax,sections{j}(1:2),'mypoints');
%         plotLine(ax,sections{j},num2str(j));
%     end
end
%%
function [clinedir,clines,cumlen] = clineProperties(clines,idL)
    %for each point from idL to the end use the centre-line coordinates and
    %direction to find the lengths anad directions along the centre-line
%     nlines = length(clines);
%     dx = cell(nlines,1);  dy = dx; slen = dx; cumlen = dx; clinedir = dx; 
%     delx = zeros(nlines,1); dely = delx;
    nrec = length(clines{1,1});                  %length of first line
    
    for i=1:length(clines)
        lp = clines{1,i};
        nl = length(lp);
        if idL>nrec                              %start point not in line
            nrec = nrec+nl;
            continue;
        elseif ~exist('dx','var')                  %start point in line
            dx = diff([lp(idL:end-1).x]);   %omit trailing NaN
            dy = diff([lp(idL:end-1).y]);   
            clines{1,i} = lp(idL:end);           %crop line to start point
        else                                     %subsequent lines
            dx = diff([lp(1:end-1).x]);     %omit trailing NaN
            dy = diff([lp(1:end-1).y]);  
        end
        dx = [dx(1),dx,dx(end)];               %pad to make same length as lines
        dy = [dy(1),dy,dy(end)];
        delx(i) = mean(abs(dx));             %mean length of dx
        dely(i) = mean(abs(dy));

        slen = hypot(dx,dy);           %length between points
        cumlen{i} = cumsum(slen);            %cumulative length
        theta = atan(dy./dx);             %direction between points
        clinedir{i}(1) = theta(1);
        for j=2:length(theta)                   %mean direction at point
            clinedir{i}(j) = (theta(j-1)+theta(j))/2;
        end    
    end
%     avdelx = mean(delx);
%     avdely = mean(dely);    
end

%%
% function [cumlen,clinedir,clines] = clineProperties(clines,idL)
%     %for each point from idL to the end use the centre-line coordinates and
%     %direction to find the lengths anad directions along the centre-line
%     nlines = length(clines);
%     dx = cell(nlines,1);  dy = dx; slen = dx; cumlen = dx; clinedir = dx; 
%     delx = zeros(nlines,1); dely = delx;
%     nrec = length(clines{1,1});                  %length of first line
%     
%     for i=1:length(clines)
%         lp = clines{1,i};
%         nl = length(lp);
%         if idL>nrec                              %start point not in line
%             nrec = nrec+nl;
%             continue;
%         elseif all([dx{:}]==0)                   %start point in line
%             [dx{i}] = diff([lp(idL:end-1).x]);   %omit trailing NaN
%             [dy{i}] = diff([lp(idL:end-1).y]);   
%             clines{1,i} = lp(idL:end);           %crop line to start point
%         else                                     %subsequent lines
%             [dx{i}] = diff([lp(1:end-1).x]);     %omit trailing NaN
%             [dy{i}] = diff([lp(1:end-1).y]);  
%         end
%         dx{i} = [dx{i}(1),dx{i}];               %pad to make same length as lines
%         dy{i} = [dy{i}(1),dy{i}];
%         delx(i) = mean(abs(dx{i}));             %mean length of dx
%         dely(i) = mean(abs(dy{i}));
% 
%         slen{i} = hypot(dx{i},dy{i});           %length between points
%         cumlen{i} = cumsum(slen{i});            %cumulative length
%         theta = atan(dy{i}./dx{i});             %direction between points
%         clinedir{i}(1) = theta(1);
%         for j=2:length(theta)                   %mean direction at point
%             clinedir{i}(j) = (theta(j-1)+theta(j))/2;
%         end        
%     end
% %     avdelx = mean(delx);
% %     avdely = mean(dely);    
% end

%%
function points = lines2points(lines)
    %convert cell array of lines to a cell array of points
    nrec = length(lines);
    points = [];
    for i=1:nrec
        points = [points,lines{1,i}]; %#ok<AGROW> 
    end
end

%%
function clines = points2lines(cpoints)
    %convert cell array of points to a cell array of lines
    idN = [0,find(isnan([cpoints(:).x]))];  
    for i=1:length(idN)-1
        idL =idN(i)+1:idN(i+1);    
        %cell so that number of points can vary
        clines{1,i} = cpoints(idL); %#ok<AGROW> 
    end
end

%%
function spoints = setSectionLines(clines,clinedir,maxlen)
    spoints = [];
    for i=1:length(clines)
        sldirpos = [clinedir{i}]+pi()/2;
        sldirneg = [clinedir{i}]-pi()/2;
        lpnts = clines{1,i};
        pospnts = sectionEndPoints(lpnts,sldirpos,maxlen);
        %plotPoints(ax,pospnts,'endpoints');
        negpnts = sectionEndPoints(lpnts,sldirneg,maxlen);
        %plotPoints(ax,negpnts,'endpoints');      
        %sections along ith length of centre-line
        nanpnts.x = NaN; nanpnts.y = NaN;        %line termination
        for j=1:length(pospnts)
            spoints = [spoints,pospnts(j),negpnts(j),nanpnts]; %#ok<AGROW> 
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
function [isNear,idL] = isPointNearLine(xData,yData,xPoint,yPoint,tol)
    % Calculate the distance from the point to the line
    distances = sqrt((xData - xPoint).^2 + (yData - yPoint).^2);
    [~,idL] = min(distances,[],'omitnan');
    isNear = any(distances < tol); % Threshold for proximity
end

%%
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
function plotLine(ax,line,nlinetxt)
    %plot a section line
    if nargin<3
        nlinetxt = num2str(length(findobj(ax,'Tag','mylines'))+1); 
    end
    hpts = findobj(ax,'Tag','mypoints');    
    hold on
    %plot(ax,[points(:).x],[points(:).y],'-r','PickableParts','none','Tag','mysection')
    hpts(1).MarkerSize = 7;                %make start point marker bigger
    text(ax,line(1).x,line(1).y,sprintf('%s',nlinetxt),...
        'HorizontalAlignment','center','PickableParts','none','FontSize',6,'Tag','mytext');
    hold off
end

%% Background graphics
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
    elseif isfield(cobj.Data,'GeoImage')
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
function ax = plotCLine(ax,points)
    %plot the clentreline as non-pickable linework
    hold on
    plot(ax,[points(:).x],[points(:).y],'+k','MarkerSize',4,...
                              'Tag','clines');
    plot(ax,[points(:).x],[points(:).y],'ok','MarkerSize',3,...
                              'Tag','clines');
    hold off
end

%%
function ax = toggle_view(ax,slines)
    %switch a line of the lines and points defined with the start point
    %of each line emphasised with a red circle marker.
    hline = findobj(ax,'Tag','mylines');
    if isempty(hline)           %toggle line and points on
        for j=1:length(slines)
            aline = (slines{1,j});
            ax = plotPoints(ax,aline,'mypoints');     %call one at a time
            sline = gd_setlines(ax,'',aline);         %to order numbering
            plotLine(ax,sline{1},num2str(j));     
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



