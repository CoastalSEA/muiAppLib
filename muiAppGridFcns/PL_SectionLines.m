classdef PL_SectionLines < PLinterface
%
%-------class help---------------------------------------------------------
% NAME
%   PL_SectionLines.m
% PURPOSE
%   Class to create a set of cross-sections lines based on the points in a
%   centre-line and any enclosing boundary
% USAGE
%  
% SEE ALSO
%   inherits PLinterface
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2025
%--------------------------------------------------------------------------
%    
    properties (Transient)
        %set Abstract properties for PLinterface
        outPoints = [] %not used but has to be set        
        isXYZ = false  %not used but has to be set
        outLines = []
        Set            %instance of PL_Sections class, or struct with 
                       %same field names
        cLines         %updated plines version of Set.ChannelLine (modified
                       %when sections are set but not saved)
       xLength         %half length of cross-sections
    end
       
    methods
        function obj = PL_SectionLines(figtitle,tag,position)          
            %constructor code: 
            if nargin<3, position = []; end
            obj = setFigure(obj,figtitle,tag,position);
        end 
    end
%% 
    methods (Static)  
        function [slines,clines] = Figure(grid,promptxt,setlines,isdel)
            %
            %-------function help------------------------------------------
            % PURPOSE
            %   Figure to interactively edit points or lines
            % USAGE
            %   lines = PL_SectionLines.Figure(grid,promptxt,setlines,isdel);
            % INPUTS
            %   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
            %   promptxt- character string used for initial prompt in title
            %   setlines - struct with fields Boundary, ChannelLine, 
            %              ChannelProps and SectionLines, or an instance 
            %              of the PL_Sections class
            %   isdel - logical flag true to delete figure on completion - optional, 
            %           default is false
            % OUTPUTS
            %   lines - if lines are input, format is the same as the input, otherwise
            %           if nlines is scalar (default it type 2 if nlines empty)
            %            outype=0: array of structs with x, y and z fields defining selected points,
            %            outype=1: Nx2 or Nx3 array.
            %            outype=2: struct with x, y (and z) vector fields
            %            outype=3: table with x, y (and z) vector fields
            %            lines = [] if user closes figure, or no points defined
            % NOTES
            %   NB: if the figure window is closed the function returns lines=[] and
            %       not the input points
            %--------------------------------------------------------------
            %
            if nargin<4, isdel = false; end
            figtitle = sprintf('Extract centre-line');
            tag = 'PlotFig'; %used for collective deletes of a group
            position = [0,0.03,1,0.93];
            obj = PL_SectionLines(figtitle,tag,position);

            %plot grid and initialise axes (needed for context menus)
            obj.Axes = plotGrid(obj,grid);            
            obj.Axes.Title.String = promptxt; %initial prompt

            %define menu to be used
            mtext = {'Sections','Figure'};
            mcall = repmat({[]},1,length(mtext));
            menu = struct('label',mtext,'callback',mcall);
            %default menus are defined for Point, Line and Figure include:
            % Line: 'Add','Edit','Extend','Insert','Join','Split','Delete'
            % Figure: 'View','Undo','Save','Save & Exit','Quit'
            %and are defined in PLinterface. getCallBacks
            %Bespoke options are then added by creating tables for each
            %addtional menu option in the function setSubMenus        
            submenus = setSubMenus(obj);
            %set the menus and submenus            
            obj = setMenu(obj,menu,submenus);
            
            obj.Set = setlines;                  %PL_Sections instance
            b_lines = obj.Set.Boundary;          %Boundary coordinates   
            hold on
            plot(obj.Axes,b_lines.x,b_lines.y,'k');
            hold off;
     
            obj.cLines = gd_lines2points(obj.Set.ChannelLine);%Centre line coordinates
            gd_plotpoints(obj.Axes,obj.cLines,'clines',5);    %5= plot as centre-lines

            %handle input of existing lines
            outype = getInLines(obj,obj.Set.SectionLines);

            obj.Figure.Visible = 'on';
            obj = waitForFigure(obj);

            if isempty(obj.outLines)
                slines = []; clines = [];
            else                
                slines = gd_points2lines(obj.outLines,outype.lines);
                clines = gd_points2lines(obj.cLines,outype.lines);
            end

            %delete figure if isdel has been set by call.
            if isdel
                delete(obj.Figure)
            end
        end
    end

%--------------------------------------------------------------------------
% Additional Menu callback functions
%--------------------------------------------------------------------------
    methods (Access=protected)
        function setSections(obj,~,~)
            %use centreline to generate a set of sections    
            resetMenu(obj);
            clearGraphics(obj,{'mylines','mypoints','mytext'});
            setSectionLength(obj);

            %get user to define the mouth by defining a point near to a point on the
            %channel centre-line
            promptxt = 'Left click to set mouth point, right click to quit';
            ok = 0;
            axes(obj.Axes)
            tol = (abs(diff(xlim))+abs(diff(ylim)))/2/100; %(1e4+1e5)/2/100~O[300m]            
            while ok<1
                [mpnt,hp] = gd_setpoint(obj.Axes,promptxt,'startpoint',false);
                delete(hp)
                if isempty(mpnt)
                    return;                        %figure has been deleted
                else
                   [isNear,idL] = PLinterface.isPointNearLine(obj.cLines,mpnt,tol);                        
                   if isNear, ok = 1; end
                end
            end     
            %for each point from idL to the end use the centreline coordinates and
            %direction to define a section at right angles to the centreline
            c_cplines = gd_plines2cplines(obj.cLines);
            [clinedir,c_cplines,~] = gd_curvelineprops(c_cplines,idL(1));
        
            %generate the section lines for clinedir +pi/2 and -pi/2
            obj.pLines = setSectionLines(obj,c_cplines,clinedir,obj.xLength);
            obj.cLines = gd_cplines2plines(c_cplines);
            plotSections(obj);
            resetMenu(obj,false)  
        end

%%
        function clipLines(obj,~,~)
            %clip the extent of the sections lines to a defined boundary
            resetMenu(obj);
            c_lines = obj.cLines;
            c_lines(isnan([c_lines(:).x])) = [];   %remove Nans so only points
            ncentres = length(c_lines);            %number of centreline points 

            s_cplines = gd_plines2cplines(obj.pLines);  %section lines
            nsections = length(s_cplines);
            if nsections~=ncentres
                warndlg('Number of sections do not match number of points on centre-lines')
                return;
            else
                c_lines = gd_points2lines(c_lines,1);
            end
            %crop lines to first crossing of boundary line
            [s_lines,idc] = clipSectionLines(obj,s_cplines,c_lines);
            obj.pLines = gd_lines2points(s_lines);

            %idc is the centre-line index of deleted points excluding NaNs

            for i=1:length(idc)
                 idd = find([obj.cLines(:).x]==c_lines(idc(i),1) & ...
                                    [obj.cLines(:).y]==c_lines(idc(i),2)); 
                 if ~isempty(idd), obj.cLines(idd) = []; end
            end
            clearGraphics(obj,{'mylines','mypoints','mytext','clines'});
            ax = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);  %set line  
            obj.Axes = gd_plotpoints(ax,obj.cLines,'clines',5); %set centre-lines
            resetMenu(obj,false)  
        end

%%
        function addLine(obj,~,~)
            %order the lines to be in along-channel sequence
            resetMenu(obj);
            % prompt1 = sprintf('Add line\nSelect first point of section');   
            % prompt2 = sprintf('Add line\nSelect second point of section');
            % prompts = {prompt1,prompt2};
            promptcl= sprintf('Add line\nLeft click to select centre-line, right click on a line to quit');
            promptcp = sprintf('Add line\nLeft click to select centre-line point, right click to quit');

            clearGraphics(obj,{'mylines','mypoints','mytext','clines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.cLines,'clines',2);      %set centre-lines
            cpline = gd_getpline(obj.Axes,promptcl,'clines');              %get line to add to
            while~isempty(cpline)
                %gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);            %show sections
                gd_plotpoints(obj.Axes,cpline,'clines',5);                 %show centre-line points
                [NewPnt,Hc] = gd_setpoint(obj.Axes,promptcp,'clpoint',obj.isXYZ); 
                while ~isempty(NewPnt)  
                    cpline = insertSection(obj,NewPnt,cpline);             %insert new section
                    delete(Hc);                                            %delete digitising points   
                    gd_plotpoints(obj.Axes,NewPnt,'clines',5);             %set centre-line point
                    [NewPnt,Hc] = gd_setpoint(obj.Axes,promptcp,'clpoint',obj.isXYZ);
                end
                clearGraphics(obj,{'mylines','clines'});
                gd_plotpoints(obj.Axes,obj.cLines,'clines',2);             %set centre-lines
                cpline = gd_getpline(obj.Axes,promptcl,'clines');          %get line to add to
            end

            clearGraphics(obj,{'mylines','clines'});
            obj.Axes = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);%set line  
            obj.Axes = gd_plotpoints(obj.Axes,obj.cLines,'clines',5); %set centre-lines
            resetMenu(obj,false)  
        end

%%
        function deleteLine(obj,~,~)
            %order the lines to be in along-channel sequence
            resetMenu(obj);
            promptxt = sprintf('Delete line\nSelect line to Delete, right click on any line to quit\nRedraw to update centre-line');
            [deline,H] = gd_getpline(obj.Axes,promptxt,'mylines');         %get line to delete
            while ~isempty(deline)                
                [obj.pLines,idl] = deleteAline(obj,'pLines',deline);       %delete the line
                delete(H)
                obj.cLines(idl) = [];                                      %delete centre-line point
                clearGraphics(obj,{'mylines','clines'});
                ax = gd_plotpoints(obj.Axes,obj.pLines,'mylines',2);       %set line 
                gd_plotpoints(ax,obj.cLines,'clines',5);                   %set centreline
                deline = gd_getpline(obj.Axes,promptxt,'mylines');         %get line to delete
            end
            resetLines(obj)
            resetMenu(obj,false)       
        end

%%
        function Reset(obj,~,~)
            %Reset to the input centre-line and boundary
            clearGraphics(obj,{'mylines','mypoints','mytext','clines'});
            obj.cLines = gd_lines2points(obj.Set.ChannelLine);%Centre line coordinates
            gd_plotpoints(obj.Axes,obj.cLines,'clines',5);    %5= plot as centre-lines            
        end


%%
        function Redraw(obj,~,~)
            %redraw the current set of points and lines
            ax = obj.Axes;            
            clearGraphics(obj,{'mylines','clines'});

            if ~isempty(obj.pLines)   
                ax = gd_plotpoints(ax,obj.pLines,'mylines',2); %set line          
            end
  
            if ~isempty(obj.cLines) 
                gd_plotpoints(ax,obj.cLines,'clines',5);       %set centreline
            end
        end  

%%
        function Label(obj,~,~)
            %toggline text section markers on and off
            hpnts = findobj(obj.Axes,'Tag','mytext');
            if isempty(hpnts)           %toggle points on
                cplines = gd_plines2cplines(obj.pLines);
                for i=1:length(cplines)
                    gd_plotpoints(obj.Axes,cplines{i},num2str(i),3);  %set points
                end
            else
                clearGraphics(obj,{'mytext'});
            end   
        end
%%
%--------------------------------------------------------------------------
%  Utility functions to implement various actions
%-------------------------------------------------------------------------- 
        function plines = setSectionLines(~,cplines,clinedir,maxlen)
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
        function [s_lines,idc] = clipSectionLines(obj,s_cplines,c_lines)
            %find the nearest boundary point crossing positive and negative lines
            b_lines = obj.Set.Boundary;
            s_lines = [];  spos = []; sneg = [];  idc = []; %initialise loop variables            
            for i=1:length(s_cplines)
                s_line = gd_points2lines(s_cplines{1,i},1);
                %InterX input is 2 row vectors for x and y
                P = InterX(s_line',[b_lines.x';b_lines.y']); %intersections for 1 line
                if ~isempty(P)
                    cpnt = c_lines(i,:)';
                    %find first intersections between centre-line and initial 
                    %section end points
                    % hold on
                    % plot(obj.Axes,cpnt(1),cpnt(2),'ow','MarkerSize',8)
                    % hold off                     
                    k = 1; l = 1;
                    for j=1:size(P,2)  
                        %for each intersection point check whether point is on
                        %positive or negative side of line
                        posline = [cpnt,s_line(1,:)'];   
                        isp = ispointonline(posline,P(:,j),true,1e3); %true=check ends
                        if isp
                            spos(:,k) = P(:,j); %#ok<AGROW> 
                            k = k+1;
                        end
        
                        negline = [cpnt,s_line(2,:)'];
                        isn = ispointonline(negline,P(:,j),true,1e3); %true=check ends
                        if isn
                            sneg(:,l) = P(:,j); %#ok<AGROW> 
                            l = l+1;
                        end
                    end
                    %
                    if isp || isn
                        spos = findNearest(spos,cpnt);
                        sneg = findNearest(sneg,cpnt);
                        if ~isempty(spos) &&  ~isempty(sneg)
                            s_lines = [s_lines;spos';sneg';[NaN,NaN]];     %#ok<AGROW> 
                        else
                            idc = [idc,i]; %#ok<AGROW>
                        end
                    else
                        %no points found for centreline point - remove point
                        idc = [idc,i]; %#ok<AGROW>
                    end                    
                    % hold on
                    % pp = [spos,sneg];
                    % plot(obj.Axes,pp(1,:),pp(2,:),'ok')
                    % hold off   
                    spos = []; sneg = [];   %clear for next iteration
                else
                    idc = [idc,i]; %#ok<AGROW>
                end
            end

            %nested function-----------------------------------------------
            function point = findNearest(points,cpnt)
                if size(points,2)>1
                    p = abs(points-cpnt);
                    D = hypot(p(1,:),p(2,:));
                    [~,idx] = min(D);
                    point = points(:,idx);
                else
                    point = points;
                end
            end %----------------------------------------------------------
        end

%%
        function cline = insertSection(obj,newpnt,cline)
            %insert centre-line point on obj.cLines and cline
            tol = obj.lineLength(cline(1),cline(2));
            [ison,idP] = obj.isPointNearLine(cline,newpnt,tol);
            if ison
                %insert the new centre-line point
                [obj.cLines,~] = insertPoints(obj,'cLines',cline(idP([1,2])),newpnt);
            else
                %point not found in line so assume it is extending line
                
                if idP(1)==1 || idP(1)==length(cline)-1  %starat or end of line
                    obj.cLines = extendAline(obj,'cLines',newpnt,cline(idP(1)));
                else
                    getdialog(sprintf('Centre-line point not found\nSection not added to dataset'))
                    return;
                end                
            end
            idL = gd_findline(obj.cLines,cline(1));
            cplines = gd_plines2cplines(obj.cLines);
            cline = cplines{1,idL};
            %for each point use the centreline coordinates and
            %direction to define a section at right angles to the centreline
            if isempty(obj.xLength), setSectionLength(obj); end
            %set the length of the section lines 
            c_cplines = gd_plines2cplines(obj.cLines);
            [clinedir,c_cplines,~] = gd_curvelineprops(c_cplines,1);
            obj.pLines = setSectionLines(obj,c_cplines,clinedir,obj.xLength);
        end

%%
        function setSectionLength(obj)
            %set the length of the section lines
            inp = inputdlg({'Length of section lines from centre-line'},'Sections',...
                            1,{'1000'});
            if isempty(inp), inp{1} = '1000';  end
            obj.xLength = str2double(inp{1}); 
        end

%%
        function plotSections(obj)
            %plot the derived sections as editable lines
            ax = obj.Axes;
            cplines = gd_plines2cplines(obj.pLines);
            for j=1:length(cplines)                      %call one at a time
                aline = (cplines{1,j});                  %to order numbering
                ax = gd_plotpoints(ax,aline,'mylines',2);        %set line
                obj.Axes = gd_plotpoints(ax,aline,num2str(j),3); %set labels
            end
        end

%%
        function ax = plotGrid(obj,grid)
            %plot either the bathymetry grid or an image of it as a backdrop
            isgrid = false; isimage = false;
            hfig = obj.Figure;
            if isstruct(grid)
                %xyz grid is passed as the input variable
                ax = gd_plotgrid(hfig,grid);
                hplt = findobj(ax,'Tag','PlotGrid');
                hplt.Annotation.LegendInformation.IconDisplayStyle = 'off';  
                isgrid = true;
            else
                ax = axes(hfig);
                %image is passed to class as the 'grid' input variable
                img = grid.geoimage;     %image object
                h_im = imagesc(ax,'XData',img.XData,'YData',img.YData,'CData',img.CData);
                set(h_im, 'AlphaData', 1-isnan(img.CData)); %set Nan values to be transparent   
                colormap(img.CMap);
                clim(img.CLim);
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
        function calltable = setSubMenus(obj)  
            %user defined menus to be appended to the Figure menu 
            varnames = {'Parent','Callback','Label'};
            %Boundary menu variables
            stext = ["Set";"Clip";"Add";"Edit";"Delete";"Reset";"Label"];
            scall = {@obj.setSections; @obj.clipLines; @obj.addLine; ...
                     @obj.editLine; @obj.deleteLine; @obj.Reset; @obj.Label}; 
            nrec = length(stext);
            spart = repmat("Sections",nrec,1);
            calltable = table(spart,scall,stext,'VariableNames',varnames);
        end
    end
%--------------------------------------------------------------------------
% Static utility functions
%--------------------------------------------------------------------------
    methods (Static, Access=protected)
        %Static methods in PLinterface
        % checkDirection
        % lineLength
        % isPointNearLine
        % setLevel
        % setInterval
    end
end