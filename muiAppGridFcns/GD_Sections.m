classdef GD_Sections < handle   
%
%-------class help---------------------------------------------------------
% NAME
%   GD_Sections.m
% PURPOSE
%   Class to extract sections from a grid
% USAGE
%   obj = GD_Sections.sectionsMenu(mobj,src,gridclasses);
%         %mobj is a handle to Main UI
%         %src is handle to calling menu option
%         %gridclasses is cell array of classes to use in case selection
% SEE ALSO
%   inherits handle
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%      
    properties
        Boundary            %lines used to clip cross-section lines
        ChannelLine         %lines used to define channel network
        ChannelProps        %properties used to extract channel network
                            %fields: maxwl,dexp,cint
        SectionLines        %lines that define cross-sections
        XSections           %cross-sections obtained from grid
    end
       
    methods (Access=protected)
        function obj = GD_Sections()          
            %constructor code:            
        end 
    end
%%  
    methods (Static)  
        function sectionMenuOptions(mobj,src,gridclasses)
            %default set of menu options for use in Model UIs
            % mobj - mui model instancegetGrid
            % src - handle to calling menu option
            % gridclasses - cell array list of classes to use in case selection
            %
            % Sub-menu for using Section Tools:
            % menu.Setup(X).List = {'Boundary',...
            %                       'Channel network',...
            %                       'Section lines',...
            %                       'Sections'};                                                                        
            % menu.Setup(X).Callback = repmat({@obj.sectionsMenu},[1,4]);
            % menu.Setup(X).Separator = repmat({'off'},[1,4]);  
            %
            % Additional submenus may been be added for additional options:
            % menu.Setup(9).List = {'Generate','Load','Edit'};
            % menu.Setup(9).Callback = repmat({@obj.sectionMenuOptions},[1,3]);
            % where Load reads a shapefile, Generate uses the GD_Sections
            % methods and Edit calls gd_editlines (see EstauryDB for eg)
            %
            muicat = mobj.Cases;   %handle to muiCatalogue
            switch src.Text 
                case 'Boundary'
                    promptxt = 'Select a Case to use to define boundary:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = GD_Sections.getSections(cobj);
                    setBoundary(obj,cobj,muicat); 
                case 'Channel network'
                    promptxt = 'Select a Case to use to define channel network:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = GD_Sections.getSections(cobj);
                    setChannelNetwork(obj,cobj,mobj); 
                case 'Section lines'
                    promptxt = 'Select a Case to use to define section lines:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = GD_Sections.getSections(cobj);
                    setSectionLines(obj,cobj,muicat);
                case 'Sections'
                    promptxt = 'Select a Case to use to extract sections:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = GD_Sections.getSections(cobj);
                    setSections(obj,cobj,muicat);
                 case 'View Sections'
                    promptxt = 'Select a Case to use to extract sections:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj.Sections)
                        warndlg('No Section data available to view');
                        return;
                    end  
                    obj = cobj.Sections;
                    viewSections(obj,cobj,muicat);
            end
        end

%%
        function obj = getSections(cobj)
            %retrieve existing class object or create a new instance
            if nargin<1, cobj = []; end

            if isa(cobj.Sections,'GD_Sections')
                obj = cobj.Sections;
            else                
                obj = GD_Sections;
            end
        end

%%
        function loadLines(mobj,src,gridclasses)
            %load linework for selected line type 
            %can be Boundary, Channel Network or Section Lines
            muicat = mobj.Cases;   %handle to muiCatalogue
            linetype = src.Text;
            switch linetype
                case 'Boundary'
                    type = 'Boundary';
                case 'Channel network'
                    type = 'ChannelLine';
                case 'Section lines'
                    type = 'SectionLines';
            end
            promptxt = sprintf('Select a Case to load %s from shapefile:',linetype);           
            [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);

            [fname,path,nfiles] = getfiles('MultiSelect','off',...
                'FileType',{'*.shp;'},'PromptText','Select shape file):');
            if nfiles~=1, return; end

            %get file and read the data
            shp = gd_readshapefile(path,fname);
            if isrow(shp.x)
                shp = structfun(@transpose,shp,'UniformOutput',false);
            end
            cobj.Sections.(type) = shp;  %overwrites any existing lines
            casedesc = muicat.Catalogue.CaseDescription(cobj.CaseIndex);
            getdialog(sprintf('Data loaded for %s in Case: %s',linetype,casedesc),[],3)
        end

%%
        function editLines(mobj,src,gridclasses)
            %edit linework for selected line type from a shapefile
            %can be Boundary, Channel Network or Section Lines
            muicat = mobj.Cases;   %handle to muiCatalogue
            linetype = src.Parent.Text;
            switch linetype
                case 'Boundary'
                    type = 'Boundary';
                case 'Channel network'
                    type = 'ChannelLine';
                case 'Section lines'
                    type = 'SectionLines';
                otherwise
                    return
            end
            
            promptxt = sprintf('Select a Case to edit %s:',linetype);
            [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);

            if ~isa(cobj.Sections,'GD_Sections')  && strcmp(src.Text,'Edit')
                warndlg(sprintf('No %s data available to edit',linetype));
                return;
            elseif strcmp(src.Text,'Digitise')
                clines = 2;  %2 define outype as struct array of xy
            else
                clines = cobj.Sections.(type);
            end

            %eidt the selected line and save to case
            
            paneltxt = sprintf('Edit the %s',linetype);
            grid = getGrid(cobj,1);   %grid for selected year
            clines = gd_editlines(grid,paneltxt,clines,true);
            if isempty(clines), return; end
            answer = questdlg(sprintf('Save the edited %s',linetype),linetype,'Yes','No','Yes');
            if strcmp(answer,'Yes')
                cobj.Sections.(type) = clines;
                getdialog(sprintf('Edits saved for %s in Case: %s',linetype,grid.desc),[],3)
            end
        end
    end
%%
    methods
        function setBoundary(obj,cobj,~)
            %load a grid and interactively define a boundary line with
            %options to autogenerate and edit or digistise the line.

            if ~isfield(cobj.Data,'Grid')
                warndlg('No grid for selected case'); return;
            end
            grid = getGrid(cobj,1);             %grid for estuary

            blines = 2;                         %output format when new line
            if ~isempty(obj.Boundary)           %channel line exists
                answer = questdlg('A centre-line exists. Extend existing or create new one?',...
                                  'Channel','Extend','New','Extend');
                if strcmp(answer,'Extend')
                    blines = obj.Boundary;      %existing lines
                end                    
            end 

            %clean-up shoreline
            paneltxt = 'Create and smooth the shoreline boundary';
            shore_xy = gd_setboundary(grid,paneltxt,blines,true);
            if isempty(shore_xy), return; end

            %give user opportunity to edit the shoreline interactively
            answer = questdlg('Check/Edit or Save the shoreline?','Centre-line','Edit','Save','Quit','Edit');
            if strcmp(answer,'Edit')
                paneltxt = 'Edit the shoreline boundary';
                shore_xy = gd_editlines(grid,paneltxt,shore_xy,true);
                if isempty(shore_xy), return; end
            elseif strcmp(answer,'Quit')
                return;
            end

            %save results to grid class instance
            obj.Boundary = shore_xy;
            cobj.Sections = obj;
            getdialog(sprintf('Boundary added to grid for Case: %s',grid.desc),[],3)
        end
%%
        function setChannelNetwork(obj,cobj,mobj)
            %load a grid and extract the centre lines of the channels in
            %the network, or digitise them manually, and combine into a
            %vector of points
            if ~isfield(cobj.Data,'Grid')
                warndlg('No grid for selected case'); return;
            end
            grid = getGrid(cobj,1);   %grid for selected year
            
            if ~isempty(obj.ChannelLine)           %channel line exists
                answer = questdlg('A centre-line exists. Extend existing or create new one?',...
                                  'Channel','Extend','New','Extend');
            else
                answer = 'New';
            end   

            if strcmp(answer,'New')
            %get maximum water level to define limits of search
                promptxt = {'Maximum accessible water level?','Depth exponent','Sampling interval (m)'};
                defaults = {num2str(max(grid.z,[],'all')),'5','100'};
                inp = inputdlg(promptxt,'Water level',1,defaults);
                if isempty(inp), return; end %user cancelled
                props.maxwl = str2double(inp{1});
                props.dexp = str2double(inp{2});
                props.cint = str2double(inp{3});
                nlines = [];
            else
                props = obj.ChannelProps;
                c_plines = obj.ChannelLine;
                nlines(:,1) = c_plines.x; nlines(:,2) = c_plines.y;  %matrix of xy points
            end

            nlines = gd_setcentrelines(grid,mobj,props,nlines,true);

            ok = 0;
            while ok<1
                %loop adding lines until user cancels without defining
                %a new line
                [nline, hf] = setSingleLine(grid,mobj,props,nlines);
                if isempty(nline), ok =1; continue; end
                nlines = setCentreLine(props,nline,nlines);
                delete(hf)
            end
            if isempty(nlines), return; end                    %user cancelled without creating any lines

            %give user opportunity to edit the centre-line interactively
            promptxt = sprintf('Check/Edit or Save the centre-lines?\nTake opportunity to clean channel connections');
            answer = questdlg(promptxt,'Centre-line','Edit','Save','Quit','Edit');
            if strcmp(answer,'Edit')
                paneltxt = 'Edit the centre-lines and location of points at nodes';
                c_plines = gd_editlines(grid,paneltxt,c_plines,true);
                if isempty(c_plines), return; end
            elseif strcmp(answer,'Quit')
                return;
            end

            %define connectivity of finished lines
            c_plines.x = nlines(:,1)'; c_plines.y =  nlines(:,2)'; %x,y struct of vector points (row vectors)
            %props.topo = gd_linetopology(grid,c_plines);

            %save results to grid class instance
            obj.ChannelLine = c_plines;
            obj.ChannelProps = props;
            cobj.Sections = obj;
            getdialog(sprintf('Channel network added to grid for Case: %s',grid.desc),[],3)
            %option to also add to grid cline property
            % cobj.Data.Grid.UserData.cline.x = clines(:,1);
            % cobj.Data.Grid.UserData.cline.y = clines(:,2);

            %nested function-----------------------------------------------
            function [nline,hf] = setSingleLine(grid,mobj,props,nlines)
                %
                nline = gd_centreline(grid,mobj,props,nlines);
                if isempty(nline), hf = []; return; end
                nline = cell2mat(struct2cell(nline)');        %struct to [Nx2] array
                hf = GD_Sections.cl_checkPlot(grid,nline,nlines);
            end

            %nested function-----------------------------------------------
            function nlines = setCentreLine(props,nline,nlines)
                %
                nline = nline(1:end-1,:);                     %remove trailing NaNs
                ansr = questdlg('Accept the centreline?','Centre-line','Yes','No','Yes');
                if strcmp(ansr,'Yes')
                    %convert format of output if required
                    clength = sum(vecnorm(diff(nline),2,2));  %nline is a column vector [Nx2]
                    cpoints = round(clength/props.cint);      %number of points in new line
                    newcline = curvspace(nline,cpoints);      %curvespace uses [Nx2]
                    nlines = [nlines;newcline;[NaN,NaN]];     %restore trailing NaNs                    
                end
                %plines = struct('x',nlines(:,1),'y',nlines(:,2));
            end %----------------------------------------------------------
        end

%%
        function setSectionLines(obj,cobj,muicat)
            %use the centre-line to generate a set of sections at right
            %angles and interactively edit them, or digitise them manually.
            %at the moment this requires a grid to be present but not needed ** 
            if isempty(obj.Boundary) || isempty(obj.ChannelLine)
                warndlg('Boundary and Channel network need to be defined to extract Section Lines')
                return
            end
            %extract the section lines that are normal to the 
            %channel centre line and extend to the bounding shoreline
            paneltxt = 'Extract the channel cross-sections';
            [slines,clines,cumlen] = gd_setsectionlines(obj,cobj,paneltxt,false);
            if isempty(slines), return; end

            %give user opportunity to edit the centre-line interactively
            grid = getGrid(cobj,1);   %grid for selected year
            answer = questdlg('Check/Edit or Save the Cross-sections?','Cross-sections','Edit','Save','Quit','Edit');
            if strcmp(answer,'Edit')
                paneltxt = 'Edit or Reorder the cross-sections';                
                slines = gd_editlines(grid,paneltxt,slines,false);
                if isempty(slines), return; end
            elseif strcmp(answer,'Quit')
                return;
            end

            obj.SectionLines = slines;
            obj.ChannelLine = clines;
            obj.ChannelProps.ChannelLengths = cumlen;
            cobj.Sections = obj;
            casedesc = muicat.Catalogue.CaseDescription(cobj.CaseIndex);
            msgtxt = sprintf('Section lines added to grid for Case: %s',...
                                                                casedesc);
            getdialog(msgtxt,[],3)
        end

%%
        function setSections(obj,cobj,muicat)
            %using the SectionLines data extract the widths as a function 
            %of elevation and distance along the centre-line from the mouth
            if isempty(obj.SectionLines) || ~isfield(cobj.Data,'Grid')
                getdialog('No data',[],1); return; 
            end
            grid = getGrid(cobj,1);
            
            promptxt = {'Maximum accessible water level (mAD)',...
                        'Sampling interval (m)','Interpolation method'};
            defaults = {'10','10','makima'};
            inp = inputdlg(promptxt,'Water level',1,defaults);
            if isempty(inp), return; end %user cancelled
            zmax = str2double(inp{1});
            sint = str2double(inp{2});
            method = inp{3};

            cplines = gd_plines2cplines(gd_lines2points(obj.SectionLines));
            %plot the elevations for the defined section
            nlines = length(cplines);            
            for i=1:nlines
                points = cplines{1,i};
                xlen = diff([points(1:end-1).x]);
                ylen = diff([points(1:end-1).y]);
                slen = hypot(xlen,ylen);            %length of section        
                spnts = 0:sint:slen;                %points along section
                xq = points(1).x+spnts/slen*xlen;
                yq = points(1).y+spnts/slen*ylen;  
        
                %check for NaNs in grid
                grid.z(isnan(grid.z)) = zmax;      %clean grid for interpolation
                grid.z(grid.z>zmax) = zmax;
                [X,Y] = meshgrid(grid.x,grid.y);
                zline = interp2(X,Y,grid.z',xq,yq,method);
                if zline(1)<zmax          %adjust start point if below zmax
                    zline = [zmax,zline]; %#ok<AGROW> 
                    spnts = [0,spnts+sint];
                end
                if zline(end)<zmax        %adjust end point if below zmax
                    zline = [zline,zmax]; %#ok<AGROW> 
                    spnts = [0,spnts+sint];
                end
                
                zplines(1,i) = struct('x',[spnts,NaN],'y',[zline,NaN]); %#ok<AGROW> 
            end            
            obj.XSections = gd_points2lines(zplines,2);
            cobj.Sections = obj;
            casedesc = muicat.Catalogue.CaseDescription(cobj.CaseIndex);
            viewAlongChannelSections(obj,casedesc);
        end

%%
        function viewSections(obj,cobj,muicat)
            %view boundary channel network and cross-sections line work
            %or along-channel sections summary plot
            answer = questdlg('Plan view or Along-channel sections?','sections',...
                              'Plan','Along-channel','Plan');
            if strcmp(answer,'Plan')
                viewPlanSections(obj,cobj);
            else
                casedesc = muicat.Catalogue.CaseDescription(cobj.CaseIndex);
                viewAlongChannelSections(obj,casedesc);
            end
        end
%%
        function viewPlanSections(obj,cobj)
            %plot boundary channel network and cross-sections line work
            hf = figure('Name','Sections','Units','normalized',...
                                        'Tag','PlotFig','Visible','off');  
            isgrid = false; isimage = false;
            if isfield(cobj.Data,'Grid')
                dst = cobj.Data.Grid;
                grid = getGrid(cobj,1);                 %grid selected
                ax = gd_plotgrid(hf,grid);
                hplt = findobj(ax,'Tag','PlotGrid');
                hplt.Annotation.LegendInformation.IconDisplayStyle = 'off';  
                isgrid = true;
            elseif isfield(cobj.Data,'GeoImage')
                dst = cobj.Data.GeoImage;
                im = dst.geoimage;                      %image object
                ax = axes(hf);
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

            type = {'Boundary','ChannelLine','SectionLines'};
            for i=1:3
                switch type{i}
                    case 'Boundary'
                        sc = 'k'; ss = '-'; sw = 0.5;
                    case 'ChannelLine'
                        sc = 'g'; ss = '-.'; sw = 1.5;
                    case  'SectionLines'
                        sc = 'r'; ss = '-'; sw = 1;
                end
                %
                if ~isempty(obj.(type{i}))
                    lines = obj.(type{i});
                    hold on
                    plot(ax,lines.x,lines.y,'Color',sc,'LineStyle',ss,...
                                     'LineWidth',sw,'DisplayName',type{i});
                    if strcmp(type{i},'ChannelLine')
                        hp = plot(ax,lines.x(1),lines.y(1),'ok','MarkerSize',4,...
                                                        'MarkerFaceColor','w');
                        hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
                    hold off
                end
            end

            %axis equal tight
            xlabel('Eastings (m)'); 
            ylabel('Northings (m)');    
            title(dst.Description);    
            legend
            hf.Visible = 'on';
        end

%%
        function viewAlongChannelSections(obj,casedesc)
            %plot along-channel sections
            hf = figure('Name','Sections','Units','normalized',...
                             'Tag','PlotFig','Visible','on');
            ax = axes(hf);
            glines = {'-','--',':','-.'};
            hold on
            plines = obj.XSections;
            cplines = gd_plines2cplines(gd_lines2points(plines));
            nlines = length(cplines);
            for i=1:nlines
                aline = cplines{1,i};
                spnts = [aline(:).x];
                zpnts =[ aline(:).y];
                nline = length(findobj(ax,'Tag','asection'));
                lname = sprintf('Section %d',nline+1);        
                plot(ax,spnts,zpnts,'LineStyle',glines{rem(nline,4)+1},...
                      'LineWidth',1,'Tag','asection','DisplayName',lname,...
                      'ButtonDownFcn',@godisplay)
            end
            hold off
            if nlines<20
                legend
            end
            title(sprintf('Along-channel sections for %s',casedesc))    
        end
    end

%%
%     methods (Static, Access=private)
%         function hf = cl_checkPlot(grid,cline,nlines)
%             %plot base map of initial grid selection and defined mask
%             % cline and nlines are both [Nx2] arrays
%             points = [cline(1,:);cline(end,:)];
%             hf = figure('Name','Thalwegs','Units','normalized','Tag','PlotFig');  
%             hf.Position = [0,0,1,1];
%             ax = gd_plotgrid(hf,grid);
%             colormap(ax,'gray');
%             lines = {'-','--',':','-.'};
%             
%             hs = findobj(ax.Children,'Type','surface');
%             hs.Annotation.LegendInformation.IconDisplayStyle = 'off';
%             hold on
%             hp = plot(ax,points(:,1),points(:,2),'ok','MarkerSize',8,...
%                                    'MarkerFaceColor','w','Tag','mypoints');
%             hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
%             plot(ax,cline(:,1),cline(:,2),'r','LineStyle',lines{1},'LineWidth',2,...
%                                            'DisplayName','New line');
%             if ~isempty(nlines)
%                 plot(ax,nlines(:,1),nlines(:,2),'g','LineStyle',lines{3},'LineWidth',2,...
%                                            'DisplayName','Accepted lines');
%             end
%             hold off
%             title('Centre-line between defined start and end points')
%             legend
%         end
%     end
end