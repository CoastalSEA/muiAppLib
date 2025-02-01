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
        %
        Boundary
        ChannelLine
        % ChannelTopo
        SectionsLines
    end
    
    properties (Transient)
        Sections
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
            % mobj - mui model instance
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
                    obj = GD_Sections;
                    setBoundary(obj,cobj,muicat); 
                case 'Channel network'
                    promptxt = 'Select a Case to use to define channel network:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = GD_Sections;
                    setChannelNetwork(obj,cobj,mobj); 
                case 'Section lines'
                    promptxt = 'Select a Case to use to define section lines:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = GD_Sections;
                    setSectionLines(obj,cobj,muicat);
                case 'Sections'
                    promptxt = 'Select a Case to use to extract sections:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = GD_Sections;
                    setSections(obj,cobj,muicat);
                 case 'View Sections'
                    promptxt = 'Select a Case to use to extract sections:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    isnot = isempty(cobj) || ~isfield(cobj.Data,'Grid') ||...
                        ~isfield(cobj.Data.Grid.UserData,'sections') ||...
                        ~isa(cobj.Data.Grid.UserData.sections,'GD_Sections');
                    %
                    if isnot
                        warndlg('No Section data available to edit');
                        return;
                    end  
                    obj = cobj.Data.Grid.UserData.sections;
                    viewSections(obj,cobj,mobj);
            end
        end

%%
        function obj = getSections()
            obj = GD_Sections;
        end

%%
        function loadLines(mobj,src,gridclasses)
            %load linework for selected line type from a shapefile
            %can be Boundary, Channel Network or Section Lines

        end

%%
        function editLines(mobj,src,gridclasses)
            %edit linework for selected line type from a shapefile
            %can be Boundary, Channel Network or Section Lines
            muicat = mobj.Cases;   %handle to muiCatalogue
            switch src.Text
                case 'Boundary'
                    promptxt = 'Select a Case to use to define boundary:';
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    isnot = isempty(cobj) || ~isfield(cobj.Data,'Grid') ||...
                        ~isfield(cobj.Data.Grid.UserData,'sections') ||...
                        isempty(cobj.Data.Grid.UserData.sections.Boundary);
                    type = 'Boundary';
                case 'Channel network'
                    promptxt = 'Select a Case to use to define channel network:';
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    isnot = isempty(cobj) || ~isfield(cobj.Data,'Grid') ||...
                        ~isfield(cobj.Data.Grid.UserData,'sections') ||...
                        isempty(cobj.Data.Grid.UserData.sections.ChannelLine);
                    type = 'ChannelLine';
                case 'Section lines'
                    promptxt = 'Select a Case to use to define section lines:';
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    isnot = isempty(cobj) || ~isfield(cobj.Data,'Grid') ||...
                        ~isfield(cobj.Data.Grid.UserData,'sections') ||...
                        isempty(cobj.Data.Grid.UserData.sections.SectionLines);
                    type = 'SectionLines';
            end
            %
            if isnot
                warndlg(sprintf('No %s data available to edit',src.Text));
                return;
            end

            %eidt the selected line and save to case
            clines = cobj.Data.Grid.UserData.sections.(type);
            paneltxt = sprintf('Edit the %s',src.Text);
            grid = getGrid(cobj,1);   %grid for selected year
            clines = gd_editlines(grid,paneltxt,clines,true);
            answer = questdlg(sprintf('Save the edited %s',src.Text),src.Text,'Yes','No','Yes');
            if strcmp(answer,'Yes')
                cobj.Data.Grid.UserData.sections.(type) = clines;
                getdialog(sprintf('Edits saved for %s in Case: %s',src.Text,grid.desc),[],3)
            end
        end
    end
%%
    methods
        function setBoundary(obj,cobj,~)
            %load a grid and interactively define a boundary line with
            %options to autogenerate and edit or digistise the line.
            grid = getGrid(cobj,1);   %grid for selected year

            %clean-up shoreline
            paneltxt = 'Create and smooth the shoreline boundary';
            shore_xy = gd_boundary(grid,paneltxt,2,true);
            if isempty(shore_xy), return; end

            %give user opportunity to edit the shoreline interactively
           answer = questdlg('Check/Edit the shoreline?','Centre-line','Yes','Save','Quit','Yes');
           if strcmp(answer,'Yes')
               paneltxt = 'Edit the shoreline boundary';        
               shore_xy = gd_editlines(grid,paneltxt,shore_xy,true);
           elseif strcmp(answer,'Quit')
               return;
           end

           %save results to grid class instance 
           if isfield(cobj.Data.Grid.UserData,'sections')  && ...
                       isa(cobj.Data.Grid.UserData.sections,'GD_Sections')
                obj = cobj.Data.Grid.UserData.sections;
           end
            obj.Boundary = shore_xy;
            cobj.Data.Grid.UserData.sections = obj;
            getdialog(sprintf('Boundary added to grid for Case: %s',grid.desc),[],3)
        end
%%
        function setChannelNetwork(obj,cobj,mobj)
            %load a grid and extract the centre lines of the channels in
            %the network, or digitise them manually, and combine into a 
            %vector of points and a directed graph to define the topology
            if ~isfield(cobj.Data,'Grid'), return; end
            grid = getGrid(cobj,1);   %grid for selected year

            %get maximum water level to define
            promptxt = {'Maximum accessible water level?','Depth exponent','Sampling interval (m)'};
            defaults = {num2str(max(grid.z,[],'all')),'5','100'};
            inp = inputdlg(promptxt,'Water level',1,defaults);
            if isempty(inp), return; end %user cancelled
            props.maxwl = str2double(inp{1});
            props.dexp = str2double(inp{2});
            cint = str2double(inp{3});
            nlines = [];  %numerical [Nx2} array
            ok = 0;
            while ok<1
                nline = gd_centreline(grid,mobj,props,nlines(1:end-1,:)); 
                if isempty(nline), ok = 1; continue; end
                nline = cell2mat(struct2cell(nline)');       %struct to [Nx2] array
                hf = cl_checkPlot(obj,grid,nline,nlines);
                answer = questdlg('Accept the centreline?','Centre-line','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    %convert format of output if required                    
                    clength = sum(vecnorm(diff(nline),2,2));  %cline is a column vector [Nx2]
                    cpoints = round(clength/cint);            %number of points in new line
                    newcline = curvspace(nline,cpoints);
                    nlines = [nlines;newcline;[NaN,NaN]];     %#ok<AGROW> 
                end
                delete(hf)
            end
            if isempty(nlines), return; end                   %user cancelled without creating any lines
            nlines = nlines(1:end-1,:);                       %remove trailing NaNs
            clines.x = nlines(:,1); clines.y =  nlines(:,2);  %x,y struct of vector points
           %give user opportunity to edit the centre-line interactively
           answer = questdlg('Check/Edit the centre-lines?','Centre-line','Yes','Save','Quit','Yes');
           if strcmp(answer,'Yes')
               paneltxt = 'Edit the centre-lines';
               clines = gd_editlines(grid,paneltxt,clines,true);
           elseif strcmp(answer,'Quit')
               return;
           end

           %save results to grid class instance
           if isfield(cobj.Data.Grid.UserData,'sections')  && ...
                   isa(cobj.Data.Grid.UserData.sections,'GD_Sections')
               obj = cobj.Data.Grid.UserData.sections;
           end
           obj.ChannelLine = clines;
           cobj.Data.Grid.UserData.sections = obj;
           getdialog(sprintf('Channel network added to grid for Case: %s',grid.desc),[],3)
           %option to also add to grid cline property
           % cobj.Data.Grid.UserData.cline.x = clines(:,1);
           % cobj.Data.Grid.UserData.cline.y = clines(:,2);
        end

%%
        function setSectionLines(obj,cobj,~)
            %use the centre-line to generate a set of sections at right
            %angles and interactively edit them, or digitise them manually.
            if isfield(cobj.Data.Grid.UserData,'sections')  && ...
                    isa(cobj.Data.Grid.UserData.sections,'GD_Sections')
                obj = cobj.Data.Grid.UserData.sections;
            end
            obj.SectionsLines = slines;
            cobj.Data.Grid.UserData.sections = obj;
        end

%%
        function setSections(obj,cobj,~)
            %using the Boundary and SectionLines data extract the widths as
            %a function of elevation and distance along the centre-line
            %from the mouth.

        end

%%
        function viewSections(obj,cobj,mobj)
            %view section line work and select sections to plot

        end
    end
%%
    methods (Access=private)
        function hf = cl_checkPlot(~,grid,cline,nlines)
            %plot base map of initial grid selection and defined mask
            % cline and nlines are both [Nx2] arrays
            points(1,:) = cline(1,:);
            points(2,:) = cline(end,:);
            hf = figure('Name','Thalwegs','Units','normalized','Tag','PlotFig');  
            hf.Position = [0,0,1,1];
            ax = gd_plotgrid(hf,grid);
            axis equal  %assume geographical projection or grid of similar dimensions
            axis tight
            colormap(ax,'gray');
            lines = {'-','--',':','-.'};
            
            hs = findobj(ax.Children,'Type','surface');
            hs.Annotation.LegendInformation.IconDisplayStyle = 'off';
            hold on
            hp = plot(ax,points(:,1),points(:,2),'ok','MarkerSize',8,...
                                   'MarkerFaceColor','w','Tag','mypoints');
            hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            plot(ax,cline(:,1),cline(:,2),'r','LineStyle',lines{1},'LineWidth',2,...
                                           'DisplayName','New line');
            if ~isempty(nlines)
                plot(ax,nlines(:,1),nlines(:,2),'g','LineStyle',lines{3},'LineWidth',2,...
                                           'DisplayName','Accepted lines');
            end
            hold off
            title('Centre-line between defined start and end points')
            legend
        end
    end
end