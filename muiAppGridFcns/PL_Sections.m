classdef PL_Sections < handle   
%
%-------class help---------------------------------------------------------
% NAME
%   PL_Sections.m
% PURPOSE
%   Class to extract sections from a grid
% USAGE
%   obj = PL_Sections.sectionsMenu(mobj,src,gridclasses);
%         %mobj is a handle to Main UI
%         %src is handle to calling menu option
%         %gridclasses is cell array of classes to use in case selection
% NOTES
%   calls classes that inherit PLinterface to manipulate different types of
%   linework used to generate, plot and analyse cross-sections
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
                            %fields: maxwl,dexp,cint,ChannelLengths,Network
        SectionLines        %lines that define cross-sections
        XSections           %cross-sections obtained from grid
    end
       
    methods (Access=protected)
        function obj = PL_Sections()          
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
            % where Load reads a shapefile, Generate uses the PL_Sections
            % methods and Edit calls gd_editlines (see EstauryDB for eg)
            %
            muicat = mobj.Cases;   %handle to muiCatalogue

            if  any(strcmp({'Layout','Sections','Network'},src.Text))
                srcText = src.Text;
                src = src.Parent;
            end

            switch src.Text
                case 'Boundary'
                    promptxt = 'Select a Case to use to define boundary:';
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = PL_Sections.getSections(cobj);
                    setBoundary(obj,cobj,muicat);
                case 'Channel Network'
                    promptxt = 'Select a Case to use to define channel network:';
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = PL_Sections.getSections(cobj);
                    setChannelNetwork(obj,cobj,mobj);
                case 'Channel Links'
                    promptxt = 'Select a Case to use to set channel connectivity:';
                    [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = PL_Sections.getSections(cobj);
                    setChannelTopology(obj,cobj,catrec);                    
                case 'Section Lines'
                    promptxt = 'Select a Case to use to define section lines:';
                    [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = PL_Sections.getSections(cobj);
                    setSectionLines(obj,cobj,catrec);
                case 'Sections'
                    promptxt = 'Select a Case to use to extract sections:';
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = PL_Sections.getSections(cobj);
                    setSections(obj,cobj,muicat);
                case 'View Sections'
                    promptxt = 'Select a Case to use to view section linework:';
                    [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj.Sections) || ~isa(cobj.Sections,'PL_Sections')
                        warndlg('No Section data available to view');
                        return;
                    end
                    obj = cobj.Sections;
                    viewSections(obj,cobj,catrec,srcText);
                case 'Waterbody'
                    promptxt = 'Select a Case to use to define waterbody';
                    [cobj,caserec] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    obj = PL_Sections.getSections(cobj);
                    setWaterbody(obj,cobj,muicat,classRec(muicat,caserec));                   
            end
        end

%%
        function obj = getSections(cobj)
            %retrieve existing class object or create a new instance
            if nargin<1, cobj = []; end

            if isa(cobj.Sections,'PL_Sections')
                obj = cobj.Sections;
            else                
                obj = PL_Sections;
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
                case 'Channel Network'
                    type = 'ChannelLine';
                case 'Section Lines'
                    type = 'SectionLines';
                case 'Waterbody'
                    type = 'WaterBody' ;
                otherwise
                    return
            end
            promptxt = sprintf('Select a Case to load %s from shapefile:',linetype);           
            [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(cobj), return; end

            [fname,path,nfiles] = getfiles('MultiSelect','off',...
                'FileType',{'*.shp;'},'PromptText','Select shape file):');
            if nfiles~=1, return; end

            %get file and read the data
            shp = gd_readshapefile(path,fname);
            if isrow(shp.x)
                shp = structfun(@transpose,shp,'UniformOutput',false);
            end    
            if strcmp(type,'WaterBody')
                cobj.WaterBody = shp;
            else
                %initialise PL_Sections instance if not already available
                obj = PL_Sections.getSections(cobj);
                obj.(type) = shp;  %overwrites any existing lines
                cobj.Sections = obj;
            end
            
            casedesc = catrec.CaseDescription;
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
                case 'Channel Network'
                    type = 'ChannelLine';
                case 'Section Lines'
                    type = 'SectionLines';
                case 'Waterbody'
                    type = 'WaterBody' ;
                otherwise
                    return
            end
            
            promptxt = sprintf('Select a Case to edit %s:',linetype);
            [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(cobj), return; end
            casedesc = catrec.CaseDescription;

            if ~isa(cobj.Sections,'PL_Sections') || isempty(cobj.WaterBody)
                warndlg(sprintf('No %s data available to edit',linetype));
                return;
            elseif strcmp(type,'WaterBody')
                clines = cobj.WaterBody;
            else
                clines = cobj.Sections.(type);
            end

            %edit the selected line and save to case            
            promptxt = sprintf('%s the %s\nSelect menu option',src.Text,linetype);
            grid = getGrid(cobj,1);   %grid for selected year

            clines = PL_Editor.Figure(grid,promptxt,clines,false,true);
            if isempty(clines), return; end
            answer = questdlg(sprintf('Save the edited %s',linetype),linetype,'Yes','No','Yes');
            if strcmp(answer,'Yes')
                if strcmp(type,'WaterBody')
                    cobj.WaterBody = clines;
                else
                    cobj.Sections.(type) = clines;                
                end
                getdialog(sprintf('Edits saved for %s in Case: %s',linetype,casedesc),[],3)
            end
        end

%%
        function deleteLines(mobj,src,gridclasses)
            %delete linework for the selected line type
            muicat = mobj.Cases;   %handle to muiCatalogue
            linetype = src.Parent.Text;
            switch linetype
                case 'Boundary'
                    type = 'Boundary';
                case 'Channel Network'
                    type = 'ChannelLine';
                case 'Section Lines'
                    type = 'SectionLines';
                case 'Waterbody'
                    type = 'WaterBody';
                otherwise
                    return
            end        
            
            promptxt = sprintf('Select a Case to delete %s:',linetype);
            [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(cobj), return; end
            
            %delete the selected line and save to case    
            qprompt = sprintf('Delete %s for Case %s',linetype,catrec.CaseDescription);
            answer = questdlg(qprompt,linetype,'Yes','No','Yes');
            if strcmp(answer,'Yes')
                if strcmp(type,'WaterBody')
                    cobj.WaterBody = [];
                else
                    cobj.Sections.(type) = [];
                end
                getdialog(sprintf('%s deleted for Case: %s',linetype,...
                                          catrec.CaseDescription),[],3)
            end
        end

%%
        function view_WBlines(mobj,src,gridclasses)
            %plot the waterbody boundary
            muicat = mobj.Cases;   %handle to muiCatalogue
            promptxt = sprintf('Select a Case to view %s:',src.Parent.Text);
            [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(cobj) || isempty(cobj.WaterBody)
                warndlg('No Waterbody data available to view');
                return
            end
            lines = cobj.WaterBody;
            hf = figure('Name','Waterbody','Units','normalized','Tag','PlotFig');                                          
            ax = PL_Sections.getGrid(cobj,hf); 
            hold on
            plot(ax,lines.x,lines.y,'-r','DisplayName','Waterbody',...
                                               'ButtonDownFcn',@godisplay);
            hold off
            %add lablels and title
            xlabel('Eastings (m)'); 
            ylabel('Northings (m)');    
            title(catrec.CaseDescription); 
        end


%%
        function exportLines(mobj,src,gridclasses)
            %export the selected line type to a text file or shape file
            %for a text file the format is an x,y struct with lines
            %separated by NaNs.
            muicat = mobj.Cases;   %handle to muiCatalogue
            promptxt = sprintf('Select a Case to Export %s:',src.Parent.Text);
            [cobj,~,catrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(cobj), return; end

            linetype = src.Parent.Text;
            isSXN = ~isempty(cobj.Sections);
            isWBD = ~isempty(cobj.WaterBody);
            lines = [];
            switch linetype
                case 'Boundary'
                    if isSXN && ~isempty(cobj.Sections.Boundary)
                        lines = cobj.Sections.Boundary;
                    end
                case 'Channel Network'
                    if isSXN && ~isempty(cobj.Sections.ChannelLine)
                        lines = cobj.Sections.ChannelLine;
                    end
                case 'Section Lines'
                    if isSXN && ~isempty(cobj.Sections.SectionLines)
                        lines = cobj.Sections.SectionLines;
                    end
                case 'Waterbody'
                    if isWBD && ~isempty(cobj.WaterBody)
                        lines = cobj.WaterBody;
                    end
                otherwise                    
                    return
            end
            %
            if isempty(lines)
                warndlg(sprintf('No %s data for case %s',linetype,catrec.CaseDescription))
                return
            end

            answer = questdlg('Export format to use:','Export','txt','mat','shp','txt');
            fname = matlab.lang.makeValidName(sprintf('%s_%s',linetype,catrec.CaseDescription));
            fname = inputdlg({'Filename for export:'},'Export',[1 60],{fname},'on');

            lines.Source = 'CoastalSEA EstuaryDB App';
            lines.Filename = mobj.Info.FileName;
            lines.CaseIndex = catrec.CaseID;
            lines.CaseDescription = catrec.CaseDescription;
            lines.LineType = linetype;
            if strcmp(answer,'txt')
                fileID = fopen([fname{1},'.txt'],'w');
                fprintf(fileID,'Source: %s\n',lines.Source);
                fprintf(fileID,'Project filename: %s\n',lines.Filename);
                fprintf(fileID,'CaseID: %s\n',num2str(lines.CaseIndex));
                fprintf(fileID,'Description: %s\n',lines.CaseDescription);
                fprintf(fileID,'Line Type: %s\n',lines.LineType);
                fprintf(fileID,'%6s %12s\n','x','y');
                A = [lines.x,lines.y];
                fprintf(fileID,'%f %f\n',A');
                fclose(fileID);
            elseif strcmp(answer,'mat')
                save(fname{1},'lines','-mat');
            else
                istoolbox = false;
                if check4toolbox('MAP_Toolbox')                              %toolbox is licensed to use
                    addons = matlab.addons.installedAddons;
                    istoolbox = any(matches(addons.Name,'Mapping Toolbox')); %toolbox is installed
                end 
                if ~istoolbox
                    warndlg('Mapping Toolbox is not available. Cannot export as a shape file')
                    return; 
                end
                s = mapshape(lines.x,lines.y);
                s.Metadata.Source = 'CoastalSEA EstuaryDB App';
                s.Metadata.Filename = mobj.Info.FileName;
                s.Metadata.CaseIndex = catrec.CaseID;
                s.Metadata.CaseDescription = catrec.CaseDescription;
                s.Metadata.LineType = linetype;
                shapewrite(s,fname{1});
            end
            getdialog(sprintf('Data written to %s',fname{1}))
        end
%%
        function ax = getGrid(cobj,hf)
            isgrid = false; isimage = false;
            if isfield(cobj.Data,'Grid')
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
                if any(isnan(im.CData))
                    set(h_im, 'AlphaData', 1-isnan(im.CData)); %set Nan values to be transparent 
                end
                colormap(im.CMap);
                clim(im.CLim);
                isimage = true;
            elseif isfield(cobj.Data,'Image') %&& isfield(cobj.Data.Image,'XData')
                dst = cobj.Data.Image;
                im = dst.image;                      %image object
                ax = axes(hf);
                imagesc(ax,im{1})
                isimage = true;
            else
                ax = axes(hf);
            end
            %
            if isgrid || isimage 
                axis equal tight
                colormap(ax,'gray');
                cb = colorbar;
                cb.Label.String = 'Elevation (mAD)';                 
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
                answer = questdlg('A boundaray line exists. Modify existing or create new one?',...
                                  'Channel','Modify','New','Copy','Modify');
                if strcmp(answer,'Modify')
                    blines = obj.Boundary;      %existing lines
                elseif strcmp(answer,'Copy') 
                    blines = cobj.WaterBody;
                end    
            else
                answer = questdlg('Copy Waterbody boundary or create a new boundary?',...
                                  'Channel','New','Copy','New');
                if strcmp(answer,'Copy') 
                    blines = cobj.WaterBody;
                end  
            end 

            %clean-up shoreline
            promptxt = sprintf('Create and smooth the shoreline boundary\nSelect menu option');
            %shore_xy = gd_setboundary(grid,paneltxt,blines,true);
            shore_xy = PL_Boundary.Figure(grid,promptxt,blines,true);
            if isempty(shore_xy), return; end

            %save results to grid class instance
            obj.Boundary = shore_xy;
            cobj.Sections = obj;
            getdialog(sprintf('Boundary added to grid for Case: %s',grid.desc),[],3)
        end
%%
        function setChannelNetwork(obj,cobj,~)
            %load a grid and extract the centre lines of the channels in
            %the network, or digitise them manually, and combine into a
            %vector of points
            if ~isfield(cobj.Data,'Grid')
                warndlg('No grid for selected case'); return;
            end
            grid = getGrid(cobj,1);   %grid for selected year
            
            if ~isempty(obj.ChannelLine)           %channel line exists
                promptxt = sprintf('A centre-line exists.\nModify existing centre-line or create a New one?');
                answer = questdlg(promptxt,'Channel','Modify','New','Modify');                                 
            else
                answer = 'New';
            end   

            if strcmp(answer,'New') || isempty(obj.ChannelProps)
            %get maximum water level to define limits of search
                promptxt = {'Maximum accessible water level?','Depth exponent','Sampling interval (m)'};
                defaults = {num2str(max(grid.z,[],'all')),'5','100'};
                inp = inputdlg(promptxt,'Water level',1,defaults);
                if isempty(inp), return; end %user cancelled
                props.maxwl = str2double(inp{1});
                props.dexp = str2double(inp{2});
                props.cint = str2double(inp{3});
                obj.ChannelProps = props;
                obj.ChannelLine = [];
            end

            promptxt = sprintf('Create a channel network of centre-lines\nSelect menu option');
            [clines,props] = PL_CentreLine.Figure(grid,promptxt,obj,true);
            if isempty(clines), return; end   %user cancelled without creating any lines
            
            %save linework so that can reload if mess up defining topology           
            obj.ChannelLine = clines;
            obj.ChannelProps = props;

            %save results to grid class instance
            cobj.Sections = obj;
            getdialog(sprintf('Channel network added to grid for Case: %s',grid.desc),[],3)

            %option to also add to grid cline property
            % cobj.Data.Grid.UserData.cline.x = clines(:,1);
            % cobj.Data.Grid.UserData.cline.y = clines(:,2);
        end

%%
        function setSectionLines(obj,cobj,catrec)
            %use the centre-line to generate a set of sections at right
            %angles and interactively edit them, or digitise them manually.
            %at the moment this requires a grid to be present but not needed ** 
            if isempty(obj.Boundary) || isempty(obj.ChannelLine)
                warndlg('Boundary and Channel network need to be defined to extract Section Lines')
                return
            end

            %use either a grid or a grid image
            if isfield(cobj.Data,'Grid')
                grid = getGrid(cobj,1);   %grid for selected year
            elseif isfield(cobj.Data,'GeoImage')
                grid = cobj.Data.GeoImage;
            else
                warndlg('No grid or image data for selected case');
            end

            %extract the section lines that are normal to the 
            %channel centre line and extend to the bounding shoreline
            promptxt = sprintf('Extract the channel cross-sections\nSelect menu option');            
            [slines,clines] = PL_SectionLines.Figure(grid,promptxt,obj,true);
            if isempty(slines), return; end

            obj.SectionLines = slines;
            obj.ChannelLine = clines;
            cobj.Sections = obj;
            casedesc = catrec.CaseDescription;
            msgtxt = sprintf('Section lines added to grid for Case: %s',...
                                                                casedesc);
            getdialog(msgtxt,[],3)
        end

%%
        function setSections(obj,cobj,~)
            %using the SectionLines data extract the widths as a function 
            %of elevation for all sections along the channels
            if isempty(obj.SectionLines) || ~isfield(cobj.Data,'Grid')
                getdialog('No data',[],1); return; 
            end
            grid = getGrid(cobj,1);
            %convert the SectionLines to cplines
            cplines = gd_plines2cplines(gd_lines2points(obj.SectionLines));
            %use grid to interpolate sections and plot them
            hwb = progressbar([],'Computing centre-line');
            xlines= gd_plotsections(grid,cplines);
            progressbar(hwb);
            if ~isempty(xlines)
                obj.XSections = xlines;
                cobj.Sections = obj;
            end
        end

%%
        function setChannelTopology(obj,cobj,catrec)
            %define the toplogy for a channel network            
            if isempty(obj.ChannelLine)
                warndlg('Channel network need to be defined to extract Section Lines')
                return
            end
            %use either a grid or a grid image
            if isfield(cobj.Data,'Grid')
                grid = getGrid(cobj,1);   %grid for selected year
            elseif isfield(cobj.Data,'GeoImage')
                grid = cobj.Data.GeoImage;
            else
                warndlg('No grid or image data for selected case');
            end
            
            if ~isempty(obj.SectionLines)
                %use section lines to update centre-line positions
                % s_plines = gd_points2lines(obj.SectionLines);
                % c_plines = gd_points2lines(obj.ChannelLine);
                plines = resetCentreLine(obj);
                if isempty(plines), return; end
                obj.ChannelLine = gd_points2lines(plines,2);
            end

            %generate the toplogy based on updated
            [cumlen,G,hf,~] = gd_linetopology(grid,obj.ChannelLine);
            obj.ChannelProps.Network = G;
            obj.ChannelProps.ChannelLengths = cumlen;
            delete(hf)            %delete the figure but not network graph
            casedesc = catrec.CaseDescription;
            msgtxt = sprintf('Section lines added to grid for Case: %s',...
                                                                casedesc);
            getdialog(msgtxt,[],3)
        end


%%
        function nc_plines = resetCentreLine(obj)
            %update the centre-line so that the points lie on the sections and
            %return updated points  
            s_plines = gd_lines2points(obj.SectionLines);
            c_plines = gd_lines2points(obj.ChannelLine);
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
                    % cp = [];
                    P = InterX(sline',cline'); %intersections for 1 line
                    %line intersection misses points at the start and end of each
                    %centre-line line. Check if the point is a least on the section.
                    if isempty(P)
                        %see if section is on the first point of the centre-line
                        ison = ispointonline(sline(1:end-1,:)',cline(1,:)',1,1e2);
                        if ison
                            cp = gd_lines2points(cline(1,:));
                            plines = [plines,cp];             %#ok<AGROW> 
                        else  
                            
                            ison = ispointonline(sline(1:end-1,:)',cline(end-1,:)',1,1e2);
                            if ison
                                cp = gd_lines2points(cline(end-1,:));
                                plines = [plines,cp];         %#ok<AGROW>
                            end
                        end
                    else
                        cp = struct('x',P(1,1),'y',P(2,1));
                        plines = [plines,cp];                 %#ok<AGROW> 
                        
                    end

                    % if ~isempty(cp)
                    % hold on, plot(cp.x,cp.y,'ok')
                    % end
                end
                cplines{1,j} = plines;                        %#ok<AGROW> 
            end
            %reinstate NaN line terminators
            nanpnts.x = NaN; nanpnts.y = NaN;                 %line termination
            nc_plines = [];
            for j=1:nlines
                nc_plines = [nc_plines,cplines{1,j},nanpnts]; %#ok<AGROW> 
            end

            %check plot
            hf = figure('Name','Topology','Units','normalized','Tag','PlotFig');  
            hf.Position = [0.05,0.05,0.9,0.9];
            ax = axes(hf);
            gd_plotpoints(ax,s_plines,'mylines',2);    %plot section lines
            hold on
            gd_plotpoints(ax,nc_plines,'mypoints',1);   %plot new centre-line as points
            hold off
            axis equal 

            answer = questdlg('Is there one point for every section?',...
                                              'Topology','Yes','No','Yes');
            if strcmp(answer,'No')
                warndlg('Try regenerating/editing centre-line and/or sections')
                nc_plines = [];   %return centre-line unchanged                
            end
            delete(hf)
        end


%%
        function setWaterbody(obj,cobj,muicat,classrec)
            %use PL_Boundary to create a polygon boundary
            if ~isfield(cobj.Data,'Grid')
                warndlg('No grid for selected case'); return;
            end
            grid = getGrid(cobj,1);             %grid for estuary

            wlines = 2;                         %output format when new line
            if ~isempty(cobj.WaterBody)         %waterbody line exists
                answer = questdlg('A waterbody polygon exists. Edit existing or create a new one?',...
                                  'Channel','Modify','New','Copy','Modify');
                if strcmp(answer,'Modify')
                    wlines = cobj.WaterBody;    %existing lines
                elseif strcmp(answer,'Copy') 
                    wlines = obj.Boundary;      %copy Sections boundary
                end 
            else
                answer = questdlg('Copy Sections boundary or create a new waterbody?',...
                                  'Channel','New','Copy','New');
                if strcmp(answer,'Copy') 
                    wlines = obj.Boundary;      %copy Sections boundary
                end 
            end 

            %clean-up shoreline
            promptxt = sprintf('Create and smooth the shoreline boundary\nSelect menu option');
            %shore_xy = gd_setboundary(grid,paneltxt,blines,true);
            wline = PL_Boundary.Figure(grid,promptxt,wlines,true);
            if isempty(wline), return; end

            cobj.WaterBody = wline;
            updateCase(muicat,cobj,classrec);
        end
%%
%--------------------------------------------------------------------------
% Plot functions
%--------------------------------------------------------------------------
        function viewSections(obj,cobj,catrec,srcText)
            %view boundary channel network and cross-sections line work
            %or along-channel sections summary plot
            casedesc = catrec.CaseDescription;
            switch srcText
                case 'Layout'
                    viewPlanSections(obj,cobj,casedesc);
                case 'Sections'                    
                    viewAlongChannelSections(obj,casedesc);
                case 'Network'
                    viewChannelNetwork(obj,casedesc);
            end
        end

%%
        function ax = viewPlanSections(obj,cobj,casedesc,hf)
            %plot boundary channel network and cross-sections line work
            if nargin<4
                hf = figure('Name','Sections','Units','normalized',...
                                         'Tag','PlotFig','Visible','off');
            end

            ax = PL_Sections.getGrid(cobj,hf); 

            type = {'Boundary','SectionLines','ChannelLine'};
            for i=1:3
                switch type{i}
                    case 'Boundary'
                        green = mcolor('green');
                        sc = green; ss = '-'; sw = 2;
                    case 'ChannelLine'
                        sc = 'b'; ss = ':'; sw = 1.5;
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
                        cplines = gd_plines2cplines(gd_lines2points(lines)); 
                        for j=1:length(cplines)
                            pline = cplines{1,j};
                            [ax,H] = gd_plotpoints(ax,pline,num2str(j),3); 
                            H(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
                        end
                    end
                    hold off
                end
            end

            %add lablels and legend
            xlabel('Eastings (m)'); 
            ylabel('Northings (m)');    
            title(casedesc);    
            legend
            hf.Visible = 'on';
        end

%%
        function viewAlongChannelSections(obj,casedesc)
            %plot along-channel sections
            hf = figure('Name','Sections','Units','normalized',...
                             'Tag','PlotFig','Visible','on');
            subplot(axes(hf));
            glines = {'-','--',':','-.'};

            plines = obj.XSections;
            cplines = gd_plines2cplines(gd_lines2points(plines));

            x = obj.ChannelProps.ChannelLengths;
            idN = [0,find(isnan(x))];
            idN = arrayfun(@(x,i) x-(i-1),idN,1:length(idN));
            nreach = length(idN)-1;
            for i=1:nreach
                si = subplot(ceil(nreach/2),2,i);
                hold on
                idx = idN(i)+1:idN(i+1);
                nlines = length(idx);
                for k=1:nlines
                    aline = cplines{1,idx(k)};
                    spnts = [aline(:).x];
                    zpnts = [aline(:).y];
                    
                    lname = sprintf('Reach %d, Section %d',i,k);        
                    plot(si,spnts,zpnts,'LineStyle',glines{rem(k,4)+1},...
                          'LineWidth',1,'Tag','asection','DisplayName',lname,...
                          'ButtonDownFcn',@godisplay)
                end
                hold off 

                title(sprintf('Reach %d',i))        
                
            end
            sgtitle(sprintf('Along-channel sections for %s',casedesc))
        end

%%
        function viewChannelNetwork(obj,casedesc)
            %plot the directed graph of the channel network
            hf = figure('Name','Network','Units','normalized',...
                             'Tag','PlotFig','Visible','on');
            ax = axes(hf);
            G = obj.ChannelProps.Network;
            nlabel = strcat(cellstr(G.Nodes.Names),{' ('},...
                                        cellstr(G.Nodes.Distance),{'m)'});
            plot(ax,G,'EdgeLabel',G.Edges.Weight,'NodeLabel',nlabel)
            title(sprintf('Channel network for %s',casedesc))   
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