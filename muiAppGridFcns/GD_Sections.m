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
        ChannelTopo
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
            end
        end
    end
%%
    methods
        function setBoundary(obj,cobj,muicat)
            %load a grid and interactively define a bouindary line with
            %options to autogenerate and edit or digistise the line.
            gridclasses = {'EDBimport'}; %add other classes if needed
            grid = getGrid(cobj,1);   %grid for selected year

            

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
            clines = [];
            ok = 0;
            while ok<1
                [cline,hf] = gd_centreline(grid,mobj,props,clines); 
                if isempty(cline), ok = 1; continue; end
                answer = questdlg('Accept the centreline?','Centre-line','Yes','No','Yes');
                if strcmp(answer,'Yes')
                     %convert format of output if required
                    cline = gd_pnt2vec(cline,1);
                    clength = sum(vecnorm(diff(cline),2,2));

                    cpoints = round(clength/cint);
                    newcline = curvspace(cline,cpoints);
                    clines = [clines;newcline;[NaN,NaN]]; %#ok<AGROW> 
                end
                delete(hf)
            end

           answer = questdlg('Check/Edit the centre-lines?','Centre-line','Yes','No','Yes');
           if strcmp(answer,'Yes')
               paneltxt = 'Edit the centre-lines';
               promptxt = {'Select line to edit'};
               clines = gd_editlines(grid,paneltxt,promptxt,clines,true);
           end

            obj.ChannelLine = clines;
            

            cobj.Data.Grid.UserData.cline.x = clines(:,1);
            cobj.Data.Grid.UserData.cline.y = clines(:,2);

        end
%%
        function setSectionLines(obj,cobj,muicat)
            %use the centre-line to generate a set of sections at right
            %angles and interactively edit them, or digitise them manually.
            %
        end
%%
        function setSections(obj,cobj,muicat)
            %using the Boundary and SectionLines data extract the widths as
            %a function of elevation and distance along the centre-line
            %from the mouth.
        end
    end
%%
%     methods (Access=private)
%         function shplines = shapeLines(~,clines)
% 
%         end
% 
%     end
end