classdef GD_ImportData < FGDinterface
%
%-------class help---------------------------------------------------------
% NAME
%   GD_ImportData.m
% PURPOSE
%   Class to import grid data, adding the results to dstable
%   and a record in a dscatlogue (as a property of muiCatalogue)
% USAGE
%   obj = GD_ImportData()
% SEE ALSO
%   inherits GDinterface and muiDataSet and is part of the grid tools used 
%   by mui Apps. FGD_ImportData is the same but inherits FGDinterface which
%   has a range of additional estuary/inlet specific tools.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%    
    properties  
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:
    end
    
    methods 
        function obj = GD_ImportData()               
            %class constructor                          
        end
    end
%%
    methods (Static)
        function loadData(muicat,~)
            %read and load a data set from a file
            % mobj - handle to modelui instance 
            % classname - name of class being loaded
            obj = GD_ImportData; %initialise instance of class
            
            [fname,path,nfiles] = getfiles('MultiSelect','on',...
                'FileType','*.txt;*.grd;*.csv;*.xyz;*.mat','PromptText','Select file(s):');
            if nfiles==0
                return;
            elseif ~iscell(fname)
                fname = {fname};   %single select returns char
            end
            
            answer = questdlg('Is grid a channel/estuary/inlet?','Type','Yes','No','No');
            ischannel = strcmp(answer,'Yes');

            if nfiles>1
                timetxt = {num2str(0:1:nfiles-1),'years'};
            else
                timetxt = {'0','years'};
            end
            ok = 0;
            promptxt = {'Define grid intervals:','Units'};
            while ok<1
                timetxt = inputdlg(promptxt,'Load grid',1,timetxt);
                if isempty(timetxt), return; end
                %
                if strcmp(timetxt{2},'years')
                    timesteps = years(str2double(split(timetxt{1})));
                else
                    timesteps = str2double(split(timetxt{1}));
                end
                %check length is correct and that all values are unique
                if length(timesteps)==nfiles && isunique(timesteps)
                    ok = 1; 
                end                
            end
            
            %see if grid needs flipping or rotating
            filename{1,1} = [path,fname{1}];
            [~,~,ext] = fileparts(fname{1});
            if strcmp(ext,'.mat')
                promptxt = {'X grid-spacing','Y grid-spacing'};
                gridints = inputdlg(promptxt,'Input',1,{'2','2'});
                data = readmatfile(filename{1,1},gridints);
            else
                data = readinputfile(filename{1,1},1);  %header defines file read format
            end
            if isempty(data), return; end 
            data{3} = setDataRange(obj,data{3});

            grid = formatGridData(obj,data);%assign data to struct for x,y,z 
            if isempty(grid), return; end         %user deleted orientGrid UI
            [grid,rotate] = orientGrid(obj,grid); %option to flip or rotate grid
            if isempty(grid), return; end         %user deleted orientGrid UI
            newgrid(1,:,:) = grid.z;

            %now load each file selected
            nhead = 1; %number of header lines in file   
            for i=2:nfiles
                filename{i,1} = [path,fname{i}];
                if strcmp(ext,'.mat')
                    data = readmatfile(filename{i,1},gridints);
                else
                    data = readinputfile(filename{i,1},nhead);
                end
                if isempty(data), continue; end
                data{3} = setDataRange(obj,data{3});
                grid = formatGridData(obj,data); %assign data to struct for x,y,z
                grid = rotateGridData(obj,grid,rotate);
                newgrid(i,:,:) = grid.z;
            end
            %default values used when not a channnel or not required
            %alternative would be to have a class that inherits GDinterface
            dims = setChannel(obj,grid,timesteps,ischannel);
                                             
            %assign metadata about data source and save grid
            meta.source = filename;
            meta.data = sprintf('Rotate option = %d',rotate);
            obj = setGrid(obj,{newgrid},dims,meta);
            
            %setDataRecord classobj, muiCatalogue obj, dataset, classtype
            setCase(muicat,obj,'grid_data');
        end  
    end
%%
    methods
        function addData(obj,~,~,muicat) 
            %add additional data to an existing user dataset (called from
            %model UI using useCase in muiCatalogue
            [fname,path,~] = getfiles('MultiSelect','off',...
                'FileType','*.txt;*.grd;*.csv;*.xyz;*.mat','PromptText','Select file:');
            nhead = 1;
            filename = [path,fname];
            [~,~,ext] = fileparts(fname);
            if strcmp(ext,'.mat')
                promptxt = {'X grid-spacing','Y grid-spacing'};
                gridints = inputdlg(promptxt,'Input',1,{'2','2'});
                data = readmatfile(filename,gridints);
            else
                data = readinputfile(filename,nhead);  
            end
            if isempty(data), return; end
            
            %existing loaded grid
            dst = obj.Data.Grid;      %selected dstable
            x = dst.Dimensions.X;
            y = dst.Dimensions.Y;
            existimes = dst.RowNames;
            ok = 1;
            while ok>0  %add timestep and ensure is unique from existing values
                timestep = inputdlg('Define grid interval:','Load grid',1,{'X'});
                if isempty(timestep), return; end
                %
                if isduration(existimes)
                    timestep = years(str2double(timestep));
                else
                    timestep = str2double(timestep);
                end
                if ~ismember(existimes,timestep), ok=0; end %ensure value entered is unique
            end
            
            %check min and maximum values of grid being loaded
            data{3} = setDataRange(obj,data{3});
            
            %plot new grid and allow user to change orientation
            %NB this does not check that changes made match existing grids
            grid = formatGridData(obj,data); %assign data to struct for x,y,z 
            [grid,~] = orientGrid(obj,grid); %option to flip or rotate grid
            newgrid(1,:,:) = grid.z;
            
            %check that dimemsnions of new and exiting grids are the same
            if length(x)~=length(grid.x) || length(y)~=length(grid.y)
                warndlg('Dimensions of grid being added do not match existing grid')
                return
            end

            answer = questdlg('Is grid a channel/estuary/inlet?','Type','Yes','No','No');
            ischannel = strcmp(answer,'Yes');
            dims = setChannel(obj,grid,timestep,ischannel);

            %add grid to existing Case table
            addGrid(obj,muicat,newgrid,timestep,dims,filename,true);
        end        
%%
        function qcData(obj,classrec,catrec,muicat) %#ok<INUSD>
            %quality control a dataset
            warndlg('No quality control defined for this format in GD_ImportData');
        end  
%%
        function tabPlot(obj,src)
            %generate plot for display on Q-Plot tab
            if ~strcmp(src.Tag,'Plot') && ~strcmp(src.Tag,'FigButton')
                if isfield(obj.Data,'Hypsometry')
                    cf_model_tabs(obj,src);   
                    return;
                end
%                 return;
            end
            
            dst = obj.Data.Grid;
            
            if height(dst)>1
                %propmpt user to select timestep
                list = dst.DataTable.Properties.RowNames;
                irec = listdlg('PromptString','Select grid:',...
                               'Name','Tab plot','SelectionMode','single',...
                               'ListSize',[160,160],'ListString',list);
                if isempty(irec), return; end
            else
                irec = 1;
            end
            
            grid = getGrid(obj,irec);
            if isempty(grid.z), return; end

            %set >Figure button and create axes
            if strcmp(src.Tag,'Plot') || strcmp(src.Tag,'FigButton')
                tabcb  = @(src,evdat)tabPlot(obj,src);            
                ax = tabfigureplot(obj,src,tabcb,false);
                ax.NextPlot = 'add';
            else
                ax = src; %user passing an axis as src rather than a uicontrol
            end

            %plot form as a contour plot
            contourf(grid.x,grid.y,grid.z');
            ax = gd_ax_dir(ax,grid.x,grid.y);
            colormap('parula')
            shading flat
            %gd_colormap([min(grid.z,[],'all'),max(grid.z,[],'all')])
            %caxis([-8,2]);
            cb = colorbar;
            cb.Label.String = 'Elevation (mAD)';
            xlabel('Length (m)'); 
            ylabel('Width (m)'); 
            
            if isduration(grid.t)
                gridnum =char(grid.t);
            else
                gridnum = num2str(grid.t);
            end
            title(sprintf('%s (%s: %s)',dst.Description,dst.RowDescription,gridnum));
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot
        end       
    end
%%
    methods (Access = private)
        function grid = formatGridData(~,data)
            %format grid to load into a dstable as single row per grid
            %saves data as z array [1,Nx,Ny] for each time step
            %X and Y assumed fixed and saved as dimensions 
            grid.x = unique(data{:,1},'stable'); %stable so that orientatioin of 
            grid.y = unique(data{:,2},'stable'); %grid is preserved to input direction
            z = data{:,3};
            try
                Nx = length(grid.x);
                Ny = length(grid.y);
                grid.z = reshape(data{:,3},Nx,Ny);
            catch
                try
                    %[grid.x,grid.y,grid.z] = xyz2grid(data{:,1},data{:,2},data{:,3});
                    [grid.x,~,xi] = unique(data{:,1},'sorted'); 
                    [grid.y,~,yi] = unique(data{:,2},'sorted'); 
                    grid.z = accumarray([xi yi],z(:),[],[],NaN); 
                    %grid.z = flipud(Z);
                catch
                    warndlg('Unable to resolve grid from data')
                    grid = [];
                end
            end
        end
%%
        function zdata = setDataRange(~,zdata)
            %display z range of data and allow user to set new bounds
            minz = num2str(min(zdata));
            maxz = num2str(max(zdata));
            defaults = {minz,maxz};
            promptxt = {'Minimum z value','Maximum z value'};
            answer = inputdlg(promptxt,'Data range',1,defaults);
            if isempty(answer), return; end %user cancelled, limits unchanged
            
            minz = str2double(answer{1});
            maxz = str2double(answer{2});
            zdata(zdata<minz) = NaN;
            zdata(zdata>maxz) = NaN;
        end
%%        
        function dims = setChannel(~,grid,timesteps,ischannel)
            %default values used when not a channnel or not required
            %alternative would be to have a class that inherits GDinterface
            ishead = false; xM = 0;  Lt = max(grid.x);
            if ischannel
                answer = questdlg('Minimum X is mouth?','Load grid',...
                                               'True','False','N/A','True');            
                if strcmp(answer,'False')
                    ishead = true;
                end

                promptxt = {'Distance to mouth from grid origin (from min x)',...
                            'Distance from mouth to head (default is max-min x)'};
                defaults = {'0',num2str(max(grid.x))};
                inp = inputdlg(promptxt,'Load grid',1,defaults);
                if ~isempty(inp)
                    xM = str2double(inp{1});
                    Lt = str2double(inp{2});
                end
            end   
            
            %assume that the X and Y dimensions and xM are the same in all files
            Rv = struct('Hr',[],'Wr',[],'Ar',[]);
            dims = struct('x',grid.x,'y',grid.y,'t',timesteps,'xM',xM,...
                                         'Lt',Lt,'Rv',Rv,'ishead',ishead);          
        end
    end
end