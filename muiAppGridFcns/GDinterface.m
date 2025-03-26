classdef (Abstract = true) GDinterface < muiDataSet
%
%-------class help------------------------------------------------------
% NAME
%   GDinterface.m
% PURPOSE
%   Abstract class providing additional methods for use with gridded data sets
%   (eg imported or model data) to extend the functionality of the muiDataSet
%   abstract class.
% NOTES
%   muiDataSet is used as a superclass to provide data handling 
%   functionality for classes that import different types of data or models
%   that need to store outputs in a consistent and documented format.
%   - inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
%
%   The grid is assigned using setGrid as a dstable held in obj.Data.Grid.
%   obj.Data.Grid is used in getGrid to recover a 'grid' struct with the 
%   following fields: x,y,z,t,irow,desc,metadata,ishead,xM,cline.
%
%   Model grids include RunParam whereas imported grids do not. A model
%   generated Grid includes copies of the classes used as input to the model.
%
%   When using reGridData, any 'New' grid definition must use an
%   Input class with a classname containing the string 'GridData',
%   'GridProps', or 'GridParams', eg similar to GD_GridProps.
% SEE ALSO
%   see ChannelForm model for example of usage
%   GD_GridProps.m, GD_ImportData.m
%
% Author: Ian Townend
% CoastalSEA (c)Jan 2022
%--------------------------------------------------------------------------
%
    methods
        function obj = setGrid(obj,griddata,dims,meta)
            %load the grid results into a dstable 
            % obj - any GDinterface subclass
            % gridddata - cell array with matrix of elevations, z and any
            % other data to be included in the dstable
            % dims - structure with dimensions of griddata
            %        x and y vector dimensions of z - must be unique
            %        t datetimes, durations or values to assign to rows in dstable
            %        cline struct of centre line co-ordinates 
            % meta - stuct with source and data fields to describe model
            %        used and additional run details, respectively
            if nargin<4
                meta.source = metaclass(obj).Name; %use classname as default
                meta.data = '';  
            end                
            if (~isfield(dims,'t') || isempty(dims.t)), dims.t = years(0); end
            
            %create dstable for grid data
            zdsp = GDinterface.setDSproperties();
            if ~isduration(dims.t)
                zdsp = dsproperties(zdsp,'Grid properties'); 
                editDSproperty(zdsp,'Row');
            end
            zdst = dstable(griddata{:},'RowNames',dims.t,'DSproperties',zdsp); 
            zdst.Dimensions.X = dims.x; %row or column vector because
            zdst.Dimensions.Y = dims.y; %dstable holds Dimensions as column vectors   
       
            %assign metadata about model
            zdst.Source = meta.source;  %source file or model class
            zdst.MetaData = meta.data;  %details any manipulation eg type of model or grid rotation
            
            %struct for xy coordinate of curvilinear transformation
            if isfield(dims,'cline')
                zdst.UserData.cline = dims.cline;
            else
                nrec = length(dims.t);
                zdst.UserData.cline(nrec,1) = struct('x',[],'y',[],...
                                            'maxiter',[],'tolerance',[]);
            end
            
            obj.Data.Grid = zdst;             
        end
%%
function grid = getGrid(obj,irow,promptxt,dsetname)
            %retrieve a grid and the index for the case and row
            % obj - any GDinterface subclass that contains a Grid dstable
            % irow - row index of the dstable to extract grid from (optional)
            % promptxt - user input prompt (optional)
            % dsetname - default is 'Grid'. Option to enable mutliple instances
            %           for a given case with a name 'GridN' wher N is numeric
            if nargin<2
                irow = [];
                promptxt = 'Select grid:';
                dsetname = 'Grid';
            elseif nargin<3
                promptxt = 'Select grid:';
                dsetname = 'Grid';
            elseif nargin<4
                dsetname = 'Grid';
            end

            if isempty(promptxt)
                promptxt = 'Select grid:';
            end
            
            dst = obj.Data.(dsetname);
            if height(dst)>1 && isempty(irow)
                %propmpt user to select timestep
                list = dst.DataTable.Properties.RowNames;
                irow = listdlg('PromptString',promptxt,...
                               'Name','Grid selection','SelectionMode','single',...
                               'ListSize',[200,200],'ListString',list);
                if isempty(irow), grid = []; return; end
            elseif isempty(irow)
                irow = 1;
            end
            
            %extract grid for selected row  
            grid.z = squeeze(dst.Z(irow,:,:));           
            grid.x = dst.Dimensions.X;
            grid.y = dst.Dimensions.Y;
            grid.t = dst.RowNames(irow); %text for selected row
 
            %add other index and description
            grid.irow = irow;
            grid.desc = dst.Description;
            grid.metadata = dst.MetaData;
            %struct for xy coordinate of curvilinear transformation
            grid.cline = dst.UserData.cline(irow);             
        end 
%%
        function addGrid(obj,muicat,newgrid,timestep,dims,sourcetxt,ismsg)            
            %add a grid to an existing Case table
            % obj - any GDinterface subclass that contains a Grid dstable
            % muicat - handle to muiCatalogue
            % newgrid - grid to add to Form table
            % timestep - value to use for RowName
            % dims - empty (not used) included to be consistent with FDGintereface
            % sourcetxt - file name or model name to identify data source
            % ismsg - true to use defualt message, false to suppress message
            dst = obj.Data.Grid;      %selected dstable
            %add as a new row to existing table
            nrec = height(dst)+1;
            dst.DataTable{nrec,1} = newgrid;
            dst.RowNames(nrec) = timestep;
            %assign metadata about data
            dst.Source{nrec} = sourcetxt; 
            %struct for xy coordinate of curvilinear transformation
            if isfield(dims,'cline')
                dst.UserData.cline(nrec) = dims.cline;
            else
                dst.UserData.cline(nrec,1) = struct('x',[],'y',[],...
                                            'maxiter',[],'tolerance',[]);
            end
            %sort in row order and load case
            dst = sortrows(dst);
            obj.Data.Grid = dst;  
            classrec = classRec(muicat,caseRec(muicat,obj.CaseIndex));
            updateCase(muicat,obj,classrec,ismsg);
        end
%%
        function deleteGrid(obj,classrec,catrec,muicat)
            %delete a grid from an existing Case table (ie a row)
            % obj - any GDinterface subclass that contains a Grid dstable
            % classrec - class record number
            % catrec - catalogue record for instance
            % muicat - handle to muiCatalogue
            dst = obj.Data.Grid;      %selected dstable

            delist = dst.DataTable.Properties.RowNames; %char rownames
            %select variable to use
            promptxt = {sprintf('Select Row(s) to delete')}; 
            row2use = 1;
            if length(delist)>1
                [row2use,ok] = listdlg('PromptString',promptxt,...
                                 'Name','Delete','SelectionMode','multiple',...
                                 'ListSize',[250,100],'ListString',delist);
                if ok<1, return; end  
            end
            promptxt = sprintf('Delete row: %s\n',delist{row2use});
            selopt = questdlg(promptxt,'Delete row',...
                                      'Yes','No','No');
            if strcmp(selopt,'No'), return; end

            dst.DataTable(row2use,:) = [];  %delete selected variable
            dst.Source(row2use) = [];
            dst.UserData.cline(row2use) = [];
            
            obj.Data.Grid = dst;   
            msgtxt = sprintf('Grid data deleted from: %s',catrec.CaseDescription);
            
            %check whether there are any property tables that need updating
            [obj,isok] = deleteProperties(obj,row2use);
            if isok
                msgtxt = sprintf('Grid data and Form Properties deleted from: %s',catrec.CaseDescription);
            end

            %update the modified Case record
            updateCase(muicat,obj,classrec,false); %suppress message 
            getdialog(msgtxt);
        end             
%% ------------------------------------------------------------------------
% Methods to manipulate grids used in Grid Tools menu:
% translateGrid, rotateGrid, reGridData, subGridData, exportGrid,
% addFormProperties, setSubGrid
%--------------------------------------------------------------------------
        function translateGrid(obj,muicat)
            %interactively translate grid
            %apply to all grids in dataset and updating existing table.
            % obj - any GDinterface subclass that contains a Form dstable
            % classrec - class record of selected instance (can be obtained
            % using selectCaseObj. 
            % muicat - handle to muiCatalogue
            formdst = copy(obj.Data.Grid);
            oldgrid.x = formdst.Dimensions.X;
            oldgrid.y = formdst.Dimensions.Y;
            oldgrid.z = squeeze(formdst.Z(1,:,:)); %first row/time
            
            grid = setGridOrigin(obj,oldgrid);
            if isempty(grid), return; end
            formdst.Dimensions.X = grid.x;
            formdst.Dimensions.Y = grid.y;
  
            %save as a new case or update existing case
            addORupdate(obj,muicat,formdst);            
        end
%%
        function rotateGrid(obj,muicat)
            %interactively flip or rotate grid to make selection and then
            %apply to all grids in dataset and updating existing table.
            % obj - any GDinterface subclass that contains a Form dstable
            % classrec - class record of selected instance (can be obtained
            % using selectCaseObj). 
            % muicat - handle to muiCatalogue
            formdst = copy(obj.Data.Grid);
            oldgrid.x = formdst.Dimensions.X;
            oldgrid.y = formdst.Dimensions.Y;
            oldgrid.z = squeeze(formdst.Z(1,:,:)); %first row/time
            
            [~,rotate] = orientGrid(obj,oldgrid);
            if rotate==0, return; end
            nrec = height(formdst);
            for i=1:nrec
                grid = getGrid(obj,i);
                for j=1:length(rotate)
                    grid = rotateGridData(obj,grid,rotate(j));
                end
                sz = num2cell(size(grid.z));
                newgrids(i,:,:) = reshape(grid.z,1,sz{:}); %#ok<AGROW>
            end
            formdst.DataTable.Z = newgrids;
            formdst.Dimensions.X = grid.x;
            formdst.Dimensions.Y = grid.y;
            
            %assign metadata about data
            if rotate>0
                rotxt = {'Fliped LR','Flipped UD','Rotated +90 deg',...
                               'Rotated -90 deg','Invert X','Inveret Y'};
                orientxt = rotxt{rotate(1)};
                for j=2:length(rotate)
                    orientxt = sprintf('%s; %s',orientxt,rotxt{rotate(j)});
                end
                if iscell(formdst.Source)
                    formdst.Source{1} = sprintf('%s grid from: %s',orientxt,formdst.Source{1});
                else
                    formdst.Source = sprintf('%s grid from: %s',orientxt,formdst.Source);
                end
            end

            %save as a new case or update existing case
            addORupdate(obj,muicat,formdst); 
        end     
%%
        function reGridData(obj,mobj,gridclasses)
            %regrid a gridded dataset to match another grid or to user
            %specified dimensions. Because the grid changes size need to
            %apply to all time steps and save as a new record
            % obj - any GDinterface subclass that contains a Form dstable
            % mobj - mui model instance     
            % gridclasses - clases to be used for case selection
            muicat = mobj.Cases;   %handle to muiCatalogue
            muiInp = mobj.Inputs;  %handle to inpuy classes used in mui
            
            %get selected grid Dimensions
            formdst = copy(obj.Data.Grid);
            casedesc = formdst.Description;
            x = formdst.Dimensions.X;
            y = formdst.Dimensions.Y;
            
            %prompt user to define source of X-Y axes
            answer = questdlg('Use Grid Parameter dimensions, or an existing Cartesian Grid?',...
                'Regrid','New','Existing','New');
            if strcmp(answer,'New')
                %use Grid Parameter definitions to
                if isempty(muiInp)
                    %no Input classes defined
                    warndlg('Set Grid Dimensions in order to create a New grid');
                    return
                else
                    inpclasses = fieldnames(muiInp);
                    idx = contains(inpclasses,{'GridData','GridProps',...
                        'GridParams'},'IgnoreCase',true);
                    if isempty(idx)                %no grid definition class found
                        warndlg('Grid definition class not found in reGridData');
                        return
                    elseif any(idx) && sum(idx)>1  %trap multiple classes
                        warndlg('Multiple grid definition classes found in reGridData');
                        return
                    elseif any(idx) && sum(idx)==1 %single grid definition class found
                        gridclass = inpclasses{idx};
                        grdobj = muiInp.(gridclass);
                        [newgrid.x,newgrid.y] = getGridDimensions(grdobj);
                        metatxt = sprintf('Regridded from %s using run parameter definition',casedesc);
                    else                           %Input classes do not include a grid definition class
                        warndlg('Grid definition class not found in reGridData');
                        return
                    end
                end
            else
                %select an existing Cartesian grid to use
                promptxt = 'Select a Model Grid to define grid X and Y ranges:';
                refobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                if isempty(refobj), return; end
                
                refdst = refobj.Data.Grid;
                griddesc = refdst.Description;
                metatxt = sprintf('Regridded from %s using %s',casedesc,griddesc);
                newgrid.x = refdst.Dimensions.X;
                newgrid.y = refdst.Dimensions.Y;
            end
            
            %create new grid
            [X,Y] = meshgrid(newgrid.x,newgrid.y);
            Nx = length(newgrid.x); Ny = length(newgrid.y);
            for i=1:height(formdst)
                z =squeeze(formdst.DataTable.Z(i,:,:));
                Z = griddata(y,x,z,Y,X,'linear');            
                impdata(i,:,:) = reshape(Z',1,Nx,Ny); %#ok<AGROW>
            end
            formdst.DataTable.Z = impdata;     %update grid data
            formdst.Dimensions.X = newgrid.x;
            formdst.Dimensions.Y = newgrid.y;
            formdst.Source = metatxt;          %assign metadata about data

            %save as a new case or update existing case
            addORupdate(obj,muicat,formdst); 
        end
%%
        function subGridData(obj,muicat)
            %interactively define a subgrid and save grid as a new Case
            % obj - any GDinterface subclass that contains a Form dstable
            % muicat - handle to muiCatalogue
            formdst = copy(obj.Data.Grid);
            casedesc = formdst.Description; 
            grid = getGrid(obj,1); %first row/time  
            %create subgrid selection plot
            [subdomain,sublimitxt] = gd_subdomain(grid);
            if isempty(subdomain), return; end
            
            %extract selected grid
            [subgrid,ixo,iyo] = gd_subgrid(grid,subdomain); %muifunction
            xo = subgrid.x; yo = subgrid.y;
            Nx = length(xo); Ny = length(yo);
            for i=1:height(formdst)
                zi =squeeze(formdst.DataTable.Z(i,:,:));
                zo = zi(min(ixo):max(ixo),min(iyo):max(iyo));
                impdata(i,:,:) = reshape(zo,1,Nx,Ny); %#ok<AGROW>                
                if isfield(formdst.UserData,'xM')
                    %adjust distance to mouth, if changed
                    xM = formdst.UserData.xM(i,1);
                    delX0 = min(grid.x)-min(xo);
                    formdst.UserData.xM(i,1) = xM+delX0;
                end
            end 
            %extract defined subgrid for each row in table
            formdst.DataTable.Z = impdata;
            formdst.Dimensions.X = xo;
            formdst.Dimensions.Y = yo;
            formdst.Source = sprintf('Subgrid of %s using %s',...
                                                casedesc,sublimitxt);
                            
            %save as a new case or update existing case
            addORupdate(obj,muicat,formdst); 
        end
%%
        function displayGridDims(obj)
            %get dimensions of a grid.
            % obj - any GDinterface subclass that contains a Form dstable
            % classrec - class record of selected instance (can be obtained
            % using selectCaseObj. 
            % muicat - handle to muiCatalogue
            grid = getGrid(obj);
            format long
            gdims = gd_dimensions(grid);
            display(gdims);    %display in long forat to the Command window
            format("default")  %restore default (tablefigure uses default)
            figtitle = 'Grid Dimensions';
            headtxt = sprintf('Grid dimensions for %s',grid.desc);
            tablefigure(figtitle,headtxt,gdims);            
        end
%%       
        function addSurface(obj,muicat)
            %add horizontal surface to an extisting grid (eg surrounding
            %land level)
            grid = getGrid(obj,1); %first row/time  
            zmax = max(grid.z,[],'all');
            promptxt = {'Level of cut surface (mAD)','Level to be added (or NaN)'};                                   
            defaults = {num2str(zmax),num2str(NaN)};
            answer = inputdlg(promptxt,'Add surface',1,defaults);
            if isempty(answer), return; end
            zcut = str2double(answer{1});
            zlevel = str2double(answer{2});

            formdst = obj.Data.Grid;
            zi = zeros(size(formdst.Z));
            for i=1:height(formdst)
                grid = getGrid(obj,i);
                idz = (grid.z>zcut) | (isnan(grid.z));
                grid.z(idz) = zlevel;
                zi(i,:,:) = grid.z;
            end
            %save resulting grid
            xynew = struct('zi',zi,'cline',formdst.UserData.cline);
            %save as a new case or update existing case
            addORupdate(obj,muicat,xynew); 
        end
%%
        function infillSurface(obj,muicat)
            %add horizontal surface to an extisting grid (eg surrounding
            %land level)
            methods = {'0 - default plate method',...
                       '1 - least squares plate method',...
                       '2 - linear system of nan elements',...
                       '3 - As 0 with higher del^4 opeartor',...
                       '4 - spring method between neighbours',...
                       '5 - 8 nearest neighbours (not rec)'};
            answer = listdlg('ListString',methods,'SelectionMode','single',...
                      'ListSize',[200,100],'PromptString','Select method');
            if isempty(answer), return; end
            method = answer-1;
            
            formdst = obj.Data.Grid;
            zi = zeros(size(formdst.Z));
            for i=1:height(formdst)
                grid = getGrid(obj,i);
                newz = inpaint_nans(grid.z,method);
                zi(i,:,:) = newz;
            end
            %save resulting grid
            xynew = struct('zi',zi,'cline',formdst.UserData.cline);
            %save as a new case or update existing case
            addORupdate(obj,muicat,xynew); 
        end
%%
        function curvilinear_xy2sn(obj,muicat)
            %map grid from cartesian to curvilinear coordinates
            cline = obj.Data.Grid.UserData.cline;
            if ~isempty(cline.tolerance)
                %a cline definition exists including tolerance settings
                %check whether to reset and get a new defitiion                
                answer1 = questdlg('Using existing centre-line or import new?',...
                               'Centrel-line','Existing','New','New');
                if strcmp(answer1,'New')    
                    %reset tolerance to force new inputs below and in
                    %gd_xy2sn
                    cline.tolerance = []; cline.maxiter = [];
                end
            end
            %
            if isempty(cline.tolerance)
                [fname,path,~] = getfiles('MultiSelect','off',...
                    'FileType','*.txt;*.csv','PromptText','Select centre-line file:');
                if fname==0, return; end %user cancelled

                nhead = 1;
                filename = [path,fname];
                data = readinputfile(filename,nhead);  
                if isempty(data), return; end
                cline.x = data{1};
                cline.y = data{2};
            end
            %obj.Data.Grid.UserData.cline = struct('x',[],'y',[]);
            %add offsets to cline to account for any offset at the mouth
            msg1 = 'To avoid slithers being created when the grid is transformed';
            msg2 = 'an offset can be added to adjust the centre line position';
            msg3 = 'in the x-direction relative to the mouth';
            answer = questdlg(sprintf('%s\n%s\n%s\nIs an offset required?',msg1,msg2,msg3),...
                               'Origin','Yes','No','Yes');
            if strcmp(answer,'Yes')
                grid = getGrid(obj,1);
                ixM = floor(grid.xM/abs(grid.x(2)-grid.x(1)))+1;
                %offset of meander from mouth to avoid skew altering mouth alingment
                offset = inputdlg('Offset distance along x-axis from mouth (m)?','xy2sn',1,{'0'});
                if isempty(offset)
                    nx0=0; 
                else
                    delx = abs(cline.x(2)-cline.x(1)); 
                    nx0 = round(str2double(offset)/delx);
                end
                ioffset = ixM+nx0;
                cline.y = [ones(ioffset,1)*cline.y(1);cline.y(1:end-ioffset)];
            end
            
            %add meander to each grid in table
            formdst = obj.Data.Grid;
            zi = zeros(size(formdst.Z));
            for i=1:height(formdst)
                grid = getGrid(obj,i);
                [grid,cline] = gd_xy2sn(grid,cline,true,true);
                if isempty(grid), return; end
                zi(i,:,:) = grid.z;
            end
            %save resulting grid
            
            xy2sn = struct('zi',zi,'cline',cline);
            addORupdate(obj,muicat,xy2sn);
        end
%%
        function curvilinear_sn2xy(obj,muicat)
            %map grid from curvilinear to cartesian coordinates
            
            %remove meander from each grid in table
            formdst = obj.Data.Grid;
            zi = zeros(size(formdst.Z));
            for i=1:height(formdst)
                grid = getGrid(obj,i);
                cline = grid.cline;
                if isempty(cline) || isempty(cline.x), continue; end
                grid = gd_sn2xy(grid,cline,true);
                zi(i,:,:) = grid.z;
            end
            %save resulting grid
            sn2xy = struct('zi',zi,'cline',formdst.UserData.cline);
            addORupdate(obj,muicat,sn2xy);
        end
%%
        function exportGrid(obj)
            %prompt user to select a Case and export grid as xyz tuples
            grid = getGrid(obj);
            if isempty(grid), return; end %user did not select a row case
            desc = [grid.desc,'_',char(grid.t)];             
            desc = matlab.lang.makeValidName(desc);
            fname = sprintf('%s%s.txt',filesep,desc);

            [Y,X] = meshgrid(grid.y,grid.x);
            N = numel(X);
            xyz = [reshape(X,[N,1]),reshape(Y,[N,1]),reshape(grid.z,[N,1])];

            % Save as a text file            
            header = '%f %f %f';
            xyztxt = sprintf('%f\t%f\t%f\n',xyz');
            path = pwd;
            fid = fopen([path,fname],'w');
            fprintf(fid,'%s\n', header);
            fprintf(fid,'%s\n', xyztxt);
            fclose(fid);
            
            getdialog(sprintf('Grid written to %s',desc));
        end

%%
function [grid,orient] = orientGrid(obj,grid0)
            %Allow user to select whether to flip the grid (180 deg) or
            %rotate by 90 deg.
            % rotate - 0 = no change, 1 = flip lr, 2 = flip ud, 
            %          3 = +90, 4 = -90, 5 = invert x, 6 = invert y
            % orient - vector of rotate values in order applied
            figtitle = sprintf('Orient data');
            promptxt = 'Flip or rotate grid?';
            tag = 'PlotFig'; %used for collective deletes of a group
            butnames = {'Accept','Flip LR','Flip UD','+90','-90','invX','invY','Reset'};
            position = [0,0.03,1,0.93];
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames,position,0.8);
            %plotGrid(obj,h_plt,grid0.x,grid0.y,grid0.z');
            gd_plotgrid(h_plt,grid0);
            axis equal tight %assume geographical projection or grid of similar dimensions
            ok = 0;
            orient = [];
            grid = grid0;
            while ok<1
                waitfor(h_but,'Tag')
                if ~ishandle(h_but)   %this handles the user deleting figure window
                    grid = []; orient = 0; return; %continue with no rotation                    
                elseif strcmp(h_but.Tag,'Flip LR')
                    rotate = 1;   
                elseif strcmp(h_but.Tag,'Flip UD')
                    rotate = 2;     
                elseif strcmp(h_but.Tag,'+90')
                    rotate = 3;
                elseif strcmp(h_but.Tag,'-90')
                    rotate = 4;
                elseif strcmp(h_but.Tag,'invX')
                    grid.x = flipud(grid.x);
                    rotate = 5;
                elseif strcmp(h_but.Tag,'invY')    
                    grid.y = flipud(grid.y);
                    rotate = 6;
                elseif strcmp(h_but.Tag,'Reset')
                    %plotGrid(obj,h_plt,grid0.x,grid0.y,grid0.z');
                    gd_plotgrid(h_plt,grid0);
                    grid = grid0; 
                    orient = [];
                    h_but.Tag = 'reset';
                    continue;
                else                %user accepts selection
                    ok = 1;
                    continue;
                end
                
                if ~any(strcmp(h_but.Tag,{'invX','invY'}))
                    grid = rotateGridData(obj,grid,rotate);
                    orient = [orient,rotate]; %#ok<AGROW> 
                else
                    orient = [orient,rotate]; %#ok<AGROW> 
                end
                gd_plotgrid(h_plt,grid);
                h_but.Tag = 'reset';
            end
            delete(h_plt.Parent);
        end
    end
%% ------------------------------------------------------------------------
% Static methods to maniputlate grids: setCombinedGrids, gridMenuOptions 
%--------------------------------------------------------------------------
    methods (Static)
        function setCombinedGrids(muicat,gridclasses,ismax,prompts)     
            %superimpose one grid on another
            % obj - any GDinterface subclass that contains a Form dstable
            % muicat - handle to muiCatalogue
            % gridclasses - clases to be used for case selection
            % ismax - true: use maximum value at each point, 
            %         false: use minimum (optional)
            % prompts - struct, with fields 'first' and 'second' to define
            %           prompts for selection of two cases to be combined (optional)
            if nargin<3                
                prompts.first = 'Select a Case to use as the ''Master'' grid';
                prompts.second = 'Select a Case to use as the ''Base'' grid';
                ismax = [];
            elseif nargin<4
                prompts.first = 'Select a Case to use as the ''Master'' grid';
                prompts.second = 'Select a Case to use as the ''Base'' grid';
            end
            
            %get first case
            obj1 = selectCaseObj(muicat,[],gridclasses,prompts.first);
            if isempty(obj1), return; end
            dst1 = copy(obj1.Data.Grid);
            grid1 = getGrid(obj1,1);

            %get second case
            obj2 = selectCaseObj(muicat,[],gridclasses,prompts.second);
            if isempty(obj2), return; end
            grid2 = getGrid(obj2);

            if isempty(ismax)
                %prompt user if the selection condition has not been given
                answer = questdlg('Use maximum or minimum values at each point in grid?',...
                                  'Add grids','Max','Min','Max');
                ismax = true;
                if strcmp(answer,'Min'), ismax = false; end
            end

            for i=1:height(dst1)
                grid1.z =squeeze(dst1.DataTable.Z(i,:,:));
                dst1.DataTable.Z(i,:,:) = GDinterface.combineGrids(grid1,...
                                                            grid2,ismax);
            end
            
            %update data range to capture combined grid range
            activatedynamicprops(dst1,{'Z'}); %calls updateVarNames which resets range
            
            %assign metadata about data
            dst1.Source = sprintf('Comined grids using %s and %s',...
                                                   grid1.desc,grid2.desc);  
            %save as a new case
            setGridObj(obj1,muicat,dst1);  
        end
%%
        function diffGridsPlot(muicat,gridclasses,prompts)     
            %generate a plot of the difference between two grids
            % obj - any GDinterface subclass that contains a Form dstable
            % muicat - handle to muiCatalogue
            % gridclasses - clases to be used for case selection
            % prompts - struct, with file 'first' and 'second' to define
            %           prompts for selection of two cases to be combined (optional)
            if nargin<3               
                prompts.first = 'Select Case for z1 (diff is z2-z1)';
                prompts.second = 'Select Case for z2 (diff is z2-z1)';
            end
            %get first case
            obj1 = selectCaseObj(muicat,[],gridclasses,prompts.first);
            if isempty(obj1), return; end
            grid1 = getGrid(obj1);

            %get second case
            obj2 = selectCaseObj(muicat,[],gridclasses,prompts.second);
            if isempty(obj2), return; end
            grid2 = getGrid(obj2);
            
            %check that grids have the same dimensions and if not select
            if numel(grid2.z)~=numel(grid1.z)
                [X1,Y1] = meshgrid(grid1.x,grid1.y);
                [X2,Y2] = meshgrid(grid2.x,grid2.y);    
                answer = questdlg('Grid dimensions differ, select grid to interpolate onto',...
                    'Grid difference',grid1.desc,grid2.desc,'Quit',grid1.desc);
                if strcmp(answer,'Quit')
                    return;
                elseif strcmp(answer,grid1.desc)
                    plotx = grid1.x;
                    ploty = grid1.y;
                    gridinterp = interp2(X2,Y2,grid2.z',X1,Y1,'linear');
                    zdiff = gridinterp'-grid1.z;
                else
                    plotx = grid2.x;
                    ploty = grid2.y;
                    gridinterp = interp2(X1,Y1,grid1.z',X2,Y2,'linear');
                    zdiff = grid2.z-gridinterp';
                end
            else
                plotx = grid1.x;
                ploty = grid1.y;
                zdiff = grid2.z-grid1.z;
            end
            
            ztest = zdiff(~isnan(zdiff)); %need test to exclude NaNs
            if all(ztest==ztest(1,1),'all')
                %contourf does not render constant surfaces
                msgbox(sprintf('Constant difference of %.3f',ztest(1,1)),'Difference Plot')
            else
                figure('Name','Grid Plot','Units','normalized','Tag','PlotFig');
                zlist = GDinterface.getContours(zdiff);                
                if isempty(zlist)
                    contourf(plotx,ploty,zdiff','EdgeColor', 'none');
                    h3 = colorbar;
                else                %use contour intervals if defined
                    contourf(plotx,ploty,zdiff',zlist,'EdgeColor', 'none');
                    h3 = colorbar('Ticks',zlist);
                    clim([zlist(1),zlist(end)]);
                end                 
                h3.Label.String = 'Change in elevation (m)';
                xlabel('Distance along x-axis (m)'); 
                ylabel('Distance along y-axis (m)'); 
                title('Difference plot')
                subtitle(sprintf('%s @ %s - %s @ %s',grid1.desc,grid1.t,...
                                                    grid2.desc,grid2.t))
            end
        end
%%
        function plotSections(muicat,gridclasses)
            %display grid and allow user to interactively define start and
            %end points of a line and plot section along line in a figure
            promptxt = {'Select Case to use (Cancel to quit):','Select timestep:'};
            [obj,~,irec] = selectCaseDatasetRow(muicat,[],...
                                                         gridclasses,promptxt,1);
            if isempty(obj) || isempty(irec)
                return;
            elseif ~isfield(obj.Data,'Grid')
                warndlg('No grid for selected Case'); return; 
            end
            grid = getGrid(obj,irec);   %grid for selected year
            pmtxt = sprintf('Define sections and plot\nSelect menu option');
            PL_PlotSections.Figure(grid,pmtxt,true);      %function plots selected sections  
        end
%%
        function points = getGridLine(muicat,gridclasses,issave)
            %interactively digitise a line and save to a file
            if nargin<3
                issave = true;
            end
            promptxt = {'Select Case to use (Cancel to quit):','Select timestep:'};
            [obj,~,irec] = selectCaseDatasetRow(muicat,[],...
                                                         gridclasses,promptxt,1);
            if isempty(obj) || isempty(irec)
                return;
            elseif ~isfield(obj.Data,'Grid')
                warndlg('No grid for selected Case'); return; 
            end
            
            %check whether z values are also to be captured
            isxyz = false;
            promptz = 'Capture z values?';
            zanswer = questdlg(promptz,'Digitize Z','Yes','No','No');
            if strcmp(zanswer,'Yes'), isxyz = true; end
            
            %desc = sprintf('%s at %s',obj.Data.Grid.Description,char(obj.Data.Grid.RowNames(irec)));
            grid = getGrid(obj,irec);   %grid for selected year
            %digitise points and return xy struct (outype=2)
            points = gd_digitisepoints(grid,{'Digitise Lines'},2,isxyz,false);
            if isempty(points), return; end

            % Save as a text file
            if issave && ~isempty(points)
                path = pwd;
                promptxt = {sprintf('Save points to file or cancel\nPath'),'File name'};
                defaults = {path,'Grid line'};
                answers = inputdlg(promptxt,'DigiLine',1,defaults);
                if isempty(answers), return; end
                xy = cell2mat(struct2cell(points'));
                fid = fopen([answers{1},'\',answers{2},'.txt'],'w');
                if isxyz  %z has been digitised as well as x,y
                    fprintf(fid,'%s\n', '%f %f %f');
                    fprintf(fid,'%f\t%f\t%f\n', xy);               
                else
                    fprintf(fid,'%s\n', '%f %f');
                    fprintf(fid,'%f\t%f\n', xy);
                end
                fclose(fid);
            end
        end

%%
        function gridImage(muicat,gridclasses)
            %generate a tiff image of a grid
             promptxt = {'Select Case to use (Cancel to quit):','Select timestep:'};
            [obj,classrec,irec] = selectCaseDatasetRow(muicat,[],...
                                                         gridclasses,promptxt,1);
            if isempty(obj) || isempty(irec)
                return;
            elseif ~isfield(obj.Data,'Grid')
                warndlg('No grid for selected Case'); return; 
            end

            grid = getGrid(obj,irec);   %grid for selected year
            xLim = minmax(grid.x);  %limits of grid data
            yLim = minmax(grid.y);
            hf = figure('Tag','PlotFig','Visible','off');
            ax = axes(hf);
            Z = grid.z';
            C = pcolor(ax,Z);
            shading interp
            axis equal tight            
            xlim(xLim);  ylim(yLim);
            cmap = gd_colormap([min(grid.z,[],'all'),max(grid.z,[],'all')]);
            cLimits = clim;
            im = struct('XData',xLim,'YData',yLim,'CData',C.CData,'CMap',cmap,'CLim',cLimits);
            
            % Use the whos function to get information about the array
            info = whos('Z');            
            if info.bytes>1e6  %extract the size in bytes (Z==C.CData)
                %option to resize the image
                answer = questdlg('Resize image?','Image','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    inp = inputdlg('Enter scaling factor (sf x dimensions)','Image',1,{'1'});
                    if ~isempty(inp) && ~strcmp(inp{1},'1')
                        factor = str2double(inp{1});
                        im.CData = imresize(C.CData,factor);
                    end
                end
            end
            ax = axes(hf);
            imagesc(ax,'XData',im.XData,'YData',im.YData,'CData',im.CData); 
            axis equal tight
            frame = getframe(ax);
            delete(hf)
            
            %check plot by overlaying image (set hf to Visible and comment
            %delete hf to activate)
            % hold(ax,'on')
            % h_im = imagesc('XData',im.XData,'YData',im.YData,'CData',im.CData);
            % h_im.AlphaData = 0.5;
            % hold(ax,'off')

            %to use this method the calling class must include a get.formatypes method
            %(see EDBimport for eg)
            answer = questdlg('Save to Case or File?','GridImage',...
                                                    'Case','File','Case');
            if strcmp(answer,'Case')
                datasetname = 'GeoImage';
                obj.DataFormats{2} = obj.formatypes{datasetname,1}{1};
                [dst,ok] = callFileFormatFcn(obj,'setData',obj,im,grid.desc);
                if ok<1 || isempty(dst), return; end          
                obj.Data.(datasetname) = dst.(datasetname);
                updateCase(muicat,obj,classrec,false);
                getdialog(sprintf('GeoImage added to Case: %s',grid.desc),[],1)   
            else
                desc = matlab.lang.makeValidName(grid.desc);
                answer = questdlg('Save as .mat or .tif file?','GridImage',...
                                                 'mat','tif','txt','txt');                
                filename = sprintf('GeoImage_%s.%s',desc,answer);
                inp = inputdlg('Filename:','GeoImage',1,{filename});
                if isempty(inp); return; end
                if strcmp(answer,'mat')
                    save(inp{1},'im',"-mat");
                elseif strcmp(answer,'tif')
                    gd_write_tiff(filename,frame.cdata,im);  
                else
                    gd_write_image(filename,desc,im)
                end
                getdialog(sprintf('Image written to %s',filename));
            end            
        end
%%
        function gridMenuOptions(mobj,src,gridclasses)
            %default set of menu options for use in Model UIs
            % mobj - mui model instance
            % src - handle to calling menu option
            % gridclasses - cell array list of classes to use in case selection
            %
            % Sub-menu for using Grid Tools:
            % menu.Setup(X).List = {'Translate Grid','Rotate Grid',...
            %                       'Re-Grid','Sub-Grid',...
            %                       'Combine Grids','Add Surface',...
            %                       'To curvilinear','From curvilinear',...                   
            %                       'Display Dimensions','Difference Plot',...
            %                       'Plot Sections','Digitise Line',...
            %                       'Export xyz Grid','User Function'};                                                                        
            % menu.Setup(X).Callback = repmat({@obj.gridMenuOptions},[1,14]);
            % menu.Setup(X).Separator = [repmat({'off'},[1,6]),...
            %                 {'on','off','on','off','off','on','on','on'}];         
            % Call using any class that inherits GDinterface:
            % e.g. in CF_FormModel call gridMenuOptions(obj,src,gridclasses);            
            muicat = mobj.Cases;   %handle to muiCatalogue
            switch src.Text 
                case 'Translate Grid'
                    promptxt = 'Select a Case to Translate (move origin):'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    translateGrid(cobj,muicat); 
                case 'Rotate Grid'
                    promptxt = 'Select a Case to Flip or Rotate:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    rotateGrid(cobj,muicat);    
                case 'Re-Grid'
                    promptxt = 'Select a Case to Regrid:'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    reGridData(cobj,mobj,gridclasses);
                case 'Sub-Grid'
                    promptxt = 'Select a Case to extract a Subgrid:'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    subGridData(cobj,muicat);                
                case 'Combine Grids'
                    GDinterface.setCombinedGrids(muicat,gridclasses);
                case 'Add Surface'
                    promptxt = 'Select a Case to add surface:'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    addSurface(cobj,muicat); 
                case 'Infill Surface'
                    promptxt = 'Select a Case to infill surface:'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    infillSurface(cobj,muicat); 
                case 'To curvilinear'
                    promptxt = 'Select a Case to convert TO a curvilinear grid:'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    curvilinear_xy2sn(cobj,muicat); 
                case 'From curvilinear'
                    promptxt = 'Select a Case to convert FROM a curvilinear grid:'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    curvilinear_sn2xy(cobj,muicat); 
                case 'Difference Plot'
                    GDinterface.diffGridsPlot(muicat,gridclasses);
                case 'Plot Sections'
                    GDinterface.plotSections(muicat,gridclasses);
                case 'Grid Image'
                    GDinterface.gridImage(muicat,gridclasses);
                case 'Digitise Line'
                    GDinterface.getGridLine(muicat,gridclasses);
                case 'Export xyz Grid'
                    promptxt = 'Select a Case to Export as xyz tuples:'; 
                    [cobj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    exportGrid(cobj);  
                case 'Display Dimensions'
                    promptxt = 'Select a Case to display grid dimensions'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    displayGridDims(cobj);
                case 'User Function'
                    promptxt = 'Select a Case for User Function'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    gd_user_function(cobj,mobj);
            end
        end
    end
%% ------------------------------------------------------------------------
% Protected methods to maniputlate grids: orientGrid,
% rotateGridData, getGridDimensions, setDSproperties
%--------------------------------------------------------------------------
    methods (Access=protected)
        function grid = setGridOrigin(~,grid0)
            %Allow user to define a new origin for the grid
            grid = grid0;
            figtitle = sprintf('Translate grid');
            promptxt = 'Translate grid?';
            tag = 'PlotFig'; %used for collective deletes of a group
            butnames = {'Accept','Edit','invX','invY','Quit'};
            position = [0.2,0.4,0.4,0.4];
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames,position,0.8);
            %plotGrid(obj,h_plt,grid0.x,grid0.y,grid0.z');
            gd_plotgrid(h_plt,grid0);
            ok = 0;
            while ok<1
                waitfor(h_but,'Tag')
                if ~ishandle(h_but) %this handles the user deleting figure window 
                    grid = []; return;             
                elseif strcmp(h_but.Tag,'Edit') 
                    grid = GDinterface.getOrigin(grid);
                    %h_but.Tag = 'reset';
                    %gd_plotgrid(h_plt,grid);
                elseif strcmp(h_but.Tag,'invX')
                    grid.x = flipud(grid.x);
                elseif strcmp(h_but.Tag,'invY')    
                    grid.y = flipud(grid.y);    
                elseif strcmp(h_but.Tag,'Quit') 
                    grid = []; 
                    delete(h_plt.Parent);
                    return; 
                else   %user accepted
                    ok = 1; 
                end
                gd_plotgrid(h_plt,grid);
                h_but.Tag = 'reset';
            end
            delete(h_plt.Parent);
            diffX = grid.x-grid0.x;
            diffY = grid.y-grid0.y;
            if all(diffX==0) && all(diffY==0)
                %no change so return empty grid
                grid = [];
            end
        end
%%
        function grid = rotateGridData(~,grid,rotate)
           %format grid to load into a dstable as single row per grid
           %saves data as z array [1,Nx,Ny] for each time step
           %X and Y assumed fixed and saved as dimensions
           % rotate - 0 = no change, 1 = flip, 2 = +90, 3 = -90
           Z = grid.z;
           if rotate>0
               temp = grid.y;
               if rotate==1           %flip the grid to vertical mirror 
                   Z = flipud(Z);
               elseif rotate==2       %flip the grid to horizontal mirror
                   Z = fliplr(Z);
               elseif rotate==3       %rotate the grid (+90)
                   Z = flipud(Z)';    
                   grid.y = grid.x; grid.x = temp;
               elseif rotate==4       %rotate the grid (-90)
                   Z = flipud(Z');
                   grid.y = grid.x; grid.x = temp;
               end
           end
           grid.z = Z;
           %check input plot
           % figure; contourf(grid.x,grid.y,Z');         
        end
%%
        function contouredSurfacePlot(~,ax,xi,yi,zi,zlimits,casedesc)
            %generate surface plot with contours and land/sea colormap
            hs = surfc(ax,xi,yi,zi,'FaceColor','interp','EdgeColor', 'none');
            hContour = hs(2); % get the handle to the contour lines
            hContour.ContourZLevel = zlimits(1); % set the contour's Z position (default: hAxes.ZLim(1)=-10)
            gd_colormap(zlimits);
            c = colorbar;
            c.Label.String = 'Elevation (mAD)';
            view(315,30);             
            hold on
            contour3(ax,xi,yi,zi,'LineColor',[0.8,0.8,0.8]);
            hold off  
            gd_ax_dir(ax,xi,yi);  %ensure grid and axes are correctly oriented
            xlabel('Length (km)'); 
            ylabel('Width (km)');  
            zlabel('Elevation (mAD)');
            title(casedesc,'FontWeight','normal','FontSize',10);
        end      
%% ------------------------------------------------------------------------
% Protected methods to Update existing gridsd that have been modified by
% one of the Grid Tools (rotate, sub-grid, etc)
%--------------------------------------------------------------------------               
        function addORupdate(obj,muicat,newdata)
            %prompt user to add or update existing record and save grid
            % newdata is either the dstable to be used, or a struct of the 
            % grid z and cline values (used by curvilinear functions).
            
            answer = questdlg('Add new case, or update/add to existing case?','Save grid','Add','Update','Quit','Add');
            
            if strcmp(answer,'Quit')
                 getdialog('Modified grid NOT saved');
                return;
            elseif strcmp(answer,'Add')
                %create new record
                if isa(newdata,'dstable')
                    fdst = newdata;
                else
                    fdst = copy(obj.Data.Grid);
                    fdst.Z(1,:,:) = newdata.zi;
                    fdst.UserData.cline = newdata.cline;
                end
                setGridObj(obj,muicat,fdst); %copies form property table to new instance
            else            
                %overwrite exisiting grid data set with new grid
                if isa(newdata,'dstable')
                    if isfield(obj.Data,'Grid') && ~isempty(obj.Data.Grid)
                        answer = questdlg('Add additional grid to case or Replace existing grid?','Save grid','Add','Replace','Add');
                        if strcmp(answer,'Add')
                            datasetnames = fieldnames(obj.Data);
                            nrec = sum(contains(datasetnames,'Grid'))+1;
                            datasetname = sprintf('Grid%d',nrec);
                            obj.Data.(datasetname) = newdata;
                        else
                            obj.Data.Grid = newdata;
                        end
                    end
                else
                    obj.Data.Grid.Z(:,:,:) = newdata.zi;  
                    obj.Data.Grid.UserData.cline = newdata.cline;
                end
                classrec = classRec(muicat,caseRec(muicat,obj.CaseIndex));
                updateCase(muicat,obj,classrec,true);
            end 
        end
%%
        function setGridObj(obj,muicat,formdst)
            %assign metadata and properties to new instance of object that
            %has been modified. 
            % obj - any GDinterface subclass with unmodified class instance
            % muicat - handle to muiCatalogue 
            % formdst - GDinterface format dstable for updated Grid
            %Note: grid manipulation functions update all grids in Grid and
            %formdst is the Grid dstable for all grids
            %          
            %This function may need to be overloaded if the grid derives 
            %from a model and RunParam and/or Form properties need updating   
            if isempty(obj.DataFormats)
                heq = str2func(metaclass(obj).Name);
                newobj = heq();                  %new instance of class object
            else
                formatfile = obj.DataFormats{2};
                heq = str2func(sprintf('@(file) %s(file)',metaclass(obj).Name));
                newobj = heq(formatfile);                  %new instance of class object
            end

            caserec = caseRec(muicat,obj.CaseIndex);
            casedef = getRecord(muicat,caserec);
            
            %assign to default struct for GDinterface data type
            newobj.Data.Grid = formdst;   %assign new grid(s)

            setCase(muicat,newobj,casedef.CaseType);
        end
    end
%% ------------------------------------------------------------------------
% Static, private methods to: getOrigin, getContours, combineGrids
%--------------------------------------------------------------------------
    methods (Static, Access=protected)
        function grid = getOrigin(grid)
            %Prompt user to define new grid origin
            promptxt = {'Change in X-origin','Change in Y-origin'};
            defaults = {'0','0'};
            title = 'Define origin';
            answer = inputdlg(promptxt,title,1,defaults);
            if ~isempty(answer)
                grid.x = str2double(answer{1})+grid.x;
                grid.y = str2double(answer{2})+grid.y;
            end
        end 
%%
        function zlist = getContours(zvar)
            %prompt user for max, min and z-interval to use in contour plot
            zmin = min(zvar,[],'all');
            zmax = max(zvar,[],'all');
            zint = (zmax-zmin)/10;
            promptxt = {'Minimum z contour','Maximum z contour','contour interval'};
            title = 'set contours';
            defaults = {num2str(zmin);num2str(zmax);num2str(zint)};
            answer = inputdlg(promptxt,title,1,defaults);
            if isempty(answer)
                zlist = [];
            else
                zmin = str2double(answer{1});
                zmax = str2double(answer{2});
                zint = str2double(answer{3});
                zlist = zmin:zint:zmax;
            end
        end
%%
        function new_z = combineGrids(grid1,grid2,ismax)
            %combine two grids based on the minimum or maximum values in
            %each grid
            x1 = grid1.x;   
            if min(x1)<0, x1 = max(x1)-x1; end
            
            %ensure that the two grids are aligned
            [X,Y] = ndgrid(x1,grid1.y);
            z2 = griddata(grid2.y,grid2.x,grid2.z,Y,X); 
            %create new grid
            if ismax
                new_z = max(grid1.z,z2);
            else
                new_z = min(grid1.z,z2);
            end
            [m,n] = size(new_z);
            new_z = reshape(new_z,1,m,n);  %dimenions for row in dstable
        end  

%% ------------------------------------------------------------------------
% Protected methods to set DSproperties
%-------------------------------------------------------------------------- 
        function dsp = setDSproperties()
            %define the metadata properties for the demo data set
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            % *note: these would be better in the gd_property functions so
            % *that the defintion is in the same file as the assignment
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...
                            'Name',{'Z'},...                  
                            'Description',{'Elevation'},...
                            'Unit',{'mAD'},...
                            'Label',{'Elevation (mAD)'},...
                            'QCflag',{'model'});   
            dsp.Row = struct(...
                            'Name',{'Time'},...
                            'Description',{'Time'},...
                            'Unit',{'y'},...
                            'Label',{'Time (yr)'},...
                            'Format',{'y'});       
            dsp.Dimensions = struct(...    
                            'Name',{'X','Y'},...
                            'Description',{'X-axis','Y-axis'},...
                            'Unit',{'m','m'},...
                            'Label',{'X-axis (m)','Y-axis (m)'},...
                            'Format',{'',''});   
        end       
    end 
end