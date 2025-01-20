classdef (Abstract = true) FGDinterface < GDinterface
%
%-------class help------------------------------------------------------
% NAME
%   FGDinterface.m
% PURPOSE
%   Abstract class providing additional methods for use with gridded data sets
%   (eg imported or model data) to extend the functionality of the muiDataSet
%   and GDintereface abstract classes to handle channel form properties
% NOTES
%   GDinterface  is used as a superclass to provide grid handling which, 
%   in turn, uses muiDataSetto provide data handling functionality for
%   classes that import different types of data or models that need to 
%   store outputs in a consistent and documented format.
%   - inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
%
%   The grid is assigned using setGrid as a dstable held in obj.Data.Grid.
%   obj.Data.Grid is used in getGrid to recover a 'grid' struct with the 
%   following fields: x,y,z,t,irow,desc,metadata,ishead,xM,cline.
%   Classes inheriting FGDinterface can also have additonal tables, held in
%   obj.Data, including: Hypsometry, SectionProps, GrossProps, Plan and 
%   WaterLevels.
%
%   Model grids include RunParam whereas imported grids do not. A model
%   generated Grid includes copies of the classes used as input to the model.
%   For example in ChannelForm this includes at least some of CF_HydroData,
%   CF_SediData,GD_GridProps,CF_TransData,RunProperties,WaterLevels.
%
%   The plan form held for a model uses the specified form (exp, power)
%   whereas the properties are derived from a a grid for imported datasets.
%   If Add Properties is called for either this overwrites the plan form
%   with values derived from the grid. Similarly models that use
%   addProperties will have values based on the grid not the model output.
%
% SEE ALSO
%   see ChannelForm model for example of usage
%   GDinterface, GD_GridProps.m, FGD_ImportData.m
%
% Author: Ian Townend
% CoastalSEA (c)Jan 2022
%--------------------------------------------------------------------------
%    
    methods
        function obj = setGrid(obj,griddata,dims,meta)
            %load the grid results into a dstable 
            % obj - any FGDinterface subclass
            % gridddata - cell array with matrix of elevations, z and any
            % other data to be included in the dstable
            % dims - structure with dimensions of griddata
            %        x and y vector dimensions of z - must be unique
            %        t datetimes, durations or values to assign to rows in dstable
            %        cline struct of centre line co-ordinates 
            %     and FGDinterface additions:
            %        xM distance to the mouth of an inlet/estuary
            %        ishead orientation of x-axis relative to mouth
            %        Lt distance from mouth to tidal limit
            %        Rv river properties hydraulic depth, width and CSA
            % meta - stuct with source and data fields to describe model
            %        used and additional run details, respectively
            % NB: if xM and cline are not arrays of same as t (no of grids)
            % they are padded with the first value provided
            
            %load grid by calling superclass function in GDinterface
            obj = setGrid@GDinterface(obj,griddata,dims,meta);
            zdst = obj.Data.Grid;
            
            %add properties that are specific to a channel form            
            %orientation of x-axis relative to mouth
            zdst.UserData.ishead = dims.ishead;
            %distance to mouth from origin of grid (min(x))
            nrec = length(dims.t);
            if length(dims.xM)==nrec
                zdst.UserData.xM = dims.xM; 
            else
                zdst.UserData.xM(1:nrec,1) = dims.xM(1);
            end
            
            %distance from mouth to tidal limit
            if length(dims.Lt)==nrec
                zdst.UserData.Lt = dims.Lt;
            else
                zdst.UserData.Lt(1:nrec,1) = dims.Lt(1);
            end
            
            %river properties Rv.Hr, Rv.Wr and Rv.Ar
            if length(dims.Rv)==nrec %struct array containing
                zdst.UserData.Rv = dims.Rv;
            else
                zdst.UserData.Rv(1:nrec,1) = dims.Rv(1);
            end
            
            obj.Data.Grid = zdst;
        end
%%
        function grid = getGrid(obj,irow,promptxt)
            %retrieve a grid and the index for the case and row
            % obj - any FGDinterface subclass that contains a Grid dstable
            % irow - row index of the dstable to extract grid from (optional)
            if nargin<2
                irow = [];
                promptxt = 'Select grid:';
            elseif nargin<3
                promptxt = 'Select grid:';
            end
            %get grid by calling superclass function in GDinterface
            grid = getGrid@GDinterface(obj,irow,promptxt);
            if isempty(grid), return; end
            dst = obj.Data.Grid;
            %add properties that are specific to a channel form
            grid.ishead = dst.UserData.ishead; %orientation of x-axis true if head
            grid.xM = dst.UserData.xM(grid.irow);   %distance to mouth from grid origin
            grid.Lt = dst.UserData.Lt(grid.irow);   %distance from mouth to tidal limit
            grid.Rv = dst.UserData.Rv(grid.irow);%river properties
        end 
%%
        function addGrid(obj,muicat,newgrid,timestep,dims,sourcetxt,ismsg)
            %add a grid to an existing Case table
            % obj - any FGDinterface subclass that contains a Grid dstable
            % muicat - handle to muiCatalogue
            % newgrid - grid to add to Form table
            % timestep - value to use for RowName
            % dims - struct containig xM, Lt and river properties 
            % sourcetxt - file name or model name to identify data source
            % ismsg - true to use defualt message, false to suppress message
            addGrid@GDinterface(obj,muicat,newgrid,timestep,dims,sourcetxt,false)
            dst = obj.Data.Grid;      %selected dstable     
            nrec = find(dst.RowNames==timestep); %addGrid sorts rows
            %properties added by FGDinterface
            dst.UserData.xM(nrec) = dims.xM;
            dst.UserData.Lt(nrec) = dims.Lt;
            dst.UserData.Rv(nrec) = dims.Rv;
            
            %load case
            obj.Data.Grid = dst;  
            classrec = classRec(muicat,caseRec(muicat,obj.CaseIndex));
            updateCase(muicat,obj,classrec,ismsg);
        end
%%
        function deleteGrid(obj,classrec,catrec,muicat)
            %delete a grid from an existing Case table (ie a row)
            % obj - any FGDinterface subclass that contains a Grid dstable
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

            dst.DataTable(row2use,:) = [];  %delete selected  grid
            dst.Source(row2use) = [];
            dst.UserData.cline(row2use) = [];
            %remove propoerties added by FGDinterface
            %NB: ishead is not deleted because it applies to all grids
            dst.UserData.xM(row2use) = [];
            dst.UserData.Lt(row2use) = [];
            dst.UserData.Rv(row2use) = [];
                        
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
% Methods to set grid properties: setProperties,addProperties,delProperties
%--------------------------------------------------------------------------
        function obj = setProperties(obj,zwl,Wz,limits,histint)
            %initialise the tables for derived properties of gridded data set
            % obj - any class that contains a Grid dstable   
            % zwl - struct or table, of water levels at zhw,zmt,zlw
            % Wz - struct or table of width plan form data at 3 elevations
            % limits - upper and lower limits for vertical range:
            %          0=use grid; 1=prompt user; [x,y]=limits to use           
            % histint - vertical interval for hypsometry  - optional user
            %           prompted if not defined
            if nargin<5 || isempty(histint) || histint==0
                histint = setHypsometryInterval(obj);
            end
            
            %to initialise tables use first row of Grid table, irow=1
            [pt,mtxt] = getPropTables(obj,1,zwl,Wz,limits,histint);
            if isempty(pt), return; end  %wls or plan areas not set
            
            %load the water levels into obj.Data class instance
            obj = setPropTable(obj,pt,'WaterLevels',mtxt.wl);
            %load the plan form widths into obj.Data class instance
            obj = setPropTable(obj,pt,'Plan',mtxt.pl);  
            %load the derived hypsometry properties into obj.Data class instance           
            obj = setPropTable(obj,pt,'Hypsometry',mtxt.hy);          
            %load the derived prism-area properties into obj.Data class instance
            obj = setPropTable(obj,pt,'SectionProps',mtxt.sp);
            %load the derived gross properties into obj.Data class instance
            obj = setPropTable(obj,pt,'GrossProps',mtxt.gp);
        end
%%
        function obj = addProperties(obj,irow,zwl,Wz,limits,histint)
            %add properties to an existing set of grid property tables
            % obj - any class that contains a Grid dstable   
            % irow - row of Grid table to use
            % zwl - struct or table, of water levels at zhw,zmt,zlw
            % Wz - struct or table of width plan form data at 3 elevations
            % limits - upper and lower limits for vertical range:
            %          0=use grid; 1=prompt user; [x,y]=limits to use
            % meta - stuct with source and data fields to describe model
            %        used and additional run details, respectively. Only
            %        needed to initialise tables (ie first grid)
            % histint - vertical interval for hypsometry  - optional user
            %           prompted if not defined. Ignored if property tables
            %           already exist and value for 1st grid used
            if ~isfield(obj.Data,'Hypsometry') || ...
                              isempty(obj.Data.Hypsometry) || irow==1
                %no properties defined,uses first row of Grid to create tables
                %NB: if irow=1 and tables exist this overwrites property tables
                obj = setProperties(obj,zwl,Wz,limits,histint);
                return;
            else
                %add irow grid properties to the tables, using the same
                %z-values for each grid. This is defined when setting
                %the first grid.
                histint = obj.Data.Hypsometry.UserData.histint;
                limits = obj.Data.Hypsometry.UserData.limits;  
                [pt,~] = getPropTables(obj,irow,zwl,Wz,limits,histint);
                if isempty(pt), return; end  %wls or plan areas not set
                
                %add the water levels into obj.Data class instance
                obj.Data.WaterLevels = [obj.Data.WaterLevels;pt.WaterLevels];
                %add the plan form widths into obj.Data class instance
                obj.Data.Plan = [obj.Data.Plan;pt.Plan];  
                %add the derived hypsometry properties into obj.Data class instance           
                obj.Data.Hypsometry = [obj.Data.Hypsometry;pt.Hypsometry];          
                %add the derived prism-area properties into obj.Data class instance
                obj.Data.SectionProps = [obj.Data.SectionProps;pt.SectionProps];
                %add the derived gross properties into obj.Data class instance
                obj.Data.GrossProps = [obj.Data.GrossProps;pt.GrossProps];
            end
        end
%%
        function [obj,isok] = delProperties(obj,row2use)
            %delete row of properties from existing grid property tables
            % obj - any FGDinterface subclass with Grid Property tables
            % row2use - row(s) to be deleted from ALL tables
            isok = false;
            fnames = fieldnames(obj.Data);
            fnames(strcmp(fnames,'Grid')) = [];
            ntables = length(fnames);
            if ntables>1
                for i=1:ntables
                    %delete the selected rows from the property tables
                    obj.Data.(fnames{i}).DataTable(row2use,:) = [];
                end
                isok = true;
            end    
        end    
%% ------------------------------------------------------------------------
% Methods to update grid properties: addFormProps, delFormProps,
% setModelFormProps and editGridInletData
%--------------------------------------------------------------------------
        function addFormProps(obj,muicat)
            %add form properties to a gridded data set (needs to be a
            %suitable dtm, ie inlet, estuary, valley)     
            % obj - any FGDinterface subclass that contains Property tabes
            % muicat - handle to muiCatalogue
            dst = obj.Data;
            if isfield(dst,'Hypsometry') && ~isempty(dst.Hypsometry) 
                %one or more tables exists and is not empty
                fnames = fieldnames(dst);
                casedesc = obj.Data.Grid.Description;
                questxt = sprintf('%s ',fnames{:});
                questxt = sprintf('%s\nTables: %s\nalready exist. Overwrite?',...
                                                        casedesc,questxt);
                answer = questdlg(questxt,'Add Form Properties',...
                                            'Overwrite','Quit','Quit');
                if strcmp(answer,'Quit'), return; end
            end
            %check whether to update the water levels for each grid
            allgrids = questdlg('Use the same water levels for all grids?',...
                                'Water levels','Yes','No','Yes');
            if strcmp(allgrids,'No'), isnew = true; else, isnew = false; end
            
            %get the limits to use for the hypsometry
            hypslimit = questdlg('Use grid to define limits, or define?',...
                                 'Limits','Use grid','Define','Use grid');
                             
            if strcmp(hypslimit,'Use grid')
                limits = 0;
            else
                limits = setHypsometryLimits(obj);                
            end
            histint = setHypsometryInterval(obj);
            
            %get grid water levels and add first row of property tables
            zwl = getGridWaterLevels(obj,true); %force reset of water levels
            if isempty(zwl), return; end
            grid = getGrid(obj,1);
            Wz = gd_plan_form(grid,zwl);
            if isempty(Wz), return; end           
            obj = setProperties(obj,zwl,Wz,limits,histint);

            %add rows to property tables for each additional grid
            for irow=2:height(dst.Grid)
                if isnew
                    %force reset of water levels if user selects 
                    zwl = getGridWaterLevels(obj,isnew);
                end
                obj = addProperties(obj,irow,zwl,[],limits,histint);
                if height(obj.Data.Hypsometry)~=irow
                    msg = sprintf('Water levels or plan area not defined\nProperty tables not updated');
                    warndlg(msg);
                    return;
                end
            end
            %assign updated dstable
            classrec = classRec(muicat,caseRec(muicat,obj.CaseIndex));
            updateCase(muicat,obj,classrec,true);
        end
%%
        function delFormProps(obj,muicat)
            %delete ALL property tables associated with a selected gridded
            %data set
            % obj - any FGDinterface subclass that contains Property tables
            % muicat - handle to muiCatalogue
            answer = questdlg('Delete ALL property tables?','Delete Properties',...
                                                        'Yes','No','No');
            if strcmp(answer,'No'), return; end
            
            fnames = fieldnames(obj.Data);
            fnames(strcmp(fnames,'Grid')) = [];
            ntables = length(fnames);
            if ntables>0
                for i=1:ntables
                    %delete the selected rows from the property tables
                    obj.Data.(fnames{i}) = [];
                end
            end  
            %assign updated dstable            
            classrec = classRec(muicat,caseRec(muicat,obj.CaseIndex));
            updateCase(muicat,obj,classrec,false);
            msgtxt = sprintf('Form Properties deleted from: %s',obj.Data.Grid.Description);
            getdialog(msgtxt);
        end      
%%
        function obj = setModelFormProps(obj)  %replaces setHyps_SP_GP(obj,meta)
            %add a set of Hypsometry, Section and Gross form properties 
            %function used from models such as CF_TransModel
            grdobj = obj.RunParam.GD_GridProps;
            histint = grdobj.histint;
            allZ = obj.Data.Grid.Z;  %all grids as single array
            limits = [max(allZ,[],'all')+histint,min(allZ,[],'all')-histint];
            
            %add properties to model instance
            for irow=1:height(obj.Data.Grid)
                obj = addProperties(obj,irow,obj.zWL(irow,:),obj.tPlan(irow,:),...
                                                        limits,histint);
                if height(obj.Data.Hypsometry)~=irow
                    msg = sprintf('Water levels or plan area not define\nProperty tables not updated');
                    warndlg(msg);
                    return;
                end
            end
        end
%%
        function editGridInletData(obj)
            %provide option to edit, definition of channel head, x-distance 
            %to mouth and definition of centre-line (if used)
            dst = obj.Data.Grid;

            if height(dst)>1 %select row in table if more than 1
                list = dst.DataTable.Properties.RowNames;
                idx = listdlg('PromptString','Select timestep:',...
                              'Name','Time step','SelectionMode','single',...
                              'ListSize',[200,200],'ListString',list);
                if isempty(idx), return; end %user cancelled
            else
                idx = 1;
            end
            
            defaults{1,1} = char(dst.RowNames(idx));
            defaults{2,1} = num2str(dst.UserData.ishead);
            defaults{3,1} = num2str(dst.UserData.xM(idx));
            
            if isempty(dst.UserData.cline(idx).x)  %no centre-line defined
                defaults{4,1} = '';
                defaults{5,1} = '';
            else
                defaults{4,1} = num2str(dst.UserData.cline(idx).x');
                defaults{5,1} = num2str(dst.UserData.cline(idx).y');
            end
            promptxt = {'Timestep','Minimum X is head of channel (1/0)','X-distance to mouth from grid origin',...
                          'Centre-line X','Centre-line Y'};
            ok=0;
            while ok==0         
                answer = inputdlg(promptxt,'Edit Grid Props',1,defaults);
                if isempty(answer), return; end %user cancelled
                dst.RowNames(idx) = str2duration(answer{1});
                dst.UserData.ishead = logical(str2double(answer{2}));
                dst.UserData.xM(idx) = str2double(answer{3});
                dst.UserData.cline(idx).x = str2num(answer{4})';  %#ok<ST2NM> input of vector
                dst.UserData.cline(idx).y = str2num(answer{5})';  %#ok<ST2NM>
                if length(dst.UserData.cline(idx).x)==length(dst.UserData.cline(idx).y)
                    ok = 1;
                else
                    hw = warndlg('Centre-line must have same number of x and y values');
                    waitfor(hw)
                end
            end
            obj.Data.Grid = dst;
            if ~strcmp(defaults{1},answer{1})
                warndlg('If Timestep changed, the Grid Properties need to be updated')
            end
        end
%%
        function histint = setHypsometryInterval(obj)
            %prompt user to defince the vertical interval for hypsometry            
            if isprop(obj,'RunParam')  && ~isempty(obj.RunParam) && ...
                 isfield(obj.RunParam,'GD_GridProps') && ...
                 ~isempty(obj.RunParam.GD_GridProps.histint)
                defaultval = {num2str(obj.RunParam.GD_GridProps.histint)};
            else
                defaultval = {'0.1'};
            end
            
            prmptxt = {'Interval for hypsometry:'};
            dlgtitle = 'Hypsometry';            
            answer = inputdlg(prmptxt,dlgtitle,1,defaultval);
            if isempty(answer), answer = defaultval; end
            histint = str2double(answer{1});
        end          
%%
        function limits = setHypsometryLimits(obj)
            %prompt user to define the limits to be used for hypsometry estimates
            %uses the range of the grid as default and if user Cancels
            allZ = obj.Data.Grid.Z;
            %this uses max(z) to give the hypsometry of the basin, whereas
            %the model is limited to HW
            limits = [max(allZ,[],'all')+0.1,min(allZ,[],'all')-0.1];
            prmptxt = {'Upper limit for hypsometry:','Lower limit for hypsometry:'};
            dlgtitle = 'Hypsometry';           
            defaultvals = {num2str(limits(1)),num2str(limits(2))};
            answer = inputdlg(prmptxt,dlgtitle,1,defaultvals);
            if ~isempty(answer)
                limits(1) = str2double(answer{1});
                limits(2) = str2double(answer{2});
            end
        end
%%
        function wl = getGridWaterLevels(obj,isnew)
            %get the water levels used in the model, or prompt user to
            %define if none found. isnew used to force new definition of
            %water levels when adding properties to an imported grid. 
            % wl is a struct containing zhw, zmt, zlw
            if nargin<2, isnew = false; end
            dst = obj.Data;
            gridX = dst.Grid.Dimensions.X;
            nx = length(gridX);
            if isprop(obj,'RunParam')  && ~isempty(obj.RunParam)
                %instance of CF_HydroData defining water levels used
                if ~isempty(obj.RunParam.CF_HydroData.zhw) && ...
                                 ~isscalar(obj.RunParam.CF_HydroData.zhw)
                    wlobj = obj.RunParam.CF_HydroData;
                    %wl = struct('zhw',wlobj.zhw,'zmt',wlobj.zmt,'zlw',wlobj.zlw);
                    wl= table(wlobj.zhw,wlobj.zmt,wlobj.zlw,...
                                      'VariableNames',{'zhw','zmt','zlw'});
                elseif isfield(dst,'WaterLevels') && ~isempty(dst.WaterLevels)...
                           && ~isempty(dst.WaterLevels.zhw) && ~isnew                    
                   % wl = table2struct(dst.WaterLevels.DataTable);   
                     wl = dst.WaterLevels.DataTable;  
                else 
                    wl = FGDinterface.setTidalLevels(nx,gridX);
                end
            elseif isfield(dst,'WaterLevels') && ~isempty(dst.WaterLevels)...
                           && ~isempty(dst.WaterLevels.zhw) && ~isnew
                   % wl = table2struct(dst.WaterLevels.DataTable);
                wl = dst.WaterLevels.DataTable;  
            else
                %instance of GD_GridProps defining grid used
                %get the water levels to use for the derived properties
                wl = FGDinterface.setTidalLevels(nx,gridX);
            end
        end
    end
%%
    methods (Access=protected)    
%% ------------------------------------------------------------------------
% Protected methods to create the Property tables, or assign the Property
% tables to a class instance
%--------------------------------------------------------------------------         
       function [PropTables,mtxt] = getPropTables(obj,irow,zwl,Wz,limits,histint)
            %add the derived properties to a gridded data set
            % obj - any class that contains a Grid dstable   
            % irow - row of Grid table to use
            % zwl - struct or table, of water levels at zhw,zmt,zlw
            % Wz - struct or table of width plan form data at 3 elevations
            % limits - upper and lower limits for vertical range:
            %          0=use grid; 1=prompt user; [x,y]=limits to use
            % histint - vertical interval for hypsometry  - optinal user
            %           prompted if not defined
            PropTables = []; mtxt = [];
            
            grid = getGrid(obj,irow);  %struct of x,y,z,t values
            
            %load the water levels into a dstable
            if isempty(zwl)
                zwl = getGridWaterLevels(obj);
                if isempty(zwl), return; end
                mtxt.wl = ': user input';
            else
                mtxt.wl = ': using model output';    
            end
            wldst = FGDinterface.setDStable(zwl,grid,'WaterLevels');
            
            %load the plan form widths into a dstable
            if isempty(Wz)
                Wz = gd_plan_form(grid,zwl);
                if isempty(Wz), return; end
                mtxt.pl = ': derived from grid';
            else
                mtxt.pl = ': using model output';
            end
            pldst = FGDinterface.setDStable(Wz,grid,'Plan');

            %load the derived hypsometry properties into dstable     
            mtxt.hy = ': using gd_basin_hypsometry';
            [hypsdst,histdst] = gd_basin_hypsometry(grid,zwl,histint,limits);

            %load the derived prism-area properties into dstable
            mtxt.sp = ': using gd_section_properties';
            spdst = gd_section_properties(grid,zwl,histdst);

            %load the derived gross properties into dstable
            mtxt.gp = ': using gd_gross_properties';
            gpdst = gd_gross_properties(grid,zwl,spdst);
            
            PropTables = struct('WaterLevels',wldst,'Plan',pldst,...
                                'Hypsometry',hypsdst,'SectionProps',spdst,...
                                'GrossProps',gpdst); 
        end
%%
       function obj = setPropTable(obj,PropTables,tablename,mtxt)
           %add the property table, vdst, to the grid instance
           vdst = PropTables.(tablename);
           vdst.Source = [obj.Data.Grid.Description,mtxt];
           vdst.MetaData = obj.Data.Grid.MetaData;
           vdst.Description = obj.Data.Grid.Description;
           obj.Data.(tablename) = vdst;
       end   
%% ------------------------------------------------------------------------
% Protected methods to Update existing gridsd that have been modified by
% one of the Grid Tools (rotate, sub-grid, etc). Includes updating of
% RunParam and Form properties.
%--------------------------------------------------------------------------     
        function caseid = setGridObj(obj,muicat,formdst)
            %assign metadata and properties to new instance of object that
            %has been modified. 
            % obj - any FGDinterface subclass with unmodified class instance
            % muicat - handle to muiCatalogue 
            % formdst - GDinterface format dstable for updated Grid
            %Note: grid manipulation functions update all grids in Grid and
            %formdst is the Grid dstable for all grids
            %
            %NB: this method overloads GDinterface.setGridObj
            heq = str2func(metaclass(obj).Name);
            newobj = heq();                  %new instance of class object
            caserec = caseRec(muicat,obj.CaseIndex);
            casedef = getRecord(muicat,caserec);
            
            %assign to default struct for GDinterface data type
            newobj.Data.Grid = formdst;   %assign new grid(s)

            oldgrid = getGrid(obj,1);     %first row of old Grid instance
            if ~isempty(obj.RunParam)
                rpinput = fieldnames(obj.RunParam);
                for i=1:length(rpinput)
                    newobj.RunParam.(rpinput{i}) = copy(obj.RunParam.(rpinput{i}));     
                end
            end
            %update the RunParam classes if the grid is from a model
            newobj = updateRunParamClasses(newobj,oldgrid);
            
            %update any Selection properties if used
            if isprop(obj,'Selection')
                newobj.Selection = obj.Selection;
            end
            
            %update any other property tables included in the source table
            fnames = fieldnames(obj.Data);
            fnames(strcmp(fnames,'Grid')) = [];
            ntables = length(fnames);
            if ntables>1
                for i=1:ntables
                    %copy the source tables
                    newobj.Data.(fnames{i}) = copy(obj.Data.(fnames{i}));
                end
                %update the source tables to the dimensions of new grid
                newobj = updateFormProperties(newobj,oldgrid);
            end
            caseid = setCase(muicat,newobj,casedef.CaseType);
        end
%%
        function obj = updateRunParamClasses(obj,oldgrid)
            %update the RunParam classes 
            % obj - any grid class derived from model that uses RunParam
            % oldgrid - original grid definition
            formdst = obj.Data.Grid;
            %if RunParam have been set, update these to the new grid
            if isprop(obj,'RunParam') && ~isempty(obj.RunParam) &&...
                                      isfield(obj.RunParam,'GD_GridProps')
                x = formdst.Dimensions.X;
                y = formdst.Dimensions.Y;
                grdobj = obj.RunParam.GD_GridProps;
                obj.RunParam.GD_GridProps =  setGridDimensions(grdobj,x,y);
                obj.RunParam.CF_HydroData = interpHydroData(obj,oldgrid,true); %true return obj
            end            
        end
%%        
        function newzwl = interpHydroData(obj,oldgrid,isobj,irow)
            %interpolate model water levels along x-axis onto modified grid
            %water levels in CF_HydroData are transient so may not exist
            % isobj - logical flag, true to return obj of type updated,
            %         otherwise returns a table fo water levels zhw,zmt,zlw
            % irow - empty or not used to check obj.RunParam.CF_HydroData
            %        or irow grid to retrieve obj.Data.WaterLevels
            if nargin<4, irow = []; end
            
            if isempty(irow)
                hydobj = obj.RunParam.CF_HydroData;    %current water levels
            elseif ~isfield(obj.Data,'WaterLevels')  && ...
                                 isfield(obj.RunParam,'CF_HydroData')
                hydobj = obj.RunParam.CF_HydroData;    %current water levels 
            else
                hydobj = obj.Data.WaterLevels(irow,:); %saved water levels
            end
            
            if isempty(hydobj.zhw) || isscalar(hydobj.zhw)
               newzwl = hydobj;
               return; 
            end

            newX = obj.Data.Grid.Dimensions.X;
            oldX = oldgrid.x;
            zhw = interp1(oldX,hydobj.zhw,newX,'linear')';
            zmt = interp1(oldX,hydobj.zmt,newX,'linear')';
            zlw = interp1(oldX,hydobj.zlw,newX,'linear')';

            %determine how to assign results based on isobj and obj used
            if isobj
                if isa(hydobj,'CF_HydroData')
                    %populate properties of object
                    hydobj.zhw = zhw; hydobj.zmt = zmt; hydobj.zlw = zlw;
                else
                    %assign to row of existing table  
                    hydobj.DataTable{irow,:} = [zhw,zlw,zmt];
                end
                newzwl = hydobj;
            else
                newzwl = table(zhw,zmt,zlw,'VariableNames',{'zhw','zmt','zlw'});
            end
        end        
%%
        function obj = updateFormProperties(obj,oldgrid)
            %update the addional property tables for Form, Plan,
            %SectionProps, GrossProps and WaterLevels if dimensions are
            %modified eg by rotating, sub-sampling etc
            % obj - any FGDinterface subclass with Props dstables to update
            % oldgrid - original grid defintion
            
            %delete existing table
            % fnames = fieldnames(obj.Data);
            % fnames(strcmp(fnames,'Grid')) = [];
            % ntables = length(fnames);
            % if ntables>0
            %     for i=1:ntables
            %         %delete the selected rows from the property tables
            %         obj.Data.(fnames{i}) = [];
            %     end
            % end  

            %retrieve the vertical interval and limits used for existing
            %Form properties
            histint = obj.Data.Hypsometry.UserData.histint;
            limits = obj.Data.Hypsometry.UserData.limits;
            
            %if WaterLevel properties are being used update to new grid 
            %add rows to property tables for each additional grid
            for irow=1:height(obj.Data.Grid)
                zwl = interpHydroData(obj,oldgrid,false,irow); %false retuns table
                grid = getGrid(obj,irow);
                Wz = gd_plan_form(grid,zwl);
                obj = addProperties(obj,irow,zwl,Wz,limits,histint);
                if height(obj.Data.Hypsometry)~=irow
                    msg = sprintf('Water levels or plan area not defined\nProperty tables not updated');
                    warndlg(msg);
                    return;
                end
            end    
        end

    end   
%    
    methods (Static)    
%% ------------------------------------------------------------------------
% Methods to adjust grids: addValleyBase and addShore
%--------------------------------------------------------------------------
        function addValleyBase(muicat,gridclasses)
            %use an xyz definition of the channel thalweg to create a valley base
            %this is used to modify an existing valley defined above MHW to include
            %the valley form down to the pre-Holocene surface (ie below the
            %existing channel bed)
            promptxt = {'Select Case to use (Cancel to quit):','Select timestep:'};
            [obj,~,irec] = selectCaseDatasetRow(muicat,[],...
                                                         gridclasses,promptxt,1);
            if isempty(obj) || isempty(irec), return; end
            grid = getGrid(obj,irec);   %grid for selected year
            gd = gd_dimensions(grid);
            
            %load the thalweg centre-line data
            
            [fname,path,~] = getfiles('MultiSelect','off',...
                'FileType','*.txt;*.csv','PromptText','Select centre-line file:');
            if fname==0, return; end %user cancelled
            filename = [path,fname];
            data = readinputfile(filename,1);  
            if isempty(data), return; end
            vline.x = data{1}; vline.y = data{2}; vline.z = data{3};
            %cline used to save centre line details in metadata - only 
            %updated with z-values if not already defined because the
            %existing cline may define a conformal mapping. vline used to
            %update valley
            cline = obj.Data.Grid.UserData.cline;
            if isempty(cline.x)
                cline = vline;
            end
            
            %plot selected line on grid and give user option to quit
            gd_lineongrid_plot(grid,vline,'Selected valley base line');
            answer = questdlg('Accept base line and continue?','ValleyBase','Accept','Quit','Accept');
            if strcmp(answer,'Quit'), return; end
            
            %define upper and lower bounds for the channel mask
            prmptxt = {'Upper mask cutoff','Offshore mask width','Valley base width'};
            dlgtitle = 'Valley Base';
            defaults = {'5','1000','100'};
            answer = inputdlg(prmptxt,dlgtitle,1,defaults);
            if isempty(answer)
                return; 
            else
                upcut = str2double(answer{1});      %z value of upper limit
                nin_offw = round(str2double(answer{2})/gd.dely/2); %number of nodes in half-offshore-valley
                valwidth = str2double(answer{3});   %half-width of valley base
                nin_valw = round(valwidth/gd.dely); %number of nodes in half base
            end

            %define an area for the channel from the mouth seawards
            xmouth = grid.xM+min(grid.x);
            [~,idy]= min(abs(vline.x-xmouth));  %x-index of mouth on centre-line
            ymouth = vline.y(idy);              %x-index of mouth on centre-line
            [~,idmx] = min(abs(grid.x-xmouth)); %x-index of mouth on axis
            [~,idmy] = min(abs(grid.y-ymouth)); %y-index of mouth on axis
            idoffshore = false(gd.xint,gd.yint); idcutoff = idoffshore;
            if idmx>1
                idoffshore(1:idmx,idmy-nin_offw:idmy+nin_offw) = true;
            end
            %define the area of the channel landward of the mouth
            idcutoff(idmx:end,:) = grid.z(idmx:end,:)<upcut;  %channel mask
            if size(idcutoff)==size(idoffshore)
                idx = logical(idcutoff+idoffshore);
            else
                warndlg('Offshore mask width too large for grid')
                return;
            end
            
            hd = waitbar(0,'Computing new valley base');
            %distances along thalweg from first point
            sdata = [0;cumsum(hypot(diff(vline.x),diff(vline.y)))]; 
            %interpolate thalweg xy data onto x of grid
            xthal = grid.x;
            ythal = interp1(vline.x,vline.y,xthal,'pchip','extrap');
            %find nearest y of grid to define xy of thalweg
            [X,Y] = meshgrid(grid.x,grid.y);
            N = numel(X);
            XY = [reshape(X,[N,1]),reshape(Y,[N,1])];
            k = dsearchn(XY,[xthal,ythal]);
            ythal = Y(k);
            
            %distances along interpolated thalweg from first point
            sthal = [0;cumsum(hypot(diff(xthal),diff(ythal)))]; 
            
            %start of imported thalweg line and nearest point on the grid
            %may not be the same so add difference to thalweg line length
            offset = hypot((xthal(1)-vline.x(1)),(ythal(1)-vline.y(1)));
            if sdata(1)+offset~=0
                sdata = [0;sdata+offset];      %added 0 is to ensure that base is 
                depths = [vline.z(1);vline.z]; %interpolated to edge of grid
                                               %assumes base line starts near
                                               %edge of grid 
            else
                depths = vline.z;
            end
            %interpolate thalweg elevations along line
            zthal = interp1(sdata,depths,sthal,'linear');
            %replicate the centre-line on the y nodes either side of thalweg
            xthal = repmat(xthal,nin_valw,1);
            zthal = repmat(zthal,nin_valw,1);
            y0 = ythal-valwidth/2-gd.dely;
            ythal = [];
            for i=1:nin_valw
                yi = y0+i*gd.dely;
                %assumes valley approximately follows x-axis
                ythal = [ythal;yi];  %#ok<AGROW>
            end
            %create a grid from scattered data points
            [xl,yl,zl] = gd_grid_line(xthal,ythal,zthal,0); %zero used for missing grid values
            Z = griddata(yl,xl,zl,Y,X,'linear')';           %interpolate onto valley grid
            Z(isnan(Z)) = 0;
            FGDinterface.valleybase_checkplot(grid,X,Y,Z,idx);

            %create new grid and set existing NaN values to maximum z
            newgrid = grid;
            newgrid.z(isnan(newgrid.z)) = gd.zmax; %clean valley grid for interpolation
            %set channel values of grid to zero, add thalweg and then set
            %undefined values in channel to NaN
            newgrid.z(idx) = 0;      %set channel values to zero
            idc = Z~=0;              %mask of thalweg channel
            newgrid.z(idc) = 0;      %remove any grid values in thalweg channel
            newgrid.z = newgrid.z+Z; %add valley & thalweg
            ids = idx & ~idc;        %mask of channel excluding thalweg
            newgrid.z(ids) = NaN;    %points still to be defined
            % figure; contourf(X,Y,newgrid.z'); gd_ax_dir(gca,grid.x); 
            %assume valley runs predominantly along x-axis. iterate along x
            %and interpolate each cross-section
            sz = num2cell(size(newgrid.z));
            zi = zeros(sz{:});
            waitbar(0.5,hd,'Adding new valley base')
            for i=1:gd.xint
                idz = ~isnan(newgrid.z(i,:));
                if length(newgrid.y(idz))<2
                    zi(i,:) = newgrid.z(i,:); %not sure about this!!!!!****************
                else
                    zi(i,:) = interp1(newgrid.y(idz),newgrid.z(i,idz),...
                                                     newgrid.y,'linear');     
                end
                waitbar(0.5+i/2/gd.xint,hd)
            end   
            delete(hd);
            figure; contourf(X,Y,zi'); gd_ax_dir(gca,grid.x);
            %save resulting grid
            xynew = struct('zi',zi,'cline',cline);
            %save resulting grid
            addORupdate(obj,muicat,xynew);            
        end
%%
        function valleybase_checkplot(grid,X,Y,Z,idx)
            figure; 
            subplot(2,1,1)
            contourf(X,Y,Z'); 
            gd_ax_dir(gca,grid.x);%check whether axis direction is reversed  
            title('Valley channel base')
            subplot(2,1,2)
            imagesc(idx')
            title('Mask for valley')
            sgtitle('Check plot')
        end
%%
        function addShore(muicat,gridclasses)
            %add a strip to shore side of a grid by extrapolating from the
            %neighbouring strip. 
            promptxt = 'Select a Case to add Shore:'; 
            [obj,~] = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(obj), return; end
            formdst = copy(obj.Data.Grid);

            grid = getGrid(obj,1);   %grid for first year in dataset
            gd = gd_dimensions(grid);

            promptxt = {'Width of shore to add','Offshore depth (mAD)','Exponent'};
            answer = inputdlg(promptxt,'Shore',1,{'5000','-60','1'});
            if isempty(answer), return; end
            NewShoreWidth = str2double(answer{1});
            OffDepth = str2double(answer{2});
            xponent = str2double(answer{3});

            %setup new x values relative to first x value and account for
            %x-axis direction
            newx = 0:gd.delx:NewShoreWidth-gd.delx;            
            nint = length(newx);
            %now interpolate a shoreline for each grid in record
            nrec = height(formdst);
            
            %because grid changes dimensions when adding shoreline need to
            %adjust all grids in the dataset
            for irec=1:nrec
                grid = getGrid(obj,irec);
                newz = zeros(nint,gd.yint);
                for iy = 1:gd.yint
                    ShoreDepth = grid.z(1,iy);  %current depth on edge of grid
                    depth2width = abs(ShoreDepth-OffDepth)/NewShoreWidth^xponent;
                    newz(:,iy) = OffDepth+(newx.^xponent)*depth2width;                                        
                end  
                [X,Y] = meshgrid(newx,grid.y);
                figure; contourf(X,Y,newz');
                newZ = [newz;grid.z];            %extend the irec grid             
                sz = num2cell(size(newZ));       %reshape grid as a row cell
                newgrids(irec,:,:) = reshape(newZ,1,sz{:});  %#ok<AGROW>
            end          
            formdst.DataTable.Z = newgrids;
            formdst.UserData.xM = grid.xM+(nint)*gd.delx;

            newX = grid.x(1)-gd.xsgn*fliplr(newx+gd.delx);
            formdst.Dimensions.X = [newX';grid.x];
            formdst.Dimensions.Y = grid.y;
            
            %need to update seaward water levels in obj.Data.WaterLevels
            %hence cannot use addORupdate function
            answer = questdlg('Add or update existing?','Add Mods',...
                                  'Add','Update','Cancel','Add');
            if strcmp(answer,'Add')
                %create new record
                caseid = setGridObj(obj,muicat,formdst); %copies form property table to new instance
                cobj = getCase(muicat,caseRec(muicat,caseid));
                cf_offset_wls(cobj,false);  %translate wls, false maintains vector length
                classrec = classRec(muicat,caseRec(muicat,caseid)); 
                updateCase(muicat,cobj,classrec,false); %false=no message
            elseif strcmp(answer,'Update')       
                %overwrite exisitng form data set with new form             
                obj.Data.Grid = formdst;  
                obj = cf_offset_wls(obj,true);  %translate wls, true extends vector
                classrec = classRec(muicat,caseRec(muicat,obj.CaseIndex)); 
                updateCase(muicat,obj,classrec,true);
            end   
        end
%% ------------------------------------------------------------------------
% Method to call Grid Form Tools (see also gridMenuOptions in GDinterface
% for Grid tools)
%--------------------------------------------------------------------------         
        function formMenuOptions(mobj,src,gridclasses)
            %default set of menu options for use in Model UIs
            % mobj - mui model instance
            % src - handle to calling menu option
            % gridclasses - cell array list of classes to use in case selection
            %
            % Sub-menu for using Grid Form Tools:
            % menu.Setup(7).List = {'Add Thalweg to Valley','Extrapolate Shore',...
            %                       'Add Properties','Delete Propeties,...
            %                       'Edit Inlet Definition'};                                                                        
            % menu.Setup(5).Callback = repmat({@obj.gridFormOptions},[1,5]);
            % menu.Setup(5).Separator = {'off','off','on','off','on'}; 
            %
            % Call using any class that inherits FGDinterface:
            % e.g. in CF_FormModel call
            % gridFormOptions(obj,src,gridclasses);            
            muicat = mobj.Cases;   %handle to muiCatalogue
            switch src.Text 
                case 'Add Thalweg to Valley'
                    gridclasses([1,3]) = [];
                    FGDinterface.addValleyBase(muicat,gridclasses);
                case 'Extrapolate Shore'                    
                    FGDinterface.addShore(muicat,gridclasses);
                case 'Add Properties'
                    promptxt = 'Select a Case to Add Form properties'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    addFormProps(cobj,muicat); 
                case 'Delete Properties'
                    promptxt = 'Select a Case to Delete Form properties '; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    delFormProps(cobj,muicat);
                case 'Edit Inlet Definition'
                    promptxt = 'Select a Case to Edit Inlet Data'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    editGridInletData(cobj);
                case 'User Function'
                    promptxt = 'Select a Case for User Function'; 
                    cobj = selectCaseObj(muicat,[],gridclasses,promptxt);
                    if isempty(cobj), return; end
                    gd_user_function(cobj,mobj);
            end  
        end        
    end
%% ------------------------------------------------------------------------
% Static, private methods to: setDStable, setDSproperties, setTidalLevels
%--------------------------------------------------------------------------
    methods (Static, Access=protected)
        function dst = setDStable(var,grid,tablename)     
            %create the dstabel for the var property
            dsp = FGDinterface.setDSproperties(tablename);
            dst = dstable(var,'RowNames',grid.t,'DSproperties',dsp); 
            dst.Dimensions.X = grid.x;
        end
%%
        function dsp = setDSproperties(type)
            %define the metadata properties for the demo data set
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            % *note: these would be better in the gd_property functions so
            % *that the defintion is in the same file as the assignment
            
            %struct entries are cell arrays and can be column or row vectors
            switch type    
               case 'Plan'
                    dsp.Variables = struct(...
                                    'Name',{'Whw','Wmt','Wlw'},...                  
                                    'Description',{'HW Width','MT Width','LW Width'},...                                                
                                    'Unit',{'m','m','m'},...
                                    'Label',{'Width (m)','Width (m)','Width (m)'},...                                            
                                    'QCflag',repmat({'model'},[1,3]));   
                    dsp.Row = struct(...
                                    'Name',{'Time'},...
                                    'Description',{'Time'},...
                                    'Unit',{'y'},...
                                    'Label',{'Time (yr)'},...
                                    'Format',{'y'});       
                    dsp.Dimensions = struct(...    
                                    'Name',{'X'},...
                                    'Description',{'Distance'},...
                                    'Unit',{'m'},...
                                    'Label',{'Distance (m)'},...
                                    'Format',{''});   
                                
                case 'WaterLevels'
                    dsp.Variables = struct(...
                                    'Name',{'zhw','zmt','zlw'},...                  
                                    'Description',{'High water level',...
                                            'Mean water level','Low water level'},...
                                    'Unit',{'mAD','mAD','mAD'},...
                                    'Label',{'Water level (mAD)',...
                                            'Water level (mAD)','Water level (mAD)'},...
                                    'QCflag',repmat({'model'},[1,3]));   
                    dsp.Row = struct(...
                                    'Name',{'Time'},...
                                    'Description',{'Time'},...
                                    'Unit',{'y'},...
                                    'Label',{'Time (yr)'},...
                                    'Format',{'y'});       
                    dsp.Dimensions = struct(...    
                                    'Name',{'X'},...
                                    'Description',{'Distance'},...
                                    'Unit',{'m'},...
                                    'Label',{'Distance (m)'},...
                                    'Format',{''});          
            end
        end
%%
        function zwl = setTidalLevels(nx,gridX)
            %load file of water levels that maps onto x-axis or prompt user
            %to define horiontal elevations for water levels
            % nx - number of points used for x-axis co-ordinates (optional)
            if nargin<2, gridX = []; end
            answer = questdlg('Load levels from file or define?','Water levels',...
                                'Load file','Define','Define');
            if strcmp(answer,'Load file') 
                %user provides an input file of values
                %format: row 1 - header; row2-N: 4 columns
                %of numerical values for distance, zhw, zmt, zlw                
                zwl = FGDinterface.readTidalLevels(gridX);
            else
                %user inputs single values using UI
                zwl = FGDinterface.setTidalLevels_UI();
                if isempty(zwl)
                    return;
                elseif nx>1
                    zwl.zhw = ones(1,nx)*zwl.zhw;
                    zwl.zmt = ones(1,nx)*zwl.zmt;
                    zwl.zlw = ones(1,nx)*zwl.zlw;
                end
            end
            zwl = struct2table(zwl);
        end
%%
        function wl = setTidalLevels_UI()
            %prompt user for levels to used to compute volumes and areas at high
            %and low water
            prmptxt = {'High water level:','Mean tide level:',...
                       'Low water level:'};
            dlgtitle = 'Water levels';
            defaultvals = {'1','0','-1'};
            answer = inputdlg(prmptxt,dlgtitle,1,defaultvals);
            if isempty(answer), wl = []; return; end
            
            wl.zhw = str2double(answer{1});
            wl.zmt = str2double(answer{2});
            wl.zlw = str2double(answer{3});
        end  
%%
        function wl = readTidalLevels(gridX)
            %read tide levels from a file using the dimensions consistent with
            %the x-axis of the grid. Format: row 1 - %f %f %f %f; row2-N: 4 columns 
            %of numerical values for distance, zhw, zmt, zlw
            wl = [];
            [fname,path] = getfiles('PromptText', 'Select water level file',...
                                                        'FileType', '*.txt');
            dataSpec = '%f %f %f %f';
            data = readinputfile([path,'/',fname],1,dataSpec);%get the name of the file and read it
            if isempty(data), return; end  %user aborted or file not read

            if data{1}(1)>gridX(1) || data{1}(end)<gridX(end)
                %data does not cover full extent of grid
                wanrdlg('Water level data does not cover full extnet of x-axis')
                return;
            end
            %interpolate water levels from file onto model grid
            wl.zhw = interp1(data{1},data{2},gridX','linear')';
            wl.zmt = interp1(data{1},data{3},gridX','linear')';
            wl.zlw = interp1(data{1},data{4},gridX','linear')';
        end   
    end
end