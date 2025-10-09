classdef ctBeachProfileData < muiDataSet
%
%-------class help---------------------------------------------------------
% NAME
%   ctBeachProfileData.m
% PURPOSE
%   Class to hold beach profile data
% USAGE
%   obj = ctBeachProfileData() 
% SEE ALSO
%   inherits muiDataSet and uses dstable and dscatalogue
%   format files used to load data of varying formats (variables and file format)
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%   
    properties  
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        % importing data requires muiDataSet propertiesm DataFormats and
        % FileSpec to be defined in class constructor.
        %Additional properties:  
    end
     
    methods  
        function obj = ctBeachProfileData
            %constructor to initialise object
            
            %initialise list of available input file formats. Format is:
            %{'label 1','function name 1';'label 2','function name 2'; etc}
            obj.DataFormats = {'Channel Coastal Observatory UK','profiles_cco_format';...
                               'Chainage profiles [site id date Ch z flag]','profiles_chainage_format'};          
            %define file specification, format is: {multiselect,file extension types}
            obj.FileSpec = {'on','*.txt;*.csv'};                       
        end
    end
%%
%--------------------------------------------------------------------------
%   functions to read data from file and load as a dstable
%   data formats are defined in the format files
%--------------------------------------------------------------------------
    methods  (Static)
        function obj = loadData(muicat,classname,newflg)
            %load user data set from one or more files
            % mobj - handle to modelui instance 
            % classname - name of class being loaded
            % newflag - in call from addData to indicate use of existing
            %           record
            if nargin<3
                newflg = 'Yes';
            end
            
            obj = ctBeachProfileData; 
            %get file format from class definition        
            idf = getFileFormatID(obj,muicat);  %existing uses of format
            if strcmp(newflg,'Yes') || isempty(idf) || length(idf)>1
                %no format set or multiple formats set, so get selection
                ok = setFileFormatID(obj);
                if ok<1, return; end
            else
                obj.idFormat = idf;             %use existing format
            end
            idf = obj.idFormat;  %local value to reset additional instances               

            [fname,path,nfiles] = getfiles('MultiSelect',obj.FileSpec{1},...
                'FileType',obj.FileSpec{2},'PromptText','Select file(s):');
            if ~iscell(fname)
                fname = {fname};   %single select returns char
            end
         
            funcname = 'getData';
            hw = waitbar(0, 'Loading data. Please wait');
            %load file and create master collection which can
            %have multiple profiles (ie locations saved as id_rec)            
            for jf=1:nfiles
                filename = [path fname{jf}];
                [newprofiles,ok] = callFileFormatFcn(obj,funcname,obj,filename);
                if ok<1 || isempty(newprofiles), continue; end
                
                %newprofiles is a struct of dstables with profile_id as the fieldname
                profid = fieldnames(newprofiles);
                idx = strcmp(muicat.Catalogue.CaseClass,classname);
                existprofs = muicat.Catalogue.CaseDescription(idx);
                %loop round and add each profile as a new record
                for ip=1:length(profid)
                    aprof = newprofiles.(profid{ip}); %dstable of profile data
                    %if the profile does not exist save as new record                    
                    if isempty(existprofs) || ...
                                        all(~strcmp(existprofs,profid{ip}))
                        %cumulative list of files names used to load data
                        aprof.Source{jf} = filename;             
                        %add first data set to empty classobj, or
                        %add new profile to classobj                        
                        setDataSetRecord(obj,muicat,aprof,'data',profid(ip),true);
                        obj = ctBeachProfileData;
                        obj.idFormat = idf;    
                    else
                        %profile exists - add data to existing record
                        classrec = find(strcmp(existprofs,profid{ip})); 
                        localObj = muicat.DataSets.(classname)(classrec);                        
                        localObj = addBPdataFile(localObj,aprof);                       
                        updateCase(muicat,localObj,classrec,false);
                    end  
                end
                clear existprofs newprofiles aprof profid
                waitbar(jf/nfiles)
            end

            close(hw);
            if nfiles>1
                getdialog(sprintf('Data loaded in class: %s',classname)); 
            end
        end
%%
        function obj = addData(muicat,classname)
            %add additional data to an existing user data timeseries
            % obj - returns class handle
            newflg = 'No'; %indicates that format selection already defined
            obj = ctBeachProfileData.loadData(muicat,classname,newflg);
        end
   end
%%
    methods             
        function tabPlot(obj,src)
            %generate plot for display on Plot tab (default timeseries plot)
            %data is retrieved by GUIinterface.getTabData            
            %get data for variable
            dataset = getDataSetName(obj);
            dst = obj.Data.(dataset);
            time = dst.RowNames;
           
            ht = findobj(src,'Type','axes');
            delete(ht); 
            ax = axes('Parent',src,'Tag','ProfilePlot'); 
            hold on       
            for j=1:length(time)
                elevation = dst.Elevation(j,:);
                chainage = dst.Chainage(j,:);                     
                plot(ax,chainage,elevation)                
            end  
            h1 = plot(ax,xlim, 0*[1 1],'-.k','LineWidth',0.2);
            % Exclude line from legend
            set(get(get(h1,'Annotation'),'LegendInformation'),...
                                                'IconDisplayStyle','off');   
            hold off
            if length(time)<21
                legend(datestr(time,'dd-mmm-yyyy'))
            end
            sitename = split(dst.Description,'__');
            if length(sitename)>1
                title(sprintf('Profile No: %s-%s',sitename{1},sitename{2}));
            else
                title(sprintf('Profile No: %s',sitename{1}))
            end
            ax.XLabel.String = 'Chainage (m)';
            ax.YLabel.String = 'Elevation (mOD)';
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot 
        end 
%%
        function getProfileLocationsPlot(obj,muicat)
            %get data to plot the locations of all the selected profiles
            %prompt user to select profiles to be used (multiple profiles)
            [~,caserec] = getProfileSelection(obj,muicat);
            %
            answer = questdlg('Use existing base figure or new figure?',...
                               'Profiles','New','Existing','New');
            if strcmp(answer,'New')
                hf = figure('Name','Profile Locations','Tag','PlotFig', ...
                            'NumberTitle','off','Units','normalized');                    
                ax = axes(hf);
            else
                selectedFig = select_figure({'PlotFig'});
                if isempty(selectedFig),return; end
                ax = findobj(selectedFig.Children,'Type','Axes');
            end
            plotProfileLocations(obj,muicat,ax,caserec,false);
            ax.XLabel.String = 'Eastings';
            ax.YLabel.String = 'Northings';
        end
%%
        function [idp,caserec] = getProfileSelection(~,muicat)
            %get data to plot the locations of all the selected profiles
            %prompt user to select profiles to be used (multiple profiles)
            %returns id of profile in class list and caserec
            classname = 'ctBeachProfileData';
            caserec = strcmp(muicat.Catalogue.CaseClass,classname);
            promptxt = 'Select profiles to be used:';
            caselist = muicat.Catalogue.CaseDescription(caserec)';
            classid = muicat.Catalogue.CaseID(caserec);
            [idp,ok] = listdlg('Name','Shoreline location', ...
                            'PromptString',promptxt,'ListString',caselist);                        
            if ok<1, return; end
            %convert class id selection to caserec values
            nselect = length(idp);
            caserec = zeros(1,nselect);
            for i=1:nselect
                caserec(i) = find(muicat.Catalogue.CaseID==classid(idp(i)));
            end
        end
%%
        function ax = plotProfileLocations(obj,muicat,ax,caserec,isbaseplt)
            %genernate plot of profile locations based on caserec selection
            %check that E and N are included in data set             
            dst = getDataset(muicat,caserec(1),1);            
            varnames = dst.VariableNames;
            if ~any(strcmp(varnames,'Eastings'))
                warndlg('Timeseries does not include Easting and Northing co-ordinates');
                return;
            end
            
            nprof = length(caserec);            
            nint = subsampleprofiles(obj,nprof);
            pnam = cell(nprof,1); side = pnam; 
            Es = NaN(nprof,1); Ns = Es; pX = NaN(nprof,2); pY = pX;
            txt_angle = NaN(nprof,1); 
            for ix=1:nprof
                dst = getDataset(muicat,caserec(ix),1); 
                E = dst.Eastings;
                N = dst.Northings;           
                Ch = dst.Chainage;
                [~,idminCh] = min(Ch(:,1),[],'omitnan');       %minimum in first column
                [~,idch] = max(Ch,[],'all','omitnan','linear');%maximum in Ch
                [idrow,idmaxCh]  = ind2sub(size(Ch),idch);     %row and col of max
                Es(ix,1) = E(idminCh,1); Ee = E(idrow,idmaxCh);
                Ns(ix,1) = N(idminCh,1); Ne = N(idrow,idmaxCh);
                pX(ix,:) = [Es(ix,1),Ee];
                pY(ix,:) = [Ns(ix,1),Ne];                
                pnam(ix) = {dst.Description};
                %get angle - use atan2 to avoid divide by zero
                txt_angle(ix,1) = rad2deg(atan2((Ne-Ns(ix,1)),(Ee-Es(ix,1))));
                %adjust alignment to be left or right of coordinates and
                %ensure that text is upright by keeping angle to =/-90
                if txt_angle(ix)>-90 && txt_angle(ix)<90
                    side(ix) = {'right'};                    
                else
                    side(ix) = {'left'};
                    txt_angle(ix) = txt_angle(ix)-180*sign(txt_angle(ix));
                end       
            end 
            %sort profiles into alongshore order
            [sortedE,sortedN,idd] = sortENdata2line(Es,Ns);
            pX = pX(idd,:);
            pY = pY(idd,:);
            pnam = pnam(idd);
            txt_angle = txt_angle(idd);
            side = side(idd);
            
            hold on
            for ij = 1:nint:nprof
                h1 = plot(ax,pX(ij,:),pY(ij,:));
                if isbaseplt
                    h1.Color = [0.8,0.8,0.8];
                end
                % Exclude line from legend
                set(get(get(h1,'Annotation'),'LegendInformation'),...
                                                'IconDisplayStyle','off');
                txt = text(ax,pX(ij,1),pY(ij,1),pnam{ij},'FontSize',6,...
                 'HorizontalAlignment',side{ij}); 
                txt.Rotation = txt_angle(ij);                
            end
            h2 = plot(ax,sortedE,sortedN,':k');
            set(get(get(h2,'Annotation'),'LegendInformation'),...
                                                'IconDisplayStyle','off');
            hold off   
        end
%%
        function nint = subsampleprofiles(~,nprof)
            %prompt user for the number of intervals to sample at
            nint = 1;
            if nprof>20
                questxt = sprintf(...
                'There are %g profiles\nDo you want to sub-sample?',nprof);
                answer = questdlg(questxt,'BeachProfiles','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    sub = inputdlg('What intervals?','BeachProfiles',1,{'2'});
                    if isempty(sub), sub = '2'; end
                    nint = str2double(sub);
                end
            end
        end
    end
%%      
    methods (Access=private)
        function obj = addBPdataFile(obj,addedprof)
            %load data from a single data file (each file can have multiple 
            %profiles for a given date (CCO format)
            % localObj - instance of BeachProfileData class
            % classrec - id of the record for the profile to be used 
            % addedprof - the profile to be added to the record
            existprof = obj.Data.Dataset;
            %find the start and end of the new and existing data sets
            oldlen = length(existprof.DataTable{1,1});
            newlen = length(addedprof.DataTable{1,1});
            diff = abs(newlen-oldlen);
            pad = ones(1,diff)*NaN;
            %check if new and existing profiles are same length and if not 
            %pad the shorter profile(s) with NaNs   
            nrec = height(existprof.DataTable);
            varnames = existprof.VariableNames;
            if newlen>oldlen
                for it = 1:nrec
                    for jv = 1:length(varnames)
                        existprof.(varnames{jv})(it,oldlen+1:newlen) = pad;
                    end
                end
            elseif newlen<oldlen
                for jv = 1:length(varnames)
                    addedprof.(varnames{jv})(1,newlen+1:oldlen) = pad;
                end
            end
            
            newprof = mergerows(existprof,addedprof); 
            obj.Data.Dataset= newprof;
        end
    end
end