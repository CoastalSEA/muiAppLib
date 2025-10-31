classdef (Abstract = true) waveModels < muiDataSet
%
%-------abstract class help------------------------------------------------
% NAME
%   waveModels.m
% PURPOSE
%   Abstract class wave models to define data access methods
% NOTES
%  Uses muiDataSet as superclass and is used as an abstract class for wave
%  models in CoastalTools and the WaveRayModel Apps
% SEE ALSO
%   ctWaveModel.m and WRM_WaveModel.m
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2025
%--------------------------------------------------------------------------
% 
    properties (Abstract)  %properties that all subclasses must include
        ModelType          %model used for the particular instance
    end

%%
    methods
        function [tsdst,caserec] = getWaveModelDataset(obj,mobj,type,addnames,caserec)
            %prompt user to select a model wave dataset and add Tp if inshore
            %used to access wave model data created using ctWaveModel or
            %WRM_WaveModel from methods in CT_WaveModels, CT_BeachAnalysis
            %ct_coastal_plots, etc 
            muicat = mobj.Cases;
            if nargin<4
                addnames = {'Tp'};  %default is to add Tp
            end
            %
            if nargin<5   %no caserec to prompt for selection
                [wvobj,wvdst,ok] = selectClassInstance(obj,'ModelType',type);
                if ok<1, tsdst = []; return; end
                caserec = caseRec(muicat,wvobj.CaseIndex);
            else
                wvobj = getCase(muicat,caserec);
                wvdst = wvobj.Data.Dataset;
            end

            %if inshore wave dataset add variables requested, otherwise just
            [tsdst,timerange] = getSubSet(obj,wvdst);          %allow user to extract a subset 
            if strcmp(type,'Inwave_model')
                inpwavecid = wvobj.RunParam.ctWaveData.caseid; %source dataset (offshore)
                inpwaverec = caseRec(muicat,inpwavecid);       %case record
                srcdst = getDataset(muicat,inpwaverec,1);      %dst used to create inshore waves
                %check that source is the same length as selected wave dataset
                if height(srcdst)~=height(wvdst)    
                    msg = sprintf('Selected wave case and source dataset are different lengths\nUnable to add additional variables in getWaveModelDataset');
                    warndlg(msg); tsdst = []; return;
                end
                inpdst = removerows(srcdst,find(~timerange));  %match to selection
                dstnames = inpdst.VariableNames;               %source variables
                varnames = tsdst.VariableNames;                %inshore wave variables
                for i=1:length(addnames)                    
                    if any(strcmp(varnames,addnames{i}))
                        %variable already included
                        continue;
                    elseif any(contains(varnames,addnames{i})) && ...
                                      ~any(strcmp(varnames,addnames{i}))  
                        %source dataset contains 'addname' but not an exact 
                        %match. eg 'Tpi' instead of 'Tp'.if multiple matches 
                        %all are added - use with care!!
                        idx = contains(varnames,addnames{i});
                        tsdst.VariableNames{idx} = addnames{i};
                    elseif any(strcmp(dstnames,addnames{i})) && ...
                                            ~isempty(inpdst.(addnames{i}))
                        %variable to be added exists in source dataset
                        tsdst = addvars(tsdst,inpdst.(addnames{i}),...
                                           'NewVariableNames',addnames{i});
                    else                        
                        warndlg('Variable %s not found so not added to wave dataset',...
                                                           addnames{i});
                        tsdst = [];
                    end
                end
            end
        end
    end

%%
     methods (Access = protected)
        function [tsdst,meta] = getInputData(obj,mobj)
            %prompt user to select wave and water level data and return in
            %input dstable of data and metadata for inputs used
            tsdst = []; meta.iselvar = false;
            muicat = mobj.Cases;
            wvclassops = {'ctWaveData','muiUserModel'};
            promptxt = 'Select input wave data set:';           
            [wv_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                           'CaseClass',wvclassops,'ListSize',[300,100]);                                    
            if ok<1, return; end
            wvdst = getDataset(muicat,wv_crec,1);    %1 selects first dataset in struct
                                                     %ie Dataset or Spectra in most cases
            inputxt = sprintf('%s used for offshore waves',wvdst.Description);
            [wvdst,timerange] = getSubSet(obj,wvdst);%allow user to extract a subset 
            wvtime = wvdst.RowNames;               
            
            if isfield(wvdst.Dimensions,'freq')
                meta.source = 'Measured spectra'; 
                tsprops = getDataset(muicat,wv_caserec,2);  %2 selects Properties dataset from spectra input
                tsprops = removerows(tsprops,find(~timerange));
                wvdst = horzcat(wvdst,tsprops);               
                wvdst = activatedynamicprops(wvdst);              
            else
                meta.source = 'Measured waves';
                %check whether default variable names are not used and selection needed
                varnames = wvdst.VariableNames;
                if ~any(strcmp(varnames,'Hs'))
                    wvdst = extract_wave_data(wvdst);
                    if isempty(wvdst), return; end
                    meta.iselvar = true;  %variables selected (non-standard names)
                end
            end
            
            wlclassprops = {'ctWaterLevelData','ctTidalAnalysis','muiUserModel'};                                              
            promptxt = 'Select input water level data set (Cancel to use SWL=0):';           
            [wl_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                             'CaseClass',wlclassprops,'ListSize',[300,100]);                                 
            swl = zeros(size(wvtime));               
            if ok<1 || isempty(wl_crec)
                getdialog('Using SWL=0');
                inputxt = sprintf('%s, 0mOD used for water level',inputxt);
                wl_crec = 0;     %assign a null value if no water level data available
            else
                wldst = getDataset(muicat,wl_crec,1); 
                %check that there is water level data for period of interest
                [idst,idnd] = ts2_endpoints_in_ts1(wvdst,wldst);
                if isempty(idst)
                    getdialog('Data do not overlap. Using SWL=0');
                    inputxt = sprintf('%s, 0mOD used for water level',inputxt);
                    wl_crec = 0; %assign a null value if no water level data available
                else 
                    %select a variable from the water level dataset
                    varnames = wldst.VariableNames;
                    idx = 1;
                    if length(varnames)>1
                        [idx,ok] = listdlg('Name','WL options', ...
                            'PromptString','Select a variable:', ...
                            'SelectionMode','single','ListSize',[200,100],...
                            'ListString',varnames);
                        if ok<1, idx = 1; end
                    end
                    wldata = wldst.(varnames{idx});
                    wltime = wldst.RowNames;
                    swltime = wvtime(idst:idnd);                                        
                    %now interpolate water levels onto wave height times
                    swl(idst:idnd,1) = interp1(wltime,wldata,swltime,'linear','extrap');
                    swl(isnan(swl)) = 0;

                    inputxt = sprintf('%s, %s used for water levels',...
                                                inputxt,wldst.Description);
                end
            end
            % 
            % [wvdst,timerange] = getSubSet(obj,wvdst);      %allow user to extract a subset 
            % swl = swl(timerange);

            tsdst = addvars(wvdst,swl,'NewVariableNames','swl');
            meta.inptxt = inputxt;

            %assign the run parameters to the model instance
            if wl_crec==0
                meta.caserecs = {wv_crec};
            else
                meta.caserecs = {wv_crec,wl_crec};    
            end
        end

%%
        function [subdst,timerange] = getSubSet(~,tsdst)
            %subsample the record based on user defined start and end date
            subdst = getsampleusingrange(tsdst);
            times = tsdst.RowNames;
            timerange = ismember(times,subdst.RowNames);                      
        end
     end
end