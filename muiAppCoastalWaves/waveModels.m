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
            % mobj - model instance
            % type - 
            % addnames - 
            % caserec - id for specific case to use 

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
            [tsdst,timerange] = waveModels.getSubSet(wvdst);   %allow user to extract a subset 
            if strcmp(type,'Inwave_model')
                %inshore wave model data set so add variables requested
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
            else
                %offshore wave model data set - nothing to add
            end
        end
    end

%%
    methods (Static)
        function [cobj,tsdst,meta] = getCaseInputParams(mobj,varargin)
            %get the Case, Dataset and Input parameters
            %this combines calls to getCaseDataset and getInputParams
            %second input can be any value. If exists the limit classes to
            %select from to ctWaveData
            % varargin can include class selection optionss and dataset index
            % eg: {'ctWaveSpectrumData'},1
            meta = []; 
            if nargin>1
                %limit the classes that can be selected and specify dataset
                %as first one
                [cobj,tsdst,dsnames] = waveModels.getCaseDataset(mobj,varargin{:});
                % if isempty(dsnames) || ~any(contains(dsnames,dtype))
                %     cobj = []; tsdst = []; return; 
                % end
            else
                [cobj,tsdst,dsnames] = waveModels.getCaseDataset(mobj);
            end
            if isempty(tsdst), return; end
            [tsdst,meta] = waveModels.getInputParams(cobj,tsdst,dsnames);  %extract required variables
            if isempty(tsdst), return; end
            tsdst(1).DataTable = rmmissing(tsdst(1).DataTable);%remove nans
        end

%%
        function [cobj,tsdst,dsnames] = getCaseDataset(mobj,classops,idd)
            %get selection and load case. option to limit classopt in call
            %classops and idd optional - allows specific class and dataset 
            %to be selected in the call
            if nargin<2
                idd = [];
                classops = {'ctWaveData','ctWaveSpectrumData','ctWindData',...
                                          'WRM_WaveModel','muiUserModel'};
            elseif nargin<3
                idd = [];                
            end
            promptxt = 'Select input data to use:';
            [cobj,~,dsnames,idd] = selectCaseDataset(mobj.Cases,...
                                          [],classops,promptxt,idd);
            if isempty(cobj) || isempty(idd)
                tsdst = []; dsnames = []; return; 
            end
            tsdst = cobj.Data.(dsnames{idd});            
        end

%%        
function [xtsdst,meta] = getInputParams(cobj,tsdst,dsnames)
            %check for valid variable names when timeseries wave or wind
            %data are used to define conditions - used in ctWaveSpectraPlots
            xtsdst = []; meta = []; xmeta = [];
            varnames = tsdst.VariableNames;
            if any(strcmp(dsnames,'sptSpectrum')) && any(strcmp(varnames,'Kurt'))
                inptype = 'Spectrum';
                iselvar = false;  %variables selected (non-standard names)
                xtsdst = tsdst;
                tsprops = cobj.Data.(dsnames{2});  %2 selects Properties dataset from spectra input
                xtsdst(2) = tsprops;
            elseif any(strcmp(dsnames,'sptProperties'))  
                inptype = 'Spectrum';
                iselvar = false;  %variables selected (non-standard names)
                tsspec = cobj.Data.(dsnames{1});  %1 selects Spectra dataset from spectra input
                xtsdst = tsspec;
                xtsdst(2) = tsdst;
            elseif any(contains(dsnames,'Spectra')) && any(strcmp(varnames,'Kurt'))
                inptype = 'Spectrum';
                iselvar = false;  %variables selected (non-standard names)
                xtsdst = tsdst;
            else
                iselvar = true;  %variables selected (non-standard names)
                if isa(cobj,'ctWaveData')
                    inptype = 'Wave';
                    [xtsdst,xmeta] = extract_wave_data(tsdst); %returns 1xN array if multi-modal                    
                elseif isa(cobj,'WRM_WaveModel') 
                    inptype = 'Wave';
                    [xtsdst,xmeta] = extract_wave_data(tsdst); %returns 1xN array if multi-modal
                elseif isa(cobj,'ctWindData')
                    inptype = 'Wind';
                    [xtsdst,xmeta] = extract_wind_data(tsdst,1); %isfetch=true
                elseif isa(cobj,'muiUserModel') && contains(varnames,'Hs')  %NOT TESTED
                    inptype = 'Wave';
                    [xtsdst,xmeta] = extract_wave_data(tsdst); %returns 1xN array if multi-modal
                else
                    warndlg('Selection not yet handled in waveModels.getInputParams')
                    return
                end 
                if isempty(xmeta) || isempty(xmeta.selection), iselvar = false; end
            end
            if isempty(xtsdst), return; end

            meta = struct('inptype',inptype,'variables',xmeta,'iselvar',iselvar);                                                                
        end

%%
        function [tsdst,meta] = getInputData(mobj)
            %prompt user to select wave and water level data and return in
            %input dstable of data and metadata for inputs used in
            %CT_WaveModels - this function adds water levels if available
            % (could possibly be merged with getInputParams)
            meta.iselvar = false;
            classops = {'ctWaveData','muiUserModel','ctWindData','ctWaveSpectrumData'};  
            [cobj,tsdst,dsnames] = waveModels.getCaseDataset(mobj,classops);
            if isempty(tsdst), return; end

            [tsdst,meta] = waveModels.getInputParams(cobj,tsdst,dsnames);
            [tsdst(1),timerange] = waveModels.getSubSet(tsdst(1));%allow user to extract a subset  
            for i=2:numel(tsdst)
                tsdst(i) = removerows(tsdst(i),find(~timerange));
            end

            meta.inptxt = sprintf('%s used for offshore waves',tsdst(1).Description); 
            meta.caserecs = caseRec(mobj.Cases,cobj.CaseIndex);
            %add water levels to wave dataset (only adds to first dstable
            %if tsdst is an array - eg for spectrum case)
            %adds fields for inputxt and caserecs to meta struct
            [tsdst,meta] = addwaterlevels2waves(tsdst,mobj,meta); 
            meta.caserecs = num2cell(meta.caserecs);
        end

%%
        function [subdst,timerange] = getSubSet(tsdst)
            %subsample the record based on user defined start and end date
            subdst = getsampleusingrange(tsdst);
            times = tsdst.RowNames;
            timerange = ismember(times,subdst.RowNames);                      
        end
    end
end