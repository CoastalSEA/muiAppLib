classdef ctWaveSpectraPlots < ctWaveSpectrum
%
%-------class help------------------------------------------------------
% NAME
%   ctWaveSpectraPlots.m
% PURPOSE
%   Class for analysing wave spectra data held as spectral density as a
%   function of direction and frequency, or loaded from a file.
% NOTES
%   Spectral data are imported to the ctWaveData class using functions
%   such as wave_cco_spectrum. Spectral records hold two datasets named as
%   Spectra and Properties. Wave data can be simple unimodal or multimodal
%   descriptions of the sea state.
% SEE ALSO
%   see ct_costal_plots and SpectralTransfer for examples of use
%   inherits ctWaveSpectrum class
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2025
%--------------------------------------------------------------------------
%    
    properties 
        % inherits properties from ctWaveSpectrum
        hFig = gobjects(0)         % figure handle
        pType                      % figure type selection - used in getMultiPlot
    end

    methods
        function obj = ctWaveSpectraPlots
            %class constructor
            % default properties set in ctWaveSpectrum constructor
        end
    end

%% ------------------------------------------------------------------------
% Plotting and analysis call functions
%--------------------------------------------------------------------------
    methods (Static)
        function obj = getPlotOption(mobj)
            %user selected plotting options (used in CoastalTools)
            listxt = {'Plot a spectrum using Case data',...
                      'Plot a spectrum using a model',...
                      'Plot a spectrum loaded from file',...
                      'Plot comparison of Case spectra',...
                      'Plot Measured against Modelled',...
                      'Animation of spectrum timeseries',...
                      'Bimodal analysis of Measured spectrum',...
                      'Fit model to Measured spectra timeseries',...
                      'JONSWAP gamma from Measured timeseries',...
                      'Compare Modelled and Measured timeseries',...                      
                      'Subsample SPT spectrum timeseries',...
                      };

            selection = listdlg("ListString",listxt,"PromptString",...
                            'Select option:','SelectionMode','single',...
                            'ListSize',[260,180],'Name','Wave spectra');
            if isempty(selection), return; end
            switch selection
                case 1                  %Plot a spectrum using case data
                    obj = ctWaveSpectraPlots.plotCaseSpectrum(mobj);
                case 2                  %Plot a spectrum using a model
                    obj = ctWaveSpectraPlots.plotModelSpectrum();
                case 3                 %Plot a spectrum loaded from file
                    obj = ctWaveSpectraPlots.plotFileSpectrum();
                case 4                  %Plot a comparison of spectra for any 2 cases 
                    obj = ctWaveSpectraPlots.cfCaseSpectra(mobj);
                case 5                  %Fit model to measured spectrum using spectrum params
                    obj = ctWaveSpectraPlots.cfSpectrum2Model(mobj);
                case 6                  %Timeseries animation
                    ctWaveSpectraPlots.animateCaseSpectrum(mobj);
                case 7                  %decompose measured spectrum to bimodal form
                    ctWaveSpectraPlots.bimodalSpectrum(mobj);
                case 8                  %Fit a model to a timeseries of measured spectra
                    ctWaveSpectraPlots.cfModelled2Measured_ts(mobj,false);
                case 9                 %estimate JONSWAP gamma from spectra timeseries'
                    ctWaveSpectraPlots.estimateSpectrumGamma(mobj);
                case 10                  %compare fit of modelled wave data to measured spectra
                    ctWaveSpectraPlots.cfModelled2Measured_ts(mobj,true);                
                case 11                 %subsample SPT format spectrum timeseries
                    ctWaveSpectraPlots.subsampleSpectrum(mobj);
            end
        end

%%
        function runPlotOption(src,~,mobj)
            %main menu callback as alternative to getPlotOption (used in
            %WaveRayModel)
            switch src.Text
                case 'Case'                   %Plot a spectrum using case data
                    ctWaveSpectraPlots.plotCaseSpectrum(mobj);
                case 'Model'                  %Plot a spectrum using a model
                    ctWaveSpectraPlots.plotModelSpectrum();
                case 'SPT File'               %Plot a spectrum loaded from file
                    ctWaveSpectraPlots.plotFileSpectrum();
                case 'cf Cases'               %Plot a comparison of spectra for any 2 cases 
                    ctWaveSpectraPlots.cfCaseSpectra(mobj);
                case 'SPT v Model'            %Fit model to measured spectrum using spectrum params
                    ctWaveSpectraPlots.cfSpectrum2Model(mobj);
                case 'Animation'              %Timeseries animation
                    ctWaveSpectraPlots.animateCaseSpectrum(mobj);
                case 'Bimodal analysis'       %decompose measured spectrum to bimodal form
                    ctWaveSpectraPlots.bimodalSpectrum(mobj);    
                case 'FitModel2Measured'      %Fit model to a timeseries of measured spectra
                    ctWaveSpectraPlots.fitModel2Measured_ts(mobj);
                case 'Estimate JONSWAP gamma'  %estimate JONSWAP gamma from spectra timeseries'                    
                    ctWaveSpectraPlots.estimateSpectrumGamma(mobj);
                case 'cfModel&&Measured'       %Model case timeseries and cf to measured spectra timesereis
                    ctWaveSpectraPlots.cfModelled2Measured_ts(mobj);                
                case 'Subsample timeseries'   %subsample SPT format spectrum timeseries
                    ctWaveSpectraPlots.subsampleSpectrum(mobj);
            end
        end

%% ------------------------------------------------------------------------
% Static plotting methods
%--------------------------------------------------------------------------
        function [obj,ptype] = plotCaseSpectrum(mobj,vis,~)
            %Plot a spectrum using case data
            if nargin<2, vis = 'on'; end
            ptype = ctWaveSpectraPlots.plotType();

            %get Case dataset to be used
            [~,tsdst,meta] = waveModels.getCaseInputParams(mobj);
            if isempty(tsdst), obj = []; return; end
        
            %select record from dateset get spectrum and plot
            dates = tsdst(1).DataTable.Properties.RowNames;
            ok = 0;
            spselect = [];
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                obj = ctWaveSpectraPlots;             %initialise class object
                obj.spModel.selection = spselect;     %add any existing selection
                obj = getSpectrumObject(obj,meta.inptype,tsdst,irow);
                if isempty(obj.Spectrum.SG), return; end
                if ~isempty(obj.spModel)
                    spselect = obj.spModel.selection; %update model selection
                end

                %call plot function   
                getPlot(obj,ptype,vis);
                if nargin==3, ok = 1; end
            end
        end

%%
        function obj = plotModelSpectrum()
            %Plot a spectrum using a model and user defined wave/wind parameters
            ptype = ctWaveSpectraPlots.plotType();
     
            ok = 0;
            while ok<1    
                obj = ctWaveSpectraPlots;
                obj = setForcingConditions(obj); %UI to input wave/wind conditions
                if isempty(obj.inpData), ok = 1; continue; end

                obj = getSpectrum(obj);          %define the model to be used (Jonswap etc)
                if isempty(obj.Spectrum.SG), ok = 1; continue; end

                inputMessage(obj);                %display inputs in command window   
                display(obj.Params);              %table of combined, wind-wave and swell

                %call plot function   
                getPlot(obj,ptype);
            end
        end

%%
        function obj = plotFileSpectrum()
            %Plot a spectrum loaded from file
            ptype = ctWaveSpectraPlots.plotType();
            obj = ctWaveSpectraPlots();
            
            ok = 0;
            while ok<1 
                obj = setSpectrumFile(obj);        %load spectrum data from file
                if isempty(obj),return; end

                obj = getMeasuredSpectrum(obj);    %compute spectrum based on measured form
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                
                %call plot function
                getPlot(obj,ptype);
            end
        end

%%
        function obj = cfCaseSpectra(mobj)
            %plot a comparison of any two Cases in a single plot
            obj(2) = ctWaveSpectraPlots;
            summary = [];
            for i=1:2
                obj(i) = ctWaveSpectraPlots.plotCaseSpectrum(mobj,'off',1);
                obj(i).Params.Properties.RowNames = padRowNames(obj,i);
                summary = vertcat(summary,obj(i).Params); %#ok<AGROW>
            end
        
            %call plot function
            hf = getMultiPlot(obj);
            %add button to access wave parameters
            addDataButton(obj,hf,summary);
        end

%%
        function obj = cfSpectrum2Model(mobj)
            %plot a comparison of measured and modelled spectra
            ptype = ctWaveSpectraPlots.plotType();

            %get the measured wave spectrum to be modelled
            [~,tsdst,meta] = waveModels.getCaseInputParams(mobj,'sptSpectrum'); 
            if isempty(tdst), getdialog('Select wave spectra data'); return; end
            if ~contains(meta.inptype,'Spectrum')
                warndlg('Measured spectrum required for this option');
                obj = []; return;
            end
            propstable = tsdst(2).DataTable;  %extract sptProperties
            tsdst = tsdst(1);                 %assign sptSpectrum as tsdst

            %select record from dateset get spectrum and plot
            dates = tsdst.DataTable.Properties.RowNames;
            ok = 0; isfirst = true;
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                %construct measured spectrum
                obj = ctWaveSpectraPlots;     %initialise class object
                obj = getSpectrumObject(obj,meta.inptype,tsdst,irow);

                %add spectrum input properties
                obj.inpData.properties = propstable(irow,:);                
                %construct model spectrum
                obj(2) = setSpectrumModel(ctWaveSpectraPlots); %define the model to be used (Jonswap etc)
                obj(2) = getWaveModel(obj(2),obj(1).Params);
                if isempty(obj(2).Spectrum.SG), return; end

                %plot results
                obj(2) = setModelInputText(obj(2));   
                summary = compareProperties(obj);
                hf = plotObsModel(obj,ptype);
                %add button to access wave parameters
                addDataButton(obj,hf,summary);
                %plot omni-directional spectrum
                sax = omniSpectrumPlot(obj);
                subtitle(sax,obj(2).Plotxt.stxt)

                %plot skill
                skill = getSkillParameters(obj(1),mobj);      %initialise skill model
                metatxt = {obj(1).Plotxt.ttxt, obj(2).Plotxt.stxt};
                refvar = obj(1).Spectrum.SG;
                testvar = obj(2).Spectrum.SG;               
                option = 'Add';
                if isfirst                  
                    option = 'New'; isfirst = false;
                end
                taylor_plot(refvar,testvar,metatxt,option,3,skill);
            end
        end

%%
        function animateCaseSpectrum(mobj)
            %create animation from spectral data set or wave timeseries
            ptype = ctWaveSpectraPlots.plotType();

            %get the wave spectrum Case dataset to be animated
            [~,tsdst,meta] = waveModels.getCaseInputParams(mobj);
            if isempty(tsdst), return; end

            tsdst(1) = getsampleusingrange(tsdst(1));  
            if isempty(tsdst(1)), return; end         %invalid selection

            anobj = ctWaveSpectraPlots();
            if strcmp(meta.inptype,'Spectrum')
                obj = getMeasuredTS(anobj,tsdst(1));    %compute spectrum based on measured form
                obj(1).Plotxt.ttxt = tsdst(1).Description;
            else  
                obj = getModelTS(anobj,tsdst,meta.inptype); %compute spectrum for specified conditions
                if isempty(obj), return; end
                obj(1) = setModelInputText(obj(1));
                obj(1).Plotxt.ttxt = sprintf('%s\n%s',tsdst.Description,...
                                         obj(1).Plotxt.stxt);
            end
            
            isfixed = true;
            answer = questdlg('Allow z-axis scale to vary?','Animation','Yes','No','Yes');
            if strcmp(answer,'Yes'), isfixed = false; end
            wrm_single_animation(obj,mobj,ptype,isfixed);
        end

%%
        % function fitModel2Measured_ts(mobj)
        %     %create a measured timeseries of spectra and use the measured
        %     %properties to create a model spectrum to assess model skill
        %     [~,tsdst,meta] = waveModels.getCaseInputParams(mobj,'sptSpectrum'); 
        %     if isempty(obsdst), getdialog('Select wave spectra data'); return; end
        %     if isempty(tsdst), return; end
        %     if ~contains(meta.inptype,'Spectrum')
        %         warndlg('Measured spectrum required for this option');
        %         return;
        %     else
        %         tsdst = tsdst(1);                %assign sptSpectrum as tsdst
        %     end
        % 
        %     ok = 0;                              %with current workflow need to 
        %     while ok<1                           %limit size of timeseries
        %         tsdst = getsampleusingrange(tsdst);  
        %         if isempty(tsdst), return; end       %invalid selection
        %         nrec = height(tsdst);
        %         if nrec>10000
        %             qtxt = sprintf('Timeseries has %d records. Do you want to use shorter record?',nrec);
        %             answer = questdlg({qtxt},'Spectrum','Yes','No','Yes');
        %             if strcmp(answer,'No'), ok = 1; end                         
        %         else
        %             ok = 1;
        %         end
        %     end
        % 
        %     %running long timeseries of measured spectra can be limited by
        %     %memory so run inidividually in loop
        %     obj = ctWaveSpectraPlots;
        %     obj = setSpectrumModel(obj);
        %     [dir,freq] = spectrumDimensions(obj);
        %     obj.Spectrum.freq = freq;
        %     obj.Spectrum.dir = dir;
        %     skill = getSkillParameters(obj,mobj);    %get the parameters for skill model
        %     mtime = tsdst.RowNames;
        %     meta.source{1} = tsdst.Description;
        %     meta.source{2} = obj.spModel.form;
        % 
        %     hpw = PoolWaitbar(nrec, 'Processing skill statistics');  %and increment(hpw);
        %     parfor i=1:height(tsdst)                          %parfor loop
        %         itsdst = getDSTable(tsdst,i,[]);              %selected record
        %         obsobj = copy(obj);
        %         obsobj = setInputParams(obsobj,itsdst,'Spectrum');
        %         obsobj = getMeasuredSpectrum(obsobj);         %compute spectrum based on measured form
        %         params = wave_spectrum_params(obsobj);        %integral properties of spectrum    
        %         indst = dstable(params,'RowNames',mtime(i));
        %         modobj = copy(obj);
        %         modobj = setInputParams(modobj,indst,'Wave');
        %         modobj = getModelSpectrum(modobj);                
        %         stats(i) = get_spectrum_skill_stats(obsobj,modobj,skill);
        %         obsprops(i,:) = params;
        %         modprops(i,:) = wave_spectrum_params(modobj);
        %         increment(hpw);
        %     end
        %     delete(hpw)
            %code to run measured and model as complete arrays
            % anobj = ctWaveSpectraPlots;
            % %get timeseries of measured spectra
            % obsobj = getMeasuredTS(anobj,tsdst); 
            % 
            % spectra = [obsobj(:).Spectrum]; 
            % dates = [spectra(:).date];
            % intable = vertcat(obsobj(:).Params); %concatenate tables
            % indst = dstable(intable,'RowNames',dates);
            % indst.Description = sprintf('Parameters from %s',tsdst.Description);
            % modobj = getModelTS(anobj,indst,'Wave');
            % if isempty(modobj), return; end   %user cancelled
            % %check on gamma values
            % if anobj.spModel.gamma==0
            %     hf = figure('Tag','PlotFig'); ax = axes(hf);
            %     spec = vertcat(obsobj(:).Spectrum);
            %     % spgm = vertcat(modobj(:).spModel);
            %     % plot(ax,[spec(:).date],[spgm(:).T_gamma],'x')
            %     inpd = vertcat(modobj(:).inpData);
            %     plot(ax,[spec(:).date],[inpd(:).gamma],'x')
            % end

        %     %plot model skill and allow user to examine individual parameters
        %     ctWaveSpectraPlots.plotSpectrumModelSkill(stats,skill,meta)
        %     ctWaveSpectraPlots.parameterPlots(obsprops,modprops);
        % end
%%
        function bimodalSpectrum(mobj)
            %analyse measured spectrum for bi-modality and explore
            %representing this in a model
            %get the measured wave spectrum to be modelled
            [~,obsdst,meta] = waveModels.getCaseInputParams(mobj,'sptSpectrum');
            if isempty(obsdst), getdialog('Select wave spectra data'); return; end
            obsdst = obsdst(1);                %assign sptSpectrum as tsdst
            % if ~contains(meta.inptype,'Spectrum')
            %     warndlg('Measured spectrum required for this option');
            %     return;
            % else
            %     tsdst = tsdst(1);
            % end
            %select record from dateset get spectrum and plot
            dates = obsdst.DataTable.Properties.RowNames;
            obj = ctWaveSpectraPlots;             %initialise class object
            ok = 0;
            while ok<1
                irow = listdlg("PromptString",'Select event to plot',...
                    'SelectionMode','single','ListSize',[160,300],...
                    'ListString',dates);
                if isempty(irow), ok = 1; continue; end
                obj.Plotxt.ttxt = sprintf('%s (%s)',obsdst.Description,dates{irow});
                tsdstrow = ctWaveSpectrum.getDatasetRow(obsdst,irow);
                %set input parameters for selected record
                obj = setInputParams(obj,tsdstrow,meta.inptype);
                obj = getMeasuredSpectrum(obj);  %compute spectrum based on measured form                
                if isempty(obj.Spectrum.SG); return; end
                obj.Params = wave_spectrum_params(obj);
                freq = obj.Spectrum.freq;
                dir = obj.Spectrum.dir;
                SG = obj.Spectrum.SG;
                [idpks,idmn] = spectrumPeaks(obj);

                %construct the input spectrum
                if isempty(idmn)
                    %not bimodal use full spectrum to estimate parameters
                    w_params = wave_spectrum_params(SG,freq,dir,0);%0=diagnostics not required
                    s_params = w_params;  s_params{:,:} = 0;       %dummy table to allow concatanation
                    msgtxt = sprintf('Single peak at: %.1fs;',1/freq(idpks));
                else
                    SGw = SG(:,idmn:end);  %wind-wave energy
                    fw  = freq(idmn:end);  %wind-wave frequencies
                    w_params = wave_spectrum_params(SGw,fw,dir,0); %0=diagnostics not required
                    SGs = SG(:,1:idmn);    %swell energy
                    fs  = freq(1:idmn);    %swell frequencies
                    s_params = wave_spectrum_params(SGs,fs,dir,0); %0=diagnostics not required
                    npks = numel(idpks);
                    msgtxt = sprintf('%d peaks found at:',npks);
                    idpks = fliplr(idpks);  %for periods reverse the order
                    for i=1:npks
                        msgtxt = sprintf('%s %.1fs;',msgtxt,1/freq(idpks(i)));
                    end
                end
                getdialog(msgtxt,[],4);
                
                params(1) = dstable(w_params,'RowNames',obj(1).inpData.date);
                params(1).Description = obj(1).inpData.source;
                params(2) = dstable(s_params,'RowNames',obj(1).inpData.date);

                anobj = ctWaveSpectraPlots;
                anobj = setInputParams(anobj,params,'Wave');
                anobj = getSpectrum(anobj); 
                if isempty(anobj.Spectrum.SG), return; end
                obj(2) = anobj;

                if isempty(idmn)
                    summary = [w_params;obj(2).Params];
                    summary.Properties.RowNames = {'Wind-waves','Model W-W'};
                else
                    summary = [obj(1).Params;obj(2).Params(1,:);w_params;...
                           obj(2).Params(2,:);s_params;obj(2).Params(3,:)];
                    summary.Properties.RowNames = {'Combined','Model combined',...
                        'Wind-waves','Model W-W','Swell waves','Model swell'};
                end

                display(summary)
                ax = omniSpectrumPlot(obj);
                hold(ax,'on')
                plot(ax,[1,1]*w_params.Tp,ylim,'--b')
                if ~isempty(idmn)
                    plot(ax,[1,1]*s_params.Tp,ylim,'--r')
                    plot(ax,[1,1]*1/freq(idmn),ylim,'-.g')  
                    legend({'Measured','Modelled','Sea Tp','Swell Tp',...
                        'Crossover'},'Location','northeast');
                else
                    legend({'Measured','Modelled'},'Location','northeast');
                end
                hold(ax,'off')           

                %surface plot comparison
                hf = plotObsModel(obj,'XY');
                %add button to access wave parameters
                addDataButton(obj,hf,summary); 

                obj = ctWaveSpectraPlots; %reset to blank instance
            end
        end


%%
        function cfModelled2Measured_ts(mobj,iscase)
            %create a measured timeseries of spectra and compare with:
            % iscase=true: modelled spectra from wave properties defined
            %   by suitable wave data timeseries (e.g. Copernicus data) to
            %   examine spectrum model fit parameters for given sea location
            % iscase=false: use the measured properties to create a model 
            %   spectrum to assess model skill
            [~,obsdst,meta] = waveModels.getCaseInputParams(mobj,'sptSpectrum');
            if isempty(obsdst), getdialog('Select wave spectra data'); return; end
            obsdst = obsdst(1);                %assign sptSpectrum as tsdst

            % if ~contains(meta.inptype,'Spectrum')
            %     warndlg('Measured spectrum required for this option');
            %     return;
            % else
            % 
            % end

            ok = 0;                              %with current workflow need to 
            while ok<1                           %limit size of timeseries
                obsdst = getsampleusingrange(obsdst);  
                if isempty(obsdst), return; end       %invalid selection
                nrec = height(obsdst);
                if nrec>10000
                    qtxt = sprintf('Timeseries has %d records. Do you want to use shorter record?',nrec);
                    answer = questdlg({qtxt},'Spectrum','Yes','No','Yes');
                    if strcmp(answer,'No'), ok = 1; end                         
                else
                    ok = 1;
                end
            end
            
            if iscase           %cf measured and model from case
                %get Case dataset to be used
                [~,seldst,meta] = waveModels.getCaseInputParams(mobj);
                if isempty(seldst), return; end
    
                %match record to measured dataset
                [newtime,ids,ido] = ctWaveSpectraPlots.subSampleTime(seldst,obsdst);
                if ~isempty(ido)
                    subtable = obsdst.DataTable(ido,:); 
                    obsdst.DataTable = subtable;
                    obsdst.RowNames = newtime;
                end
    
                moddst = copy(seldst);
                for j=1:numel(moddst)
                    subtable = seldst(j).DataTable(ids,:);  
                    moddst(j).DataTable = subtable;
                    moddst(j).RowNames = newtime; %update times in case there is an offset
                end
            else
                moddst = [];   %needed for parfor loop
            end

            %running long timeseries of spectra using getMeasuredTS and 
            %getModelTS to create ctWaveSpectrum objects and then computing 
            %the statistice on the object arrays can be limited by memory. 
            %To avoid this compute statistics for one spectrum at a time.
            obj = ctWaveSpectraPlots;
            obj = setSpectrumModel(obj);
            [dir,freq] = spectrumDimensions(obj);
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir;
            skill = getSkillParameters(obj,mobj);    %get the parameters for skill model
            mtime = obsdst.RowNames;
            meta.source{1} = obsdst.Description;
            meta.source{2} = obj.spModel.form;
            
            hpw = PoolWaitbar(nrec, 'Processing skill statistics');  %and increment(hpw);
            parfor i=1:height(obsdst)                          %parfor loop
                itsdst = getDSTable(obsdst,i,[]);              %selected record
                obsobj = copy(obj);
                obsobj = setInputParams(obsobj,itsdst,'Spectrum');
                obsobj = getMeasuredSpectrum(obsobj);         %compute spectrum based on measured form
                obsprops(i,:) = wave_spectrum_params(obsobj);        %integral properties of spectrum  

                modobj = copy(obj);
                if iscase
                    itsdst = ctWaveSpectrum.getDatasetRow(moddst,i); %selected record
                    modobj = setInputParams(modobj,itsdst,'Wave');
                    modobj = getMultiModalSpectrum(modobj);   %compute spectrum for specified conditions
                    params = modobj.Params(1,:);    
                    params.Properties.RowNames = {};
                else
                    indst = dstable(obsprops(i,:),'RowNames',mtime(i));
                    modobj = setInputParams(modobj,indst,'Wave');
                    modobj = getModelSpectrum(modobj); 
                    params = wave_spectrum_params(modobj);
                end  
                modprops(i,:) = params;
                stats(i) = get_spectrum_skill_stats(obsobj,modobj,skill);                
                increment(hpw);
            end

            %plot model skill and allow user to examine individual parameters
            ctWaveSpectraPlots.plotSpectrumModelSkill(stats,skill,meta)
            ctWaveSpectraPlots.parameterPlots(obsprops,modprops);

            % ok = 0;                              %with current workflow need to 
            % while ok<1                           %limit size of timeseries
            %     obsdst = getsampleusingrange(obsdst);  
            %     if isempty(obsdst), return; end       %invalid selection
            %     nrec = height(obsdst);
            %     if nrec<10000
            %         ok = 1; 
            %     else
            %         getdialog(sprintf('Timeseries has %d records\nReduce to <10,000',nrec))
            %     end
            % end
            % 
            % anobj = ctWaveSpectraPlots;
            % %get timeseries of measured spectra
            % obsobj = getMeasuredTS(anobj,obsdst);   
            % 
            % %get Case dataset to be used
            % [~,seldst,meta] = waveModels.getCaseInputParams(mobj);
            % if isempty(seldst), return; end
            % 
            % %match record to measured dataset
            % [newtime,ids,ido] = ctWaveSpectraPlots.subSampleTime(seldst,obsdst);
            % 
            % if ~isempty(ido)
            %     obsobj = obsobj(ido);
            % end
            % 
            % moddst = copy(seldst);
            % for j=1:numel(moddst)
            %     subtable = seldst(j).DataTable(ids,:);  
            %     moddst(j).DataTable = subtable;
            %     moddst(j).RowNames = newtime; %update times in case there is an offset
            % end
            
            % modobj = getModelTS(anobj,moddst,meta.inptype); %compute spectrum for specified conditions
            % if isempty(modobj), return; end
            % modobj(1) = setModelInputText(modobj(1));
            % modobj(1).Plotxt.ttxt = sprintf('%s\n%s',moddst.Description,...
            %                                      modobj(1).Plotxt.stxt);
            % %check on gamma values
            % if (contains(anobj.spModel.form,'JONSWAP') || ...
            %     contains(anobj.spModel.form,'TMA')) && anobj.spModel.gamma==0                
            %     spec = vertcat(obsobj(:).Spectrum);
            %     % spgm = vertcat(modobj(:).spModel);
            %     % plot(ax,[spec(:).date],[spgm(:).T_gamma],'x')
            %     inpd = vertcat(modobj(:).inpData);
            %     nspec = numel(inpd(1).gamma);
            %     gam = reshape([inpd(:).gamma],nspec,[]);
            %     if ~all(gam==3.3,'all')
            %         hf = figure('Tag','PlotFig'); ax = axes(hf);
            %         plot(ax,[spec(:).date],gam(1,:),'x')
            %         for j=2:nspec
            %             hold on
            %             plot(ax,[spec(:).date],gam(j,:),'.')
            %             hold off
            %         end
            %         xlabel('gamma')
            %         ltxt = {'Wind','Swell 1','Swell 2'};
            %         legend(ltxt(1:nspec))
            %     end
            % end
            % 
            % %plot model skill and allow user to examine individual parameters
            % plotSpectrumModelSkill(obsobj,modobj,mobj,meta)
            % parameterPlots(obsobj,modobj);
        end

%%
        function estimateSpectrumGamma(mobj)
            %estimate JONSWAP gamma from spectra timeseries
            [cobj,obsdst,meta] = waveModels.getCaseInputParams(mobj,'sptSpectrum');
            if isempty(obsdst), getdialog('Select wave spectra data'); return; end
            spdst = obsdst(1);                %assign sptSpectrum as tsdst

            % if ~contains(meta.inptype,'Spectrum')
            %     warndlg('Measured spectrum required for this option');
            %     return;
            % else
            % 
            % end
            spdst = getsampleusingrange(spdst);  

            nrec = height(spdst);
            gamma = zeros(nrec,1);
            S = spdst.S;
            f = spdst.Dimensions.freq;            
            parfor i=1:nrec
                gamma(i) = wave_spectrum_gamma(S(i,:),f); %#ok<PFOUS>
            end

            %get wave height threshold and create plot
            defaults = {'0','1'};
            answer = inputdlg({'Wave height threshold:','Save results (0/1)'},'Gamma',1,defaults);
            if isempty(answer), answer = defaults; end
            Hthr = str2double(answer{1});
            issave = logical(str2double(answer{2}));

            Hs = obsdst(2).Hs;
            idx = Hs<=Hthr;
            Hs(idx) = NaN;
            gamma(idx) = NaN;

            mn_gamma = mean(gamma,'all','omitnan');
            mn_Hs = mean(Hs,'all','omitnan');

            hf = figure('Name','Gamma','Tag','PlotFig');
            ax = axes(hf); 
            yyaxis left
            stem(ax,spdst.RowNames,gamma,'.','LineWidth',0.1)
            ylabel('Gamma (-)')
            xlabel('Time')
            yyaxis right            
            plot(ax,spdst.RowNames,Hs,'.')
            ylabel('Wave height (m)')
            glegtxt = sprintf('gamma with mean %.2f',mn_gamma);
            hlegtxt = sprintf('Hs with mean %.2f',mn_Hs);
            legend({glegtxt,hlegtxt})

            if issave
                meta.source = sprintf('Gamma using %s',spdst.Description);
                meta.data = sprintf('Wave height threshold %.2f',Hthr);
                dsp = ctWaveSpectraPlots.setDSroperties('gamma');
                ctWaveSpectraPlots.addORupdate(mobj,cobj,dsp,{gamma},meta);
            end
        end

%%
        function subsampleSpectrum(mobj)
            %get Case dataset to be used
            [cobj,dst,~] = waveModels.getCaseInputParams(mobj);
            if isempty(dst), return; end
            %method can interp1 method or none. 'none' finds exact match 
            %with a tolerance if tolerance>0 seconds. 
            inp = inputdlg({'Method (none, linear, etc)','Tolerance (s)'},...
                                            'Subsample',1,{'none','900'});
            if isempty(inp), return; end

            newdst = subsample_spectra_ts(dst,mobj,inp{1},str2double(inp{2}));
            if isempty(newdst), return; end

            %save sumsampled dataset
            classname = metaclass(cobj).Name;
            heq = str2func(classname);
            obj = heq();  %new instance of class object
            obj.Data = newdst;    %newdst is a struct
            setCase(mobj.Cases,obj,'data');
            getdialog(sprintf('Resampled dataset saved as %s',classname));
        end
    end
%% ------------------------------------------------------------------------
% Plotting functions
%--------------------------------------------------------------------------
    methods
        function ax = getPlot(obj,ptype,vis,ax)
            %call plot function   
            if nargin<3, vis = 'on'; ax = []; end
            if nargin<4, ax = []; end
            obj.pType = ptype;

            if strcmp(ptype,'XY')                %single x-y plot
                ax = surfPlot(obj,vis,ax);
            else
                ax = polarPlot(obj,vis,ax);
            end  
        end

%%
        function ax = surfPlot(obj,vis,ax) 
            %plot spectral surface as XY function of frequency and direction
            hf = [];
            if nargin<3 || isempty(ax)
                idx = numel(obj.hFig)+1;
                hf = figure('Name','SpecTrans','Tag','PlotFig','Visible',vis);
                obj.hFig(idx) = hf; %used in getMultiPlot
                ax = axes(hf);
            end
            %extract data
            ptxt = obj.Plotxt;
            params = obj.Params;
            period = 1./obj.Spectrum.freq;
            dir = obj.Spectrum.dir; 
            
            SG = obj.Spectrum.SG;
            SGpk = max(SG,[],'all');
            %make plot
            surf(ax,period,dir,SG,'Tag','PlotFigSurface');
            view(2);
            shading interp
            axis tight
            if ~isempty(params)
                hold on
                nparams = height(params);
                if nparams>2
                    for i=2:nparams
                        plot3(ax,params.Tp(i),params.Dir(i),...
                              SGpk*[1,1],'+y','MarkerSize',12,'Tag','DirPk')
                    end
                end
                plot3(ax,xlim,params.Dir(1)*[1 1],SGpk*[1,1],'--w','Tag','DirPk')
                plot3(ax,params.Tp(1)*[1 1],ylim,SGpk*[1,1],'--w','Tag','TpPk')
                hold off  
            end
            xlabel(ptxt.xtxt)
            ylabel(ptxt.ytxt)

            %add the colorbar and labels   
            cb = colorbar;
            cb.Label.String = ptxt.vtxt;
            cb(1).Tag = 'WaveSpectrum';            

            title(ptxt.ttxt)
            if ~isempty(ptxt.stxt)
                subtitle(ptxt.stxt)
            end

            if ~isempty(hf)  %add button to access wave parameters
                addDataButton(obj,hf) 
            end
        end

%%
        function ax = polarPlot(obj,vis,ax)
            %plot spectral surface as polar function of frequency and direction          
            hf = [];
            if nargin<3 || isempty(ax)
                idx = numel(obj.hFig)+1;
                hf = figure('Name','SpecTrans','Tag','PlotFig','Visible',vis);
                obj.hFig(idx) = hf; %used in getMultiPlot
                ax = axes(hf);
            end

            %extract data
            ptxt = obj.Plotxt;
            period = 1./obj.Spectrum.freq;
            dir = obj.Spectrum.dir;
            SG = obj.Spectrum.SG;

            %make plot
            wid = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
            radrange = [1,20];
            %interpolate var(phi,T) onto plot domain defined by tints,rints
            tints = linspace(obj.Interp.tlim{:});  
            rints = linspace(obj.Interp.rlim{:}); 
            [Tq,Rq] = meshgrid(tints,rints); 
            warning('off',wid)
            vq = griddata(deg2rad(dir),period,SG.',Tq,Rq);
            vq(isnan(vq)) = 0;  %fill blank sector so that it plots the period labels
            warning('on',wid)
            color = mcolor('dark blue');

            [X,Y] = polarplot3d(vq,'plottype','surfn','TickSpacing',45,...
                'RadLabels',4,'RadLabelLocation',{20 'top'},...
                'GridColor',color,'TickColor',color,...
                'RadLabelColor',mcolor('dark grey'),...
                'RadialRange',radrange,'polardirection','cw','Axes',ax);
            view(2)
            shading interp 
            axis equal
            axis(ax,'off')   %remove axes from plot display
            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = ptxt.vtxt;

            text(max(X,[],'all')/2,max(Y,[],'all')*0.95,0,ptxt.xtxt); %radial label
            cb.Tag = 'WaveSpectrum';
            %improve positioning of colorbar
            cb.Position(3:4) = cb.Position(3:4)*0.6;

            title(ptxt.ttxt)
            if ~isempty(ptxt.stxt)
                subtitle(ptxt.stxt)
            end

            if ~isempty(hf)  %add button to access wave parameters
               addDataButton(obj,hf)  
            end
        end
    end

%% ------------------------------------------------------------------------
% Plotting utilities
%--------------------------------------------------------------------------  
    methods (Access={?WRM_WaveModel,?ctWaveSpectrum,?SpectralTransfer} )
        function hf = getMultiPlot(obj)
            %plot a set of axes in a single figure
            hf = figure('Name','SpecTrans','Tag','PlotFig');
            nplot = numel(obj);
            t = tiledlayout(hf, nplot, 1); % nplot rows, 1 column
            hfigs = [obj(:).hFig];
            pos2 = [0.55,0.05];
            for i=1:nplot
                % Move existing axes into tiles
                sax = findobj(hfigs(i).Children,'Type','axes');
                hcb = findobj(hfigs(i).Children,'Type','colorbar');
                % hpos = hcb.Position(1);
                tt = nexttile(t, i);
                delete(tt)           % removes the placeholder axes
                sax.Parent = t;   % Reparent to tiledlayout in figure
                sax.Layout.Tile = i; 
                if strcmp(obj(i).pType,'Polar')
                    %colorbar needs to be re-positioned for polar plots
                    hcb.Position = [0.83,pos2(i),0.023,0.4];
                end
            end
            delete(hfigs)
        end

%%
        function hf = plotObsModel(obj,ptype)
            %plot two spectra  - similar to off_in_plot in SpecralTransfer
            hf = figure('Name','SpecTrans','Tag','PlotFig');            
            ax = axes(hf);

            if strcmp(ptype,'XY')
                m = 3; n = 1;  %vertical plots of xy
                hf.Position(2) = hf.Position(2)/2;
                hf.Position(4) = hf.Position(4)*2;
            else
                m = 1; n = 2;  %horizontal plots if polar
                hf.Position(3) = hf.Position(3)*2;
            end

            Smax = max(max(obj(1).Spectrum.SG,[],'all'),max(obj(2).Spectrum.SG,[],'all'));

            s(1) = subplot(m,n,1,ax);
            getPlot(obj(1),ptype,'on',s(1));
            s(1).Title.String = 'Measured'; 
            s(1).CLim(2) = Smax;

            s(2) = subplot(m,n,2);
            getPlot(obj(2),ptype,'on',s(2));
            s(2).Title.String = 'Modelled';
            s(2).CLim(2) = Smax;    %force z-scale of cf plots to be same

            if strcmp(ptype,'XY')
                s(2).XLim = s(1).XLim;
                s(3) = subplot(m,n,3);
                idf = ismembertol(obj(2).Spectrum.freq,obj(1).Spectrum.freq,1e-3);
                spdiff = obj(1).Spectrum.SG-obj(2).Spectrum.SG(:,idf);
                spmax = max(spdiff,[],'all');
                idx = spdiff<spmax/100 & spdiff>-spmax/100;
                spdiff(idx) = 0; %remove small differences
                diffobj = copy(obj(1));
                diffobj.Spectrum.SG = spdiff;
                diffobj.Plotxt.vtxt = 'Spectral Energy (m^2s)';
                getPlot(diffobj,ptype,'on',s(3));  
                s(3).Title.String = 'Difference [Meas-Model]';
                s(3).Subtitle.String = [];
                s(3).XLim = s(1).XLim;
            end

            sgtitle(sprintf('Spectrum for: %s',obj(1).Plotxt.ttxt))
        end 

%%
        function ax = omniSpectrumPlot(obj,ax)
            %plot omni-direction wave spectrum
            if nargin<2 || isempty(ax)
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end
            
            ptxt = obj.Plotxt;
            per = 1./obj(1).Spectrum.freq;
            hold(ax,'on')
            for i=1:length(obj)
                Sf = getOmniDirSpectrum(obj(i));
                plot(ax,per,Sf);
            end
            hold(ax,'off')
            xlabel(ptxt(1).xtxt)
            ylabel(ptxt(1).ytxt) 
            title(ax,ptxt(1).ttxt) 
        end

%%
        function addShoreNormal(obj,shorenorm)
            %add the shorenormal to a plot      
            grey = mcolor('light grey');
            ax = findobj(obj.hFig,'Type','axes');
            zMax = [1,1]*(ax.ZLim(2)-1);
            hold on
            if strcmp(obj.pType,'Polar')
                maxT = ax.XLim(2);
                ang = compass2trig(shorenorm);
                xn = [0,maxT*cos(ang)]; 
                yn = [0,maxT*sin(ang)];
                plot3(ax,xn,yn,zMax,'Color',grey,'LineStyle',':','LineWidth',1);  
            else
                xn = ax.XLim; 
                yn = [1,1]*shorenorm; 
                plot3(ax,xn,yn,zMax,'Color',grey,'LineStyle',':');           
            end
            hold off 
        end

%%
        function out = compareProperties(obj)
            %write table of values to the command window
            inp = table2struct(obj(1).inpData.properties);
            m0 = (inp.Hs/4)^2;
            [~,idx] = max(obj(1).inpData.spectrum.S);        
            inp.Dp = obj(1).inpData.spectrum.Dir(idx);
            inp.Tp = 1/obj(1).inpData.spectrum.Dimensions.freq(idx);
            source = table(inp.Hs,m0,NaN,inp.Sp,inp.Tp,inp.Dp,NaN,NaN,NaN,inp.Tz); 
            source.Properties.VariableNames = {'Hs','m0','Dir','Sp','Tp','Dp',...
                                                    'Sfdpk','Dfdpk','Tfdpk','T2'};
            out = [source;obj(1).Params;obj(2).Params];
            out.Properties.RowNames = {'Source properties','Source spectrum',...
                                                 'Model spectrum'};
            display(out)
        end

%%
        function addDataButton(obj,hf,outtable)
            %add button to allow user to summary statistics
            if nargin<3
                nvar = 3;                         %number of variables in table
                outtable = [];
                inpdata = []; rownames = {};
                npam = 0; 
                if ~isempty(obj.Params)
                    outtable = obj.Params;
                    nvar = width(outtable);
                    npam = height(outtable);
                end
    
                if ~isempty(obj.inpData) && isfield(obj.inpData,'Hs')
                    nrow = numel(obj.inpData.Hs); %number of rows to add                
                    inpdata = zeros(nrow,nvar);
                    rownames = cell(nrow,1);
                    for i=1:nrow
                        inpdata(i,1:5) = [obj.inpData.Hs(i),0,obj.inpData.Dir(i),...
                                                         0,obj.inpData.Tp(i)];
                        rownames{i} = sprintf('Input parameters %d',i);
                    end
                end
                outtable = [outtable;num2cell(inpdata)];
                if isempty(outtable), return; end
                outtable.Properties.RowNames(npam+1:end) = rownames;
            end

            pvtcb =  @(src,evdat)ctWaveSpectraPlots.dataFigure(src,evdat);
            uicontrol('Parent',hf,'Style','pushbutton',...
                'String','Results','Tag','FigButton',...
                'TooltipString','Access wave parameters',...
                'Units','normalized','Position',[0.88 0.96 0.10 0.04],...
                'UserData',outtable,...
                'Callback',pvtcb);           
        end 

%%
        function rnames = padRowNames(obj,idx)
            %ensure that row names are all unique by adding selection idx
            rnames = obj(idx).Params.Properties.RowNames;
            for j=1:numel(rnames)
                rnames{j} = sprintf('Case %d: %s',idx,rnames{j});
            end
        end

%% ------------------------------------------------------------------------
% Skill functions - not wave spectra specific - copied from ModelSkill
%--------------------------------------------------------------------------    
        function skill = getSkillParameters(obj,mobj)
            %extract Skill parameters using muiSkill_RunParams class for input
            x = obj(1).Spectrum.dir;              %direction
            y = obj(1).Spectrum.freq;             %frequency
            robj = muiSkill_RunParams.setInput(mobj);   %default or current values if user cancels

            skill.Ro = robj.maxcorr;
            skill.n  = robj.skillexponent;
            skill.Inc = true;                   %flag to include skill score
            skill.W = robj.skillwindow;
            if isempty(skill.W) || skill.W==0, skill.Inc = false; end
            subdomain = robj.skillsubdomain;
            skill.SD = ctWaveSpectraPlots.getSubDomain(x,y,subdomain);
            skill.iter = robj.skilliteration;
        end
    end

%% ------------------------------------------------------------------------
% Static plotting utilities
%--------------------------------------------------------------------------  
    methods (Static, Access=private)
%%
        function plotSpectrumModelSkill(stats,skill,meta)
            %compute the skill of model v measured spectrum data and produce Taylor
            %plot of timeseries results
            %calls muiSkill_RunParams class and taylor_plot function    
            % skill = getSkillParameters(obsobj(1),mobj);    %get the parameters for skill model
            % 
            % %get the statistics for the Taylor plot
            % stats = get_spectrum_skill_stats(obsobj,modobj,skill);
        
            %add the timeseries results to the Taylor plot
            hw = waitbar(0,'Preparing plot please wait');
            ndteststd = [stats(:).teststd]./[stats(:).refstd]; %normalised std
            rLim = ceil(max(ndteststd));                       %radial limit for the plot
            if rLim>6, rLim = 6; end
            ax = taylor_plot_figure(rLim);    
            metatxt = {'Measured','Model'};
            %local skill is not plotted even if computed but is reported in
            %the table of results on the Case list button
            ax = taylor_plot_ts(ax,stats,skill,metatxt); 
            waitbar(1,hw)
            ttxt = meta.source{2};
            if meta.iselvar
                nvar = size(meta.variables.selection,1);
                ttxt = sprintf('%s using %d sea states',ttxt,nvar);
            end
            ax.Title.String = ttxt;
            ax.Subtitle.String = sprintf('Using %s as reference',meta.source{1});
            delete(hw)
        end 

%%
        % function parameterPlots(obsobj,modobj)
        %     %allow user to select from the parameters table and plot
        %     %measured against model
        %     obsprop = vertcat(obsobj(:).Params);
        %     modprop = vertcat(modobj(:).Params);
        function parameterPlots(obsprop,modprop)
            %allow user to select from the parameters table and plot
            %measured against model
            pnames = obsprop.Properties.VariableNames;
            ok = 0;
            while ok<1                
                idl = listdlg("ListString",pnames,"PromptString",...
                            'Select option:','SelectionMode','single',...
                            'ListSize',[200,120],'Name','Wave spectra');
                if isempty(idl), ok = 1; continue; end

                x = modprop.(pnames{idl});
                y = obsprop.(pnames{idl});
                xylim = max(max(x),max(y));
                hf = figure('Name','Parameters','Tag','PlotFig');
                ax = axes(hf); %#ok<LAXES>
                plot(ax,x,y,'x')
                xlabel(sprintf('Modelled %s',pnames{idl}))
                ylabel(sprintf('Measured %s',pnames{idl}))
                hold(ax,'on')
                    plot(ax,[0,xylim],[0,xylim],'--k')
                hold(ax,'off')
            end
        end 
   

%%   
        function ptype = plotType()
            %define 
            ptype = questdlg('What type of plot','XY Polar','XY','Polar','XY'); 
        end

%%
        function dataFigure(src,~)
            %generate data table for data button used in bioInfluencePlot
            hf = figure('Name','Wave parameters','Tag','PlotFig');
            colnames = src.UserData.Properties.VariableNames;
            rownames = src.UserData.Properties.RowNames;
            values = table2cell(src.UserData);
            Table = uitable('Parent',hf, ...
                            'ColumnName',colnames,...
                            'RowName', rownames, ....
                            'ColumnWidth',{'auto'}, ...
                            'Data',values);  
            Table.Position(3:4) = Table.Extent(3:4);
            hf.Position(3:4) = Table.Extent(3:4)+2*Table.Position(1);
            color = get(hf,'Color');
            set(gca,'XColor',color,'YColor',color,'TickDir','out')
            title(sprintf('Results for Figure %d',src.Parent.Number))
        end   

%%
        function sd = getSubDomain(x,y,subdomain)
            %find the subdomain in integer grid indices defined by x,y range
            %subdomain defined as [x0,xN,y0,yN];    
            if isempty(subdomain) || length(subdomain)~=4
                subdomain = [min(x),max(x),min(y),max(y)];
            end
            ix0 = find(x<=subdomain(1),1,'last');
            ixN = find(x>=subdomain(2),1,'first');
            iy0 = find(y<=subdomain(3),1,'last');
            iyN = find(y>=subdomain(4),1,'first');
            sd.x = [ix0,ix0,ixN,ixN];
            sd.y = [iyN,iy0,iy0,iyN];
        end

%%
        function [newtime,ids,ido] = subSampleTime(seldst,obsdst)
            %find the times in seldst that match the times in obsdst within
            %a user defined tolerance (seconds)
            seltime = seldst.RowNames;     %dataset being sampled   
            obstime = obsdst.RowNames;     %dataset to define new times   
            inp = inputdlg({'Use a tolerance (s)? [0 for exact match]'},'Subsample',1,{'900'});
            if isempty(inp)
                tol = seconds(0);
            else
                tol = seconds(str2double(inp{1}));
            end
            [newtime,ids] = getMatchingTimes(obstime,seltime,tol);

            if numel(newtime)~=numel(obstime)
                %if there are records missing in newtime need to match up obstime
                [~,ido] = getMatchingTimes(newtime,obstime,tol);
            else
                ido = [];
            end 
            %nested function
            function [ntime,idx] = getMatchingTimes(time1,time2,tol)
                D = abs(time1 - time2');      % duration matrix
                [minDiff, idd] = min(D, [], 2);
                tf = minDiff <= tol;
                idx = idd(tf);
                ntime = time1(tf);
            end
        end

%%
        function addORupdate(mobj,cobj,dsp,newdata,meta)
            %prompt user to add or update existing record and save newdata
            % cobj - instance of class to be added to new class
            % dsp - DSproperties struct for the new variables
            % newdata - cell array of variables to be added
            % meta - struct for source and metadata
            muicat = mobj.Cases;
            answer = questdlg('New case, or Add to existing case?',...
                                  'Save grid','New','Add','Quit','Add');
            if strcmp(answer,'Quit')
                 getdialog('Data NOT saved');
                return;
            elseif strcmp(answer,'Add')                %add to existing
                caserec = caseRec(muicat,cobj.CaseIndex);
                addVariable2CaseDS(muicat,caserec,newdata,dsp);
            else                                       %new record
                datasetname = getDataSetName(cobj);
                mtime = cobj.Data.(datasetname).RowNames;
                newdst = dstable(newdata{:},'RowNames',mtime,'DSproperties',dsp);
                newdst.Source = meta.source;
                newdst.MetaData = meta.data;
                       
                %save results
                heq = str2func(metaclass(cobj).Name);
                obj = heq();  %instance of class object
                setDataSetRecord(obj,muicat,newdst,'model');
                getdialog('Data saved');      
            end
        end

%%
        function dsp = setDSroperties(option)
           %define the variables in the dataset
            %define the metadata properties for the demo data set
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]);   
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            switch option
                case 'gamma'
                dsp.Variables = struct(...
                    'Name',{'gamma'},...
                    'Description',{'JONSWP gamma'},...
                    'Unit',{'-'},...
                    'Label',{'Gamma (-)'},...
                    'QCflag',repmat({'model'},1,1)); 
            end

            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dsp.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});   
        end
    end
end