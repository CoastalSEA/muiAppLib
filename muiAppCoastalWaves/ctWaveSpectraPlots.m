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
                obj = getSpectrumObject(obj,meta,tsdst,irow);
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

                obj = setSpectrumModel(obj);            %define model if not measured spectrum
                if isempty(obj.spModel), return; end    %user cancelled

                obj.inpData.gamma = obj.spModel.gamma;
                obj = getSpectrum(obj,[]);        %define the model to be used (Jonswap etc)
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
            obj = []; 
            %get the measured wave spectrum to be modelled
            getdialog('Select wave spectra data'); 
            [~,tsdst,meta] = waveModels.getCaseInputParams(mobj,{'ctWaveSpectrumData'},1); 
            if isempty(tsdst), getdialog('Select wave spectra data'); return; end
            if ~contains(meta.inptype,'Spectrum')
                warndlg('Measured spectrum required for this option');
                return;
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
                obj = getSpectrumObject(obj,meta,tsdst,irow);

                %add spectrum input properties
                obj.inpData.properties = propstable(irow,:);                
                %construct model spectrum
                obj(2) = setSpectrumModel(ctWaveSpectraPlots); %define the model to be used (Jonswap etc)
                indst = dstable(obj(1).Params,'RowNames',dates(irow));
                obj(2) = setInputParams(obj(2),indst,'Wave');
                obj(2) = getMultiModalSpectrum(obj(2),indst);
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
            if isempty(tsdst(1)), return; end           %invalid selection

            anobj = ctWaveSpectraPlots();
            if strcmp(meta.inptype,'Spectrum')
                obj = getMeasuredTS(anobj,tsdst(1));    %compute spectrum based on measured form
                obj(1).Plotxt.ttxt = tsdst(1).Description;
            else  
                obj = getModelTS(anobj,tsdst,meta);  %compute spectrum for specified conditions
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
        function bimodalSpectrum(mobj)
            %analyse measured spectrum for bi-modality and explore
            %representing this in a model
            %get the measured wave spectrum to be modelled
            getdialog('Select wave spectra data'); 
            [~,obsdst,meta] = waveModels.getCaseInputParams(mobj,{'ctWaveSpectrumData'},1);
            if isempty(obsdst), getdialog('Select wave spectra data'); return; end
            obsdst = obsdst(1);                %assign sptSpectrum as tsdst

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
                anobj = setSpectrumModel(anobj);
                if isempty(anobj.spModel), return; end    %user cancelled

                anobj = setInputParams(anobj,params,'Wave');
                anobj = getSpectrum(anobj,[]); 
                if isempty(anobj.Spectrum.SG), return; end
                obj(2) = anobj;

                if isempty(idmn)
                    summary = [w_params;obj(2).Params(2,:)];
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
            getdialog('Select wave spectra data'); 
            [~,obsdst,meta] = waveModels.getCaseInputParams(mobj,{'ctWaveSpectrumData'},1);
            if isempty(obsdst), getdialog('No wave spectra data selected'); return; end
            obsdst = obsdst(1);                %assign sptSpectrum as tsdst

            ok = 0;  
            while ok<1                         %option to limit record length             
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
                getdialog('Select wave data to compare with spectra'); 
                [~,seldst,meta] = waveModels.getCaseInputParams(mobj);
                if isempty(seldst), return; end
                hwb = waitbar(0,'Loading data');
                %match record to measured dataset
                [newtime,ids,ido] = ctWaveSpectraPlots.subSampleTime(seldst,obsdst);
                if ~isempty(ido)
                    subtable = obsdst.DataTable(ido,:); 
                    obsdst.DataTable = subtable;
                    obsdst.RowNames = newtime;
                end
                waitbar(1,hwb)
                moddst = copy(seldst);
                for j=1:numel(moddst)
                    subtable = seldst(j).DataTable(ids,:);  
                    moddst(j).DataTable = subtable;
                    moddst(j).RowNames = newtime; %update times in case there is an offset
                end
                meta.model{1} = moddst(1).Description;
                delete(hwb)
            else
                moddst = [];   %needed for parfor loop
                meta.model{1} = 'Model using spectrum properties';
                meta.variables.inputs = {'spectrum Hs','spectrum Tp','spectrum Dir'};
                meta.variables.seastate = [];
            end            

            %running long timeseries of spectra using getMeasuredTS and 
            %getModelTS to create ctWaveSpectrum objects and then computing 
            %the statistice on the object arrays can be limited by memory. 
            %To avoid this compute statistics for one spectrum at a time.
            obj = ctWaveSpectraPlots;
            obj = setSpectrumModel(obj);
            if isempty(obj.spModel), return; end % user cancelled
            [dir,freq] = spectrumDimensions(obj);
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir;

            %get the parameters for skill model - cancel uses default values
            skill = getSkillParameters(obj,mobj);        
            %get time and source meta data
            mtime = obsdst.RowNames;
            meta.source{1} = obsdst.Description;
            meta.source{2} = var2range(obsdst.RowRange);
            meta.model{2} = obj.spModel.form;
            meta.model{3} = {obj.spModel.gamma,obj.spModel.nspread};

            nrec = numel(mtime);
            seas = zeros(nrec,1);                         %needed in parfor to match loop size
            if ~isempty(meta.variables) && ~isempty(meta.variables.seastate)...
                                   && isa(meta.variables.seastate,'dstable')
                seas = meta.variables.seastate.DataTable; %required for parfor loop
            end
            gamma = [];

            hpw = PoolWaitbar(nrec, 'Processing skill statistics');
            tic
            parfor i=1:nrec                                   %parfor loop
                itsdst = getDSTable(obsdst,i,[]);             %selected record
                obsobj = copy(obj);
                obsobj = setInputParams(obsobj,itsdst,'Spectrum');            
                obsobj = getMeasuredSpectrum(obsobj);         %compute spectrum based on measured form
                obsprops(i,:) = wave_spectrum_params(obsobj); %integral properties of spectrum  

                modobj = copy(obj);
                if iscase
                    indst = ctWaveSpectrum.getDatasetRow(moddst,i); %selected record
                    seastate = seas(i,:);    
                else
                    indst = dstable(obsprops(i,:),'RowNames',mtime(i));
                    seastate = indst.DataTable; %contains T1,T2 etc used to estimate gamma
                end 
                modobj = setInputParams(modobj,indst,'Wave');
                modobj = getMultiModalSpectrum(modobj,seastate);                  
                params = modobj.Params(1,:);    
                params.Properties.RowNames = {};

                modprops(i,:) = params;
                if contains(modobj.spModel.form,'JONSWAP') || ...
                                        contains(modobj.spModel.form,'TMA')
                    gamma(i) = modobj.inpData.gamma(1);                
                end
                stats(i) = get_spectrum_skill_stats(obsobj,modobj,skill);                
                increment(hpw);
            end
            delete(hpw)
            elapsedTime = toc/60;  % Stop timer
            fprintf('Run time for %d steps: %.2f minutes\n',nrec,elapsedTime);
            %plot model skill and allow user to examine individual parameters
            ctWaveSpectraPlots.plotSpectrumModelSkill(stats,skill,meta)
            ctWaveSpectraPlots.parameterPlots(obsprops,modprops);

            if (contains(obj.spModel.form,'JONSWAP') || ...
                        contains(obj.spModel.form,'TMA')) && ...
                                obj.spModel.gamma<=0   
                if ~all(gamma==gamma(1),'all') %don't plot if constant
                    hf = figure('Tag','PlotFig'); ax = axes(hf);
                    plot(ax,mtime,gamma,'.')
                    xlabel('Time')
                    ylabel('Gamma')
                    title('Gamma values used for model')
                end
            end
        end

%%
        function estimateSpectrumGamma(mobj)
            %estimate JONSWAP gamma from spectra timeseries
            getdialog('Select wave spectra data'); 
            [cobj,obsdst,meta] = waveModels.getCaseInputParams(mobj,{'ctWaveSpectrumData'},1);
            if isempty(obsdst), getdialog('Select wave spectra data'); return; end
            specdst = obsdst(1);                %assign sptSpectrum as tsdst

            specdst = getsampleusingrange(specdst); 
            mtime = specdst.RowNames;
            propdst = getsampleusingtime(obsdst(2),mtime);
            
            nrec = height(specdst);
            gamma = zeros(nrec,1);
            gammalog = gamma; 
            g_bnd = [0.1,5];                       %bounds used for gamma search
            S = specdst.S;
            f = specdst.Dimensions.freq;  
            [~,idx] = max(S,[],2);  
            fp = f(idx); %peak frequency   

            hpw = PoolWaitbar(nrec, 'Processing timeseries');  %and increment(hpw);
            parfor i=1:nrec                        %parfor loop
                %compute gamma using linear (false) and log (true) misfit                              
                gamma(i) = wave_spectrum_gamma(S(i,:),f,g_bnd,[],false);
                gammalog(i) = wave_spectrum_gamma(S(i,:),f,g_bnd,[],true);
                increment(hpw);
            end
            delete(hpw)

            %save results if required
            answer = questdlg('Save results?','Gamma','Yes','No','No');
            if strcmp(answer,'Yes')
                meta.source = sprintf('Gamma using %s',specdst.Description);
                meta.data = sprintf('Peak period = %.1fs',1/fp);
                dsp = ctWaveSpectraPlots.setDSroperties('gamma');
                ctWaveSpectraPlots.addORupdate(mobj,cobj,dsp,{gamma},meta);
            end

            ok = 0;
            while ok<1
                %get wave height threshold and create plot
                defaults = {'0','0'};
                answer = inputdlg({'Wave height threshold:','Plot Log gamma (0/1)'},...
                                                   'Gamma',1,defaults);
                if isempty(answer), answer = defaults; end
                Hthr = str2double(answer{1});
                isplt = logical(str2double(answer{2}));
                ctWaveSpectraPlots.getGammaPlots(propdst,gamma,gammalog,fp,Hthr,g_bnd,isplt)
                answer = questdlg('Plot again?','Gamma','Yes','No','No');
                if strcmp(answer,'No'), ok = 1; end
            end
        end
%%
        function getGammaPlots(propdst,gamma,gammalog,fp,Hthr,g_bnd,isplt)
            %various plots of the fitted gamma
            Hs = propdst.Hs;
            idx = Hs<=Hthr;
            Hs(idx) = NaN;
            gamma(idx) = NaN;
            gammalog(idx) = NaN;
            nrec = sum(~isnan(gamma));
            % T = propdst.Tz;  %mean zero-crossing period
            T = (1./fp);         %peak period
            steep = 2*pi*Hs./(9.81*T.^2);

            %gammalog out of bounds assigned NaN 
            idl = gammalog>g_bnd(2)-0.05;
            gammalog(idl) = NaN;
            % gammalog(idl) = gamma(idl); %alternative to assign gamma values

            mn_gamma = mean(gamma,'all','omitnan');
            sd_gamma = std(gamma,[],'omitnan');
            mn_gammalog = mean(gammalog,'all','omitnan');
            mn_Hs = mean(Hs,'all','omitnan');
            sd_Hs = std(Hs,[],'omitnan');

            if nrec<1000; msz = 6; elseif nrec<10000, msz = 4; else, msz = 2; end
            %plot of gamma and Hs as timeseries
            hf = figure('Name','Gamma','Tag','PlotFig');
            ax = axes(hf); 
            yyaxis left
            stem(ax,propdst.RowNames,gamma,'.','Color',[0.8,0.8,0.8],...
                'MarkerEdgeColor',[0,0.45,0.74],'MarkerSize',msz,'LineWidth',msz/5) %note smallest line width is ~0.5 = 1 pixel
            ylabel('Gamma (-)')
            xlabel('Time')
            yyaxis right            
            plot(ax,propdst.RowNames,Hs,'.','MarkerSize',msz)
            ylabel('Wave height (m)')
            glegtxt = sprintf('gamma with mean %.2f, std %.2f',mn_gamma,sd_gamma);
            hlegtxt = sprintf('Hs with mean %.2f, std %.2f',mn_Hs,sd_Hs);
            legend({glegtxt,hlegtxt})
            title(propdst.Description)
            subtitle(sprintf('Gamma and wave height (Hthr=%.2fm; N=%d)',Hthr,nrec));
       
            %plot gamma against wave steepness
            model = regression_selection();
            if isempty(model), model = 'Linear'; end
            ttxt = sprintf('%s for\ngamma and wave steepness (Hthr=%.2fm; N=%d)',...
                            propdst.Description,Hthr,nrec);
            metatxt = {'Wave steepness','Gamma',ttxt};
            ax = regression_plot(steep,gamma,metatxt,model);
            hs = findobj(ax.Children,'Tag','DataPoints');            
            hs.Marker = '.'; hs.MarkerSize = msz;

            if isplt
                %plot linear gamma & log gamma (estimated using log(S)) v time
                hf = figure('Name','Gamma','Tag','PlotFig');
                ax = axes(hf); 
                stem(ax,propdst.RowNames,gamma,'.','Color',[0.8,0.8,0.8],...
                    'MarkerEdgeColor',[0,0.45,0.74],'MarkerSize',msz,'LineWidth',msz/5) %note smallest line width is ~0.5 = 1 pixel
                ylabel('Gamma (-)')
                xlabel('Time')     
                hold on
                plot(ax,propdst.RowNames,gammalog,'.','MarkerSize',msz)
                hold off
                glegtxt = sprintf('gamma with mean %.2f',mn_gamma);
                hlegtxt = sprintf('gamma log(S) with mean %.2f',mn_gammalog);
                legend({glegtxt,hlegtxt})
                title(propdst.Description)
                subtitle('Linear and Log gamma estimates')
    
                %plot of log gamma as a function of linear gamma
                % hf = figure('Name','Gamma','Tag','PlotFig');
                % ax = axes(hf);
                % plot(ax,gamma,gammalog,'.')
                txt.xlabel = 'gamma using S';
                txt.ylabel = 'gamma using log(S)';
                txt.title = propdst.Description;
                txt.subtitle = 'Log gamma(Linear gamma)';
                histogram_plot(gamma,gammalog,txt);   
            end
        end

%%
        function subsampleSpectrum(mobj)
            %get spectra Case dataset to be used
            getdialog('Select wave spectra data'); 
            [cobj,dst,~] = waveModels.getCaseInputParams(mobj,{'ctWaveSpectrumData'},1);            
            if isempty(dst), return; end
            %method can be interp1 method or none. 'none' finds exact match 
            %or match with a tolerance if tolerance>0 seconds. 
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
            getdialog(sprintf('Subsampled dataset saved as %s',classname));
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
            source = table(inp.Hs,m0,NaN,inp.Sp,inp.Tp,inp.Dp,NaN,NaN,NaN,NaN,inp.Tz,NaN); 
            source.Properties.VariableNames = {'Hs','m0','Dir','Sp','Tp','Dp',...
                              'Sfdpk','Tfdpk','Dfdpk','T1','T2','T10'};
            out = [source;obj(1).Params;obj(2).Params];
            out.Properties.RowNames = {'Source properties','Source spectrum',...
                                                 'Model spectrum'};
            display(out)
        end

%%
        function addDataButton(obj,hf,outtable)
            %add button to allow user to view summary statistics
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
            %compute the skill of model v measured spectrum data and produce 
            %Taylor plot of timeseries results
            hw = waitbar(0,'Preparing plot please wait');
            ndteststd = [stats(:).teststd]./[stats(:).refstd]; %normalised std
            rLim = ceil(max(ndteststd));                       %radial limit for the plot
            if rLim>6, rLim = 6; end
            ax = taylor_plot_figure(rLim);    
            metatxt = {'Measured','Model'};
            %local skill is not plotted even if computed but is reported in
            %the table of results on the Case list button
            waitbar(0.1,hw)
            ax = taylor_plot_ts(ax,stats,skill,metatxt); 
            waitbar(0.9,hw)
            ttxt = meta.model{1};
            if meta.iselvar
                nvar = size(meta.variables.selection,1);
                ttxt = sprintf('%s using %d sea states',ttxt,nvar);
            end
            waitbar(1,hw)
            ax.Title.String = ttxt;
            ax.Subtitle.String = sprintf('Using %s as reference',meta.source{1});

            %add button to display plot metadata
            uicontrol('Parent',ax.Parent,'Style','pushbutton',...
                'String','Sources','Tag','FigButton',...
                'TooltipString','Access selection meta data',...
                'Units','normalized','Position',[0.88 0.96 0.10 0.04],...
                'UserData',meta,...
                'Callback',@(src,evdat)ctWaveSpectraPlots.metaFigure(src,evdat)); 

            delete(hw)
        end 

%%
        function metaFigure(src,~)
            %add button to display plot meta data
            % hf = figure('Name','Data selection','Tag','PlotFig','Visible','off');
            meta = src.UserData;

            colnames = {'Selection'};
            rownames = {'Spectrum Case';'Range';'Model Case';'Model Form';...
                'Model Variables 1';'Model Variables 2';'Model Variables 3'};
            nvar = size(meta.variables.inputs,1);
            vartxt = cell(nvar+1,1);
            for i=1:nvar
                vartxt{i} = sprintf('%s, %s, %s',meta.variables.inputs{i,:});
            end
            vartxt{end} = '-';
            
            modelform = sprintf('%s (gamm0=%.2f, N=%d)',meta.model{2},...
                                                        meta.model{3}{:});
            values = {meta.source{1};meta.source{2};meta.model{1};modelform};
            values = [values;vartxt];
            rownames = [rownames(1:numel(values)-1,1);{' - '}];
            figtitle = sprintf('Data selection for Figure %d',src.Parent.Number);
            [~,~,ht] = tablefigure(figtitle,[],rownames,colnames,values);
             ht.ColumnWidth{1} = 3*ht.ColumnWidth{1};
        end

%%
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
            %generate data table for data button used in Spectrum cf plots
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
            % cobj - instance of class to be added to or used to defined
            %        time for new instance
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
                dsprop = dsproperties(dsp);
                addVariable2CaseDS(muicat,caserec,newdata,dsprop);
            else                                       %new record
                datasetname = 'sptSpectrum';
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