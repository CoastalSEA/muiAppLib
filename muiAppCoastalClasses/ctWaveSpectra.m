classdef ctWaveSpectra < matlab.mixin.Copyable
%
%-------class help------------------------------------------------------
% NAME
%   ctWaveSpectra.m
% PURPOSE
%   Class for analysing wave spectra data held as spectral density as a
%   function of direction and frequency, or loaded from a file
% NOTES
%   Spectral data are imported to the ctWaveData class using functions
%   such as wave_cco_spectrum. The record holds two datasets named as
%   Spectra and Properties.
% SEE ALSO
%   see ct_costal_plots and SpectralTransfer for examples of use
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2025
%--------------------------------------------------------------------------
%    
    properties 
        %struct for Spectral density matrix and dimensions (freq and dir)
        Spectrum = struct('SG',[],'freq',[],'dir',[],'date',[],'depth',[])          
        Params %table for integral properties of the spectrum
        %struct with details of input definition depending on source
        inpData        % Wave - source,Hs,Tp,Dir,ds
                       % Wind  - source,Uw,zw,Dir,Fetch,df,ds
                       % Spectrum - source,tsdst,ds,issat
                       % depth at site, ds, to include any variation in swl
        %struct with specification for spectrum model to use               
        spModel        % form - wave spectrum model
                       % source - wind or wave data
                       % spread - relationship to define directional spreading
                       % nspread - exponent for directional spreading
                       % gamma - peakiness exponent in JONSWAP spectrum
                       % depth - site depth for saturation limit
                       % inptxt - summary text for spectrum model
                               % T_gamma - dynamically set in getModelSpectrum when
                               %           T2 and Tp used to define gamma
        %struct for interpolation settings (defaults set in constructor)
        Interp         % dir - interval used to interpolate directions (deg)
                       % freq - interval used to interpolate frequency (Hz)
                       % flim - frequency limits
        %struct for plotting of wave spectra
        Plotxt         % xtxt - label for x-axis
                       % ytxt - label for y-axis
                       % vtxt - label for variable in surface plot
                       % ttxt - label for plot title
        ModelMovie
    end

    methods
        function obj = ctWaveSpectra
            %class constructor
            obj.Interp.dir = 360/512; %interval used to interpolate directions (deg)
            obj.Interp.freq = 0.005;  %interval used to interpolate frequency (Hz)
            obj.Interp.flim = [0.025,0.58]; %frequency limits
            %default plot labels
            obj.Plotxt.xtxt = 'Wave period (s)';
            obj.Plotxt.ytxt = 'Direction (degTN)';
            obj.Plotxt.vtxt = 'Spectral Energy (m^2s)';  
            obj.Plotxt.ttxt = '';
        end
    end
%% ------------------------------------------------------------------------
% Plotting and analysis call functions
%--------------------------------------------------------------------------
    methods (Static)
        function ax = getPlotOption(mobj)
            %user selected plotting options
            listxt = {'Plot a spectrum using Case data',...
                      'Plot a spectrum using a model',...
                      'Plot a bimodal model spectrum',...
                      'Plot a spectrum loaded from file',...
                      'Plot comparison of Case spectra',...
                      'Plot measured against modelled',...
                      'Animation of spectrum timeseries',...
                      'Model v Measured timeseries skill',...
                      'Bimodal model analysis',...
                      };

            selection = listdlg("ListString",listxt,"PromptString",...
                            'Select option:','SelectionMode','single',...
                            'ListSize',[220,140],'Name','Wave spectra');
            if isempty(selection), return; end
            switch selection
                case 1                  %Plot a spectrum using case data
                    ax = ctWaveSpectra.plotCaseSpectrum(mobj);
                case 2                  %Plot a spectrum using a model
                    ax = ctWaveSpectra.plotModelSpectrum();
                case 3                  %Plot a bimodal spectrum using a model
                    ax = ctWaveSpectra.plotBimodalModelSpectrum();
                case 4                  %Plot a spectrum loaded from file
                    ax = ctWaveSpectra.plotFileSpectrum();
                case 5                  %Compare case spectra
                    ax = ctWaveSpectra.plotCFCaseSpectrum(mobj);
                case 6                  %Compare measured and modelled
                    ax = ctWaveSpectra.plotCFmodelSpectrum(mobj);
                case 7                  %Timeseries animation
                    ctWaveSpectra.animateCaseSpectrum(mobj);
                case 8                  %Fit a model to a timeseries of measured spectr
                    ctWaveSpectra.fitModel2Measured(mobj);
                case 9                  %decompose measured spectrum to bimodal form
                    ctWaveSpectra.bimodalSpectrum(mobj);
            end
        end

%%
        function ax = plotCaseSpectrum(mobj)
            %Plot a spectrum using case data
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            obj = ctWaveSpectra;            
            [obj,tsdst] = getCaseInput(obj,mobj);            
            if isempty(tsdst), ax = []; return; end
            tsdst(1).DataTable = rmmissing(tsdst(1).DataTable);%remove nans
        
            %select from dates, check for NaNs, get spectrum and plot
            dates = tsdst(1).DataTable.Properties.RowNames;
            ok = 0; j = 0; ax = gobjects(0);
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                obj = getSpectrumObject(obj,tsdst,irow);
                %call plot function   
                pax = getPlot(obj,ptype);
                ax(1,j+1) = pax;                        %returns array of plot axes
            end
        end
        
%%
        function ax = plotModelSpectrum(~)
            %Plot a spectrum using a model and user defined wave/wind parameters
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');         

            ok = 0; j = 0; ax = gobjects(0);
            while ok<1    
                obj = ctWaveSpectra;
                obj = setSpectrumModel(obj);             %define the model to be used (Jonswap etc)
                if isempty(obj.spModel),return; end

                inptype = obj.spModel.source;
                obj = setForcingConditions(obj,inptype); %UI to input wave/wind conditions
                if isempty(obj), return; end    

                obj = getModelSpectrum(obj);             %compute spectrum for specified conditions
                if isempty(obj.Spectrum.SG), continue; end

                obj.Params = wave_spectrum_params(obj);  %integral properties of spectrum
                inputMessage(obj);                       %display inputs in command window   
                display(obj.Params);                     %table of combined, wind-wave and swell

                %call plot function   
                obj.Plotxt.ttxt = getModelInputText(obj);
                pax = getPlot(obj,ptype);
                ax(1,j+1) = pax;                         %returns array of plot axes
                clear obj
            end
        end

%%
        function ax = plotBimodalModelSpectrum(~)
            %Plot a bimodal spectrum using a model and user defined waved parameters
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            inptype = 'Wave';
            ok = 0; j = 0; ax = gobjects(0);
            while ok<1    
                obj = ctWaveSpectra;
                obj = setSpectrumModel(obj);             %define the model to be used (Jonswap etc)
                if isempty(obj.spModel),return; end

                obj = setForcingConditions(obj,inptype); %UI to input wave/wind conditions
                if isempty(obj), return; end   

                obj = getBimodalModelSpectrum(obj); %compute spectrum for specified conditions
                if isempty(obj.Spectrum.SG), continue; end
                
                inputMessage(obj);                       %display inputs in command window   
                display(obj.Params);                     %table of combined, wind-wave and swell
                obj.Params([2,3],:) = [];                %only use combined for plotting

                %call plot function   
                obj.Plotxt.ttxt = getModelInputText(obj);
                pax = getPlot(obj,ptype);
                ax(1,j+1) = pax;                         %returns array of plot axes
                clear obj
            end          
        end

%%
        function ax = plotFileSpectrum(~)
            %Plot a spectrum loaded from file
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            obj = ctWaveSpectra;
            
            ok = 0; j = 0; ax = gobjects(0);
            while ok<1 
                obj = getSpectrumFile(obj);        %load spectrum data from file
                if isempty(obj),return; end

                obj = getMeasuredSpectrum(obj);    %compute spectrum based on measured form
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                
                %call plot function
                obj.Plotxt.ttxt = sprintf('File: %s',matlab.lang.makeValidName(obj.inpData.file));
                pax = getPlot(obj,ptype);
                ax(1,j+1) = pax;                   %returns array of plot axes
            end
        end

%%
        function ax = plotCFmodelSpectrum(mobj)
            %plot a comparison of measured and modelled spectra
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');
            
            obj = ctWaveSpectra;            
            classops = {'ctWaveData'};
            [cobj,tsdst,dsname] = getCaseDataset(obj,mobj,classops,1);
            if ~strcmp(dsname,'Spectra')
                warndlg('Spectral data required for this option')
                ax = []; return;
            else
                tsprops = cobj.Data.Properties;
            end            

            obj.inpData.form = 'Measured';
            obj.inpData.source = 'Spectrum';           

            dates = tsdst.DataTable.Properties.RowNames;
            ok = 0; j = 0; ax = gobjects(0); isfirst = true;
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                %construct measured spectrum
                obj(1).inpData.tsdst = getDSTable(tsdst,irow,[]);     %selected record
                obj(1).inpData.tsprops = getDSTable(tsprops,irow,[]); %selected record
                obj(1).Plotxt.ttxt = sprintf('%s (%s)',tsdst.Description,dates{irow});

                obj(1) = getMeasuredSpectrum(obj(1));         %compute spectrum based on measured form
                obj(1).Params = wave_spectrum_params(obj(1)); %integral properties of spectrum
                skill = getSkillParameters(obj,mobj);
                %construct model spectrum
                obj(2) = setSpectrumModel(ctWaveSpectra);     %define the model to be used (Jonswap etc)
                obj(2) = getWaveModel(obj(2),obj(1).Params);
                if isempty(obj(2).Spectrum.SG), return; end

                %plot results
                obj(2).Plotxt.ttxt = getModelInputText(obj(2));   
                compareProperties(obj);
                ax(j+1,:) = plotObsModel(obj,ptype);
                sax = omniSpectrumPlot(obj);
                subtitle(sax,obj(2).Plotxt.ttxt)
                
                %plot skill
                metatxt = {obj(1).Plotxt.ttxt, obj(2).Plotxt.ttxt};
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
        function sax = plotCFCaseSpectrum(mobj)
            %compare spectrua selected using case data
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');
            %select from dates, check for NaNs, get spectrum and plot
            obj(2) = ctWaveSpectra;
            sax = gobjects(0);
            for i=1:2
                [obj(i),tsdst] = getCaseInput(obj(i),mobj);
                if isempty(tsdst), sax = []; return; end
                if ~isprop(tsdst,'So')
                    tsdst(1).DataTable = rmmissing(tsdst(1).DataTable);%remove nans
                end
                dates = tsdst(1).DataTable.Properties.RowNames;
                irow = listdlg("PromptString",'Select event to plot',...
                            'SelectionMode','single','ListSize',[160,300],...
                            'ListString',dates);
                if isempty(irow), return; end
                obj(i) = getSpectrumObject(obj(i),tsdst,irow);
                %call plot function   
                sax(1,i) = getPlot(obj(i),ptype);%returns array of plot axes                   
            end

            getMultiPlot(obj,sax);
        end

%%
        function animateCaseSpectrum(mobj)
            %create animation from spectral data set or wave timeseries
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            anobj = ctWaveSpectra;     
            [obj,tsdst,~] = getCaseInput(anobj,mobj);
            if isempty(tsdst), return; end

            tsdst.DataTable = rmmissing(tsdst.DataTable);   %remove nans
            tsdst = getsampleusingrange(tsdst);            

            if any(contains(tsdst.VariableNames,'Kurt'))
                obj = getMeasuredTS(anobj,tsdst);            %compute spectrum based on measured form
            else  
                dates = tsdst.DataTable.Properties.RowNames; 
                obj = getModelTS(obj,tsdst.DataTable,dates); %compute spectrum for specified conditions                   
            end

            obj(1).Plotxt.ttxt = sprintf('%s',tsdst.Description);
            wrm_single_animation(mobj,obj,ptype);
        end

%%
        function fitModel2Measured(mobj)
            %create a measured and model timeseries of spectra to determine
            %the best model fit parameters for a given sea location
            anobj = ctWaveSpectra;    
            % [obj,tsdst,plttxt,~] = getCaseInput(anobj,mobj);
            % if isempty(tsdst), return; end
            classops = {'ctWaveData'};
            [~,tsdst,dsname] = getCaseDataset(anobj,mobj,classops,1);
            if ~strcmp(dsname,'Spectra') || isempty(tsdst)
                warndlg('Spectral data required for this option')
                return; 
            end       

            tsdst.DataTable = rmmissing(tsdst.DataTable);%remove nans
            tsdst = getsampleusingrange(tsdst); 

            %get timeseries of spectra for comparison
            obsobj = getMeasuredTS(anobj,tsdst);  
            anobj = setSpectrumModel(ctWaveSpectra);   
            anobj.inpData.source = 'Spectrum';
            spectra = [obsobj(:).Spectrum]; 
            dates = [spectra(:).date];
            modobj = getModelTS(anobj,vertcat(obsobj(:).Params),dates);
            if isempty(modobj), return; end   %user cancelled
            %check on gamma values
            if isnan(anobj.spModel.gamma)
                hf = figure('Tag','PlotFig'); ax = axes(hf);
                spec = vertcat(obsobj(:).Spectrum);
                spgm = vertcat(modobj(:).spModel);
                plot(ax,[spec(:).date],[spgm(:).T_gamma],'x')
            end

            plot_spectrum_model_skill(obsobj,modobj,mobj);
            parameterPlots(obsobj,modobj);
        end

%%
        function bimodalSpectrum(mobj)
            %analyse measured spectrum for bi-modality and explore
            %representing this in a model
            obj = ctWaveSpectra;            
            [obj,tsdst] = getCaseInput(obj,mobj);
            if isempty(tsdst), return; end
            tsdst.DataTable = rmmissing(tsdst.DataTable);%remove nans
            dates = tsdst.DataTable.Properties.RowNames;

            ok = 0;
            while ok<1
                irow = listdlg("PromptString",'Select event to plot',...
                    'SelectionMode','single','ListSize',[160,300],...
                    'ListString',dates);
                if isempty(irow), ok = 1; continue; end
                obj(1).Plotxt.ttxt = sprintf('%s (%s)',tsdst.Description,dates{irow});
                obj(1).inpData.tsdst = getDSTable(tsdst,irow,[]);  %selected record
                obj(1) = getMeasuredSpectrum(obj(1));  %compute spectrum based on measured form
                freq = obj(1).Spectrum.freq;
                dir = obj(1).Spectrum.dir;
                SG = obj(1).Spectrum.SG;
                Sf = getOmniDirSpectrum(obj(1));
                idx = freq>0.05 & freq<0.2; %range for minimum search (excludes tails)
                [minSf,~] = min(Sf(idx));   %minimum spectral density in range
                idmn = find(Sf==minSf);
                [maxSf,~] = max(Sf);        %maximum spectral density

                minpeakdist = 2;         %minimum no of points separating peaks **
                minpeakht = maxSf*0.2;   %minimum height of peaks **
                [locs,~] = peakseek(Sf,minpeakdist,minpeakht);
                idlof = locs(locs<idmn);  %index of peaks less than minimum

                if isempty(idlof)
                    %not bimodal use full spectrum to estimate parameters
                    w_params = wave_spectrum_params(SG,freq,dir,0);%0=diagnostics not required
                    s_params = w_params;  s_params{:,:} = 0;       %dummy table to allow concatanation
                else
                    SGw = SG(:,idmn:end);  %wind-wave energy
                    fw  = freq(idmn:end);  %wind-wave frequencies
                    w_params = wave_spectrum_params(SGw,fw,dir,0); %0=diagnostics not required
                    SGs = SG(:,1:idmn);    %swell energy
                    fs  = freq(1:idmn);    %swell frequencies
                    s_params = wave_spectrum_params(SGs,fs,dir,0); %0=diagnostics not required
                end
                summary = [w_params;s_params];
                summary.Properties.RowNames = {'Wind-waves','Swell waves'};
                
                %can use 
                anobj = setSpectrumModel(ctWaveSpectra); 
                obj(2) = copy(anobj);
                obj(2).spModel.nspread = anobj.spModel.nspread(1);
                obj(2).spModel.gamma = anobj.spModel.gamma(1);
                if isempty(obj(2).spModel),return; end
                obj(2) = getWaveModel(obj(2),w_params);  %compute spectrum and params for specified conditions
                if isempty(obj(2).Spectrum.SG), return; end

                if isempty(idlof)
                    summary = [w_params;obj(2).Params];
                    summary.Properties.RowNames = {'Wind-waves','Model wind'};
                else
                    obj(3) = copy(anobj);
                    if length(anobj.spModel.nspread)>1
                        obj(3).spModel.nspread = anobj.spModel.nspread(2);
                    end
                    if length(anobj.spModel.gamma)>1
                        obj(3).spModel.gamma = anobj.spModel.gamma(2);  
                    end
                    if isempty(obj(3).spModel),return; end

                    obj(3) = getWaveModel(obj(3),s_params);  %compute spectrum and params for specified conditions
                    if isempty(obj(3).Spectrum.SG), return; end  
    
                    summary = [w_params;obj(2).Params;s_params;obj(3).Params];
                    summary.Properties.RowNames = {'Wind-waves','Model wind',...
                                               'Swell waves','Model swell'};
                    obj(2).Spectrum.SG = obj(2).Spectrum.SG+obj(3).Spectrum.SG; 
                    obj(3) = [];
                end
               
                display(summary)
                ax = omniSpectrumPlot(obj);
                hold(ax,'on')
                plot(ax,[1,1]*w_params.Tp,ylim,'--b')
                if ~isempty(idlof)
                    plot(ax,[1,1]*s_params.Tp,ylim,'--r')
                    plot(ax,[1,1]*1/freq(idmn),ylim,'-.g')  
                    legend({'Measured','Modelled','Sea Tp','Swell Tp',...
                        'Crossover'},'Location','northeast');
                else
                    legend({'Measured','Modelled'},'Location','northeast');
                end
                hold(ax,'off')           

                %surface plot comparison
                obj(1).Plotxt.ttxt = sprintf('%s (%s)',tsdst.Description,dates{irow});
                obj(2).Plotxt.ttxt = getModelInputText(anobj);  
                obj(2).Plotxt.vtxt = 'Modelled Spectral Energy (m^2s)';  
                plotObsModel(obj,'XY');
                obj = ctWaveSpectra; %reset to blank instance
            end
        end
    end

%% ------------------------------------------------------------------------
% Input functions
%--------------------------------------------------------------------------
    methods
        function obj = getSpectrumObject(obj,tsdst,irow)
            %get the spectrum and wave properties for selected input
            dates = tsdst(1).DataTable.Properties.RowNames;
            w_dst = getDSTable(tsdst(1),irow,[]);      %selected record
            obj.Plotxt.ttxt =  sprintf('%s: %s',obj.inpData.desc,dates{irow});

            obj.inpData.tsdst = w_dst;                 %assign data to input
            if length(tsdst)>1
                %tsdst(1) defines wind-wave. get same record for swell
                idx = find(tsdst(2).RowNames==w_dst.RowNames);
                s_dst = getDSTable(tsdst(2),idx,[]);   %selected record
                obj.inpData.tsdst(2) = s_dst; 
            end

            obj = getSpectrum(obj);  
            if isempty(obj.Spectrum.SG), return; end
            inputMessage(obj);
            display(obj.Params);

            if ~strcmp(obj.inpData.source,'Spectrum')
                obj.Plotxt.ttxt = sprintf('%s\n%s',obj.Plotxt.ttxt,...
                                                  getModelInputText(obj));                      
            end                

            obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
            %option to return diagnostics struct to investigate
            %multi-modality
            % [params,diagn] = wave_spectrum_params(obj,true); %integral properties of spectrum
            % obj.Params = params;
        end

%%
        function obj = getSpectrum(obj)
           %return measured or model spectrum depending on data inputs in
           %calling instance of ctWaveSpectra. Pass selected datetime row 
           %as a dstable or 2 x dstables if swell is included, with spectrum 
           %model already set if wave data rather than spectrum           
            ts = obj.inpData.tsdst;                  %selected data
            if strcmp(obj.inpData.source,'Spectrum') || ...
                              strcmp(obj.inpData.source,'Measured spectra') 
                obj = getMeasuredSpectrum(obj);      %compute spectrum based on measured form
            elseif isscalar(ts) %istable(ts) && numel(ts)
                obj = setInputParams(obj,ts);        %add input needed to construct spectrum
                obj = getModelSpectrum(obj);         %compute spectrum for specified conditions    
            else
                obj = setInputParams(obj,ts);        %add input needed to construct spectrum
                obj = getBimodalModelSpectrum(obj);  %compute spectrum for specified conditions                  
            end
        end

%%
        function obj = setSpectrumModel(obj)
            %set spectrum form, data source, and parameters for wave
            %spectrum and directions spreading functions           
            sp = {'JONSWAP fetch limited','TMA shallow water',...
                   'Pierson-Moskowitz fully developed','Bretschneider open ocean'};
            src = {'Wave','Wind'};
            if ~isempty(obj.inpData) && strcmp(obj.inpData.source,'Wind')
                src = {'Wind','Wave'};
            end
            
            spr = {'SPM cosine function','Donelan secant function'};
            % nsp = string(0:1:10);
            ptxt = inputHelpTxt();
            selection = inputgui('FigureTitle','Spectrum',...
                                 'InputFields',{'Wave Spectra','Data type',...
                                    'Spread function','Spread exponent',...
                                    'Jonswap gamma','Water depth'},...
                                 'Style',{'popupmenu','popupmenu','popupmenu',...
                                         'edit','edit','edit'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{sp,src,spr,'2','0','0'},...%use nspread=0 if included in wave data
                                 'PromptText',ptxt,'Position',[0.25,0.56,0.25,0.1]);
            if isempty(selection)
                spectrum = [];
            else
                spectrum.form = sp{selection{1}};
                spectrum.source = src{selection{2}};
                spectrum.spread = spr{selection{3}};
                spectrum.nspread = str2num(selection{4}); %#ok<ST2NM> handles vectors
                spectrum.gamma = str2num(selection{5});   %#ok<ST2NM>
                spectrum.depth = str2double(selection{6});
                satxt = 'excluding depth saturation';
                if strcmp(spectrum.form,'TMA shallow water') && spectrum.depth>0
                    satxt = 'including depth saturation';
                end
                spectrum.inptxt = sprintf('%s using %s (%s)\nSpreading function is %s with an exponent of %s. Jonswap gamma = %s',...
                                 sp{selection{1}},src{selection{2}},satxt,...
                                 spr{selection{3}},selection{4},selection{5});
            end 
            obj.spModel = spectrum;

            %nested function-----------------------------------------------
            function txt = inputHelpTxt()
                %txt1 = '              Select values to use';
                txt2 = 'Selection: Spread exponent should be 0 if defined in data';
                txt3 = 'Jonswap gamma only used in Jonswap and TMA models';
                txt4 = 'Setting gamma>0 overrides built-in relationships';
                txt5 = 'If input defines sea and swell, gamma and spread can accept 2 values';
                txt6 = 'Set depth to apply saturation in TMA spectrum [0= no saturation]';

                txt = sprintf('%s\n%s\n%s\n%s\n%s',txt2,txt3,txt4,txt5,txt6);
            end
        end

%%
        function obj = setDefaultSpectrumModel(obj)
            %set spectrum form, data source, and parameters for wave
            %spectrum and directions spreading functions
            spectrum.form = 'JONSWAP fetch limited';
            spectrum.source = 'Wave';
            spectrum.spread = 'SPM cosine function';
            spectrum.nspread = [2,2];
            spectrum.gamma = NaN;
            spectrum.depth = 0;
            obj.spModel = spectrum;
        end
        
%%
        function obj = setForcingConditions(obj,sptype)
            %set the user input for wave or wind conditions, or import file
            %of spectrum data (*.spt)
            if nargin<2 || isempty(sptype)
                sptype = questdlg('Wind or Wave input?','Input',...
                                                     'Wind','Wave','Wave');
            end

            switch sptype
                case 'Quit'
                    obj = []; return;

                case 'Wave'           %define spectrum using wave parameters
                    promptxt = {'Wave height (m)','Peak period (s)',...
                                     'Wave direction (degTN)'};
                    defaults = {'1.1','8.2','185'};
                    inpt = inputdlg(promptxt,'Input conditions',1,defaults);
                    if isempty(inpt), return; end  %user cancelled
                    inp.Hs = str2num(inpt{1}); %#ok<ST2NM>
                    inp.Tp = str2num(inpt{2}); %#ok<ST2NM>
                    inp.Dir = str2num(inpt{3}); %#ok<ST2NM>
                    inp.form = 'Model';

                case 'Wind'           %define spectrum using wind parameters                     
                    promptxt = {'Wind Speed (m/s)','Wind Direction (degTN)',...
                                'Height above msl (m)','Fetch Length (m)',...
                                'Depth over fetch (m)'};                                 
                    defaults = {'20.0','185','10.0','20000','0'};
                    inpt = inputdlg(promptxt,'Input conditions',1,defaults);
                    if isempty(inpt), return; end  %user cancelled
                    inp.Uw = str2double(inpt{1});
                    inp.Dir = str2double(inpt{2});
                    inp.zW = str2double(inpt{3});
                    inp.Fetch = str2double(inpt{4});
                    inp.df = str2double(inpt{5});
                    inp.form = 'Model';

                otherwise
                    getdialog('Unknown type of source in ctWaveSpectra.getForcingConditions')
                    obj = []; return; 
            end

            inp.source = sptype;
            obj.inpData = inp;
        end

%%
        function obj = getSpectrumFile(obj)
            %load spectrum from a file
            [filename,path,~] = getfiles('MultiSelect','off',...
                'FileType',{'*.spt; *.txt'},'PromptText','Select file:');
            if filename==0, obj = []; return; end  %user cancelled
            varlist = {'',[path,filename]};
            specdst = wave_cco_spectra('getData',varlist{:});
            % tsdst = horzcat(specdst.Spectra,specdst.Properties);
            % tsdst = activatedynamicprops(tsdst,specdst.Properties.VariableNames);
            % inputs.swl = 0;
            % tsdst = addvars(tsdst,inputs.swl,'NewVariableNames','swl');
            inputs.form = 'Measured';
            inputs.source = 'Spectrum';
            inputs.file = filename;
            inputs.tsdst = specdst.Spectra;
            inputs.tsprops = specdst.Properties;
            obj.inpData = inputs;
        end

%%
        function [obj,tsdst,cobj] = getCaseInput(obj,mobj)
            %get selection and sort input based on data source
            [cobj,tsdst,dsname] = getCaseDataset(obj,mobj); %returns timeseries dst
            if isempty(cobj), tsdst = []; return; end

            obj.inpData.form = 'Measured'; 
            obj.inpData.desc = sprintf('%s(%s)',tsdst.Description,dsname);
            if contains(dsname,'Spectra')            %ts data of spectra
                obj.inpData.source = 'Spectrum';
                obj.Plotxt.vtxt = 'Measured Spectral Energy (m^2s)';

            elseif strcmp(dsname,'sptProperties')     %ts data of spectra properties
                warndlg('Select ''Spectra'' dataset when using spt input')
                tsdst = [];
            else                                   %ts data of wave/wind parameters
                if isa(cobj,'ctWaveData')
                    obj.inpData.source = 'Wave';
                elseif isa(cobj,'WRM_WaveModel') 
                    obj.inpData.source = 'Wave';
                    obj.inpData.form = 'Modelled';
                elseif isa(cobj,'ctWindData')
                    obj.inpData.source = 'Wind';
                else
                    warndlg('muiUserModel not yet handled in ctWaveSpectra.setSpectrum')
                    return
                end
                obj.Plotxt.vtxt = 'Modelled Spectral Energy (m^2s)'; 
                obj = setSpectrumModel(obj);       %define the model to be used (Jonswap etc)
                if isempty(obj.spModel),return; end

                tsdst = mapInputParams(obj,tsdst); %extract required variables
            end 
        end
%%
        function [cobj,tsdst,dsname] = getCaseDataset(~,mobj,classops,idd)
            %get selection and load case. option to limit classopt in call
            if nargin<3
                idd = [];
                classops = {'ctWaveData','ctWindData','WRM_WaveModel','muiUserModel'};
            elseif nargin<4
                idd = [];                
            end

            promptxt = {'Select Case:'};
            if isempty(idd)      %user selects the datasetdsname = datasets{idd};
                [cobj,~,datasets,idd] = selectCaseDataset(mobj.Cases,[],...
                                                            classops,promptxt);
                if isempty(idd), tsdst = []; dsname = []; return; end            
                dsname = datasets{idd};                
            else                 %dataset id defined in call
                [cobj,~] = selectCaseObj(mobj.Cases,[],classops,promptxt);
                if isempty(cobj), tsdst = []; dsname = []; return; end  
                datasets = fields(cobj.Data);
                dsname = datasets{idd};
            end
            tsdst = cobj.Data.(dsname);
        end

%%
        function xtsdst = mapInputParams(obj,tsdst)
            %check for valid variable names when timeseries wave or wind
            %data are used to define conditions
            if strcmp(obj.inpData.source,'Wave')
                xtsdst = extract_wave_data(tsdst); %returns 1x2 array if bimodal
            elseif strcmp(obj.inpData.source,'Wind')
                xtsdst = extract_wind_data(tsdst,1); %isfetch=true
            end
        end

%%
        function obj = setInputParams(obj,tsdst)
            %extract input data from dstable and assign in format needed for
            %spectrum function
            sptype = obj.inpData.source;
            inp = obj.inpData;
            if strcmp(sptype,'Wave')
                for i=1:length(tsdst)
                    inp.Hs(i) = tsdst(i).Hs;
                    inp.Tp(i) = tsdst(i).Tp;
                    inp.Dir(i) = tsdst(i).Dir;
                end
            elseif strcmp(sptype,'Wind')
                inp.Uw = tsdst.AvSpeed;
                inp.Dir = tsdst.Dir;
                inp.zW = tsdst.MetaData.zW;
                inp.Fetch = tsdst.MetaData.Fetch;
            end
            obj.inpData = inp;           
        end

%%
        function inputMessage(obj)
            %write details of the input conditions to the command window
            inp = obj.inpData;
            if strcmp(inp.source,'Wave')
                if ~isfield(inp,'Hs') && isfield(inp,'tsdst')
                    inp = inp.tsdst;
                end

                for i=1:length(inp.Hs)
                    fprintf('Input parameters: Hs=%.1f m, Tp=%.1f s, Dir=%.1f dTN\n',...
                                           inp.Hs(i),inp.Tp(i),inp.Dir(i)); 
                end
            elseif strcmp(inp.source,'Wind')
                fprintf('Input parameters: U=%.1f m/s, Dir=%.1f dTN, F=%.0f m\n',...
                                            inp.Uw,inp.Dir,inp.Fetch); 
            else
                fprintf('Wave spectrum for %s\n',inp.tsdst.DataTable.Properties.RowNames{1});
            end
        end

%%
        function modeltxt = getModelInputText(obj)
            %extract the model input to define a spectrum from spModel
            spm = obj.spModel;
            spmform = split(spm.form);
            if contains(spm.form,{'Pierson-Moskowitz fully developed','Bretschneider open ocean'})
                modeltxt = sprintf('%s and %s, spread=%d ',...
                            spmform{1},spm.spread,spm.nspread);
            else
                if spm.gamma(1)==0, spm.gamma(1) = obj.inpData.gamma(1); end
                modeltxt = sprintf('%s, gamma=%.2g, and %s, spread=%d ',...
                        spmform{1},spm.gamma(1),spm.spread,spm.nspread(1));
            end

            %handle TMA depth limit
            if contains(spm.form,'TMA') && spm.depth>0
                modeltxt = sprintf('%s, d=%.1f',modeltxt,spm.depth);
            end

            %handle bimodal options
            if length(spm.gamma)>1 || length(spm.nspread)>1
                modeltxt = sprintf('%s\n',modeltxt);
            end

            if length(spm.gamma)>1
                modeltxt = sprintf('%sswell gamma=%.2g',modeltxt,spm.gamma(2));
            end

            if length(spm.nspread)>1
                modeltxt = sprintf('%s swell spread=%d',modeltxt,spm.nspread(2));
            end
        end

%%
        function [Spectrum,Properties] = saveSpectrum(obj,mobj)
            %save an array of spectra to a dstable using the spt format
            %as used to import spectra from a wave buoy
            % input is obj.Spectrum with SG, dir and freq. 
            % saved in 64 frequency intervals
            flim = obj(1).Interp.flim;
            obsfreq = [flim(1):0.005:0.1,0.11:0.01:flim(2)];            
            nvar = length(obj);
            hpw = PoolWaitbar(nvar, 'Saving spectra');
            for i=1:nvar
                stats = setSpectrum(obj(i),obsfreq);    
                varData(i,:) = varfun(@transpose,stats);
                myDatetime(i,1) = obj(i).Spectrum.date;
                input(i,:) = [obj(i).Params.Hs,obj(i).Params.T2,obj(i).Params.Sp,NaN];  %NaN is for SST
                increment(hpw);
            end
            delete(hpw) 

            varData.Properties.VariableNames = stats.Properties.VariableNames;
            input = array2table(input);

            dsp = wave_cco_spectra('setDSproperties');
            %load the results into a dstable  
            Spectrum = dstable(varData,'RowNames',myDatetime,'DSproperties',dsp.dspec); 
            Spectrum.Dimensions.freq = obsfreq;
            
            %add properties to a dstable            
            Properties =  dstable(input,'RowNames',myDatetime,'DSproperties',dsp.dsprop);
            
            if isfield(obj(1).inpData,'tsdst')
                Spectrum.Description = obj(1).inpData.tsdst(1).Description;
                Properties.Description = obj(1).inpData.tsdst(1).Description;
            end

            %save results
            if nargin>1
                dst.Spectra = Spectrum;
                dst.Properties = Properties;
                cobj = ctWaveData;
                cobj.ModelType = 'Spectrum'; 
                setDataSetRecord(cobj,mobj.Cases,dst,'spectrum');
                getdialog(sprintf('Data loaded in class: %s',classname));
            end
        end

%% ------------------------------------------------------------------------
% Spectrum construction functions
%--------------------------------------------------------------------------
        function obj = getModelSpectrum(obj)
            %construct model spectrum from input wind or wave conditions
            % obj.inpData is a struct for most sources but a table when
            % using params of a spectrum eg when called by Model v Measured 
            % timeseries skill, or Plot measured against modelled
            sp = obj.spModel;
            [dir,freq] = spectrumDimensions(obj);
    
            %directional spreading factor for selected function  
            dir0 = obj.inpData.Dir;                   %mean wave direction 
            G = directional_spreading(dir0,dir,sp.nspread(1),sp.spread);

            %add gamma and depth from obj.spModel
            obj.inpData.gamma = sp.gamma;
            %add depth if TMA depth saturation being used
            if contains(obj.spModel.form,'TMA')
                obj.inpData.ds = sp.depth;
            end

            if istable(obj.inpData)
                params = table2struct(obj.inpData);
            else
                params = obj.inpData;
            end

            %spectral energy for selected wave spectrum
            [S,gam] = wave_spectrum(sp.form,freq,params);
            if isempty(S), obj.Spectrum.SG = []; return; end

            obj.inpData.gamma = gam;  %update gamma if modified in wave_spectrum
            obj.Spectrum.SG = G*S;
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir;
            %handle depth saturation if TMA spectrum used            
            if strcmp(sp.form,'TMA shallow water') && sp.depth>0
                obj.Spectrum.depth = sp.depth;
            end

            if strcmp(obj.inpData.form,'SpectralTransfer')
                obj.Spectrum.G = G;  %used for inshore saturation depth
            end
        end

%%
        function obj = getBimodalModelSpectrum(obj)
            %construct bimodal spectrum from input wave conditions 
            w_params = obj.inpData;
            w_params.Hs(2) = []; w_params.Tp(2) = []; w_params.Dir(2) = [];
            s_params = obj.inpData;
            s_params.Hs(1) = []; s_params.Tp(1) = []; s_params.Dir(1) = [];

            windwave = copy(obj);
            windwave.spModel.nspread = obj.spModel.nspread(1);
            windwave.spModel.gamma = obj.spModel.gamma(1);
            windwave = getWaveModel(windwave,w_params);  %compute spectrum and params for specified conditions
            if isempty(windwave.Spectrum.SG), return; end   

            swellwave = copy(obj);
            
            swellwave.spModel.nspread = obj.spModel.nspread(1);
            swellwave.spModel.gamma = obj.spModel.gamma(1);
            if length(obj.spModel.nspread)>1
                swellwave.spModel.nspread = obj.spModel.nspread(2);
            end
            if length(obj.spModel.gamma)>1
                swellwave.spModel.gamma = obj.spModel.gamma(2);
            end
            swellwave = getWaveModel(swellwave,s_params);  %compute spectrum and params for specified conditions
            if isempty(swellwave.Spectrum.SG), return; end 

            obj.Spectrum.SG = windwave.Spectrum.SG+swellwave.Spectrum.SG;
            obj.Spectrum.freq = windwave.Spectrum.freq;
            obj.Spectrum.dir = windwave.Spectrum.dir;
            obj.inpData.gamma = [windwave.inpData.gamma,swellwave.inpData.gamma];

            combined = wave_spectrum_params(obj);
            summary = [combined;windwave.Params;swellwave.Params];
            summary.Properties.RowNames = {'Combined','Wind-waves','Swell waves'};
            obj.Params = summary;
        end

%%
        function obj = getMeasuredSpectrum(obj)
            %construct measured spectrum from input spectrum data
            [dir,freq] = spectrumDimensions(obj);
            % dir_int = obj.Interp.dir;  %interval used for directions (deg)
            % f_int = obj.Interp.freq;   %interval used for frequecy (Hz)
            % f_lim = obj.Interp.flim;   %range of frequency bands (Hz)
            % 
            % %direction at dir_int degree intervals
            % dir = 0:dir_int:360-dir_int;
            inp =  obj.inpData.tsdst;
            fspec = inp.Dimensions.freq';
            dirinp = {fspec,inp.Dir,inp.Spr,inp.Skew,inp.Kurt};
            G = datawell_directional_spectra(dir,false,dirinp{:}); %isplot=false       
            SG = (G.*inp.S(:)');      %force a row vector for S

            % interpolate spectrum to equal frequency intervals       
            % freq = f_lim(1):f_int:f_lim(2);  

            SGint = zeros(size(SG,1),length(freq));
            for i = 1:size(SG,1)
                SGint(i,:) = interp1(fspec, SG(i,:), freq, 'pchip', 0);
            end 

            obj.Spectrum.SG = SGint;
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir; 
            obj.Spectrum.date = obj.inpData.tsdst.RowNames;
            obj.Params = wave_spectrum_params(obj); %integral properties of spectrum  
        end 

%%
        function obj = getMeasuredTS(inobj,tsdst)
            %use a timeseries of wave spectra input data to create spectra timeseries
            nrec = height(tsdst);
            obj(nrec,1) = inobj;      
            hpw = PoolWaitbar(nrec, 'Processing measured timeseries');  %and increment(hpw);
            parfor i=1:nrec
                anobj = copy(inobj);
                anobj.inpData.tsdst = getDSTable(tsdst,i,[]); %selected record
                anobj = getMeasuredSpectrum(anobj);           %compute spectrum based on measured form
                anobj.Params = wave_spectrum_params(anobj);   %integral properties of spectrum
                obj(i,1) = anobj;
                increment(hpw);
            end
            delete(hpw)
        end

%%
        function mod_obj = getModelTS(inobj,params,dates)
            %use a timeseries of wave data to create a timeseries of spectra  
           %dates = tsdst.DataTable.Properties.RowNames;                
            nrec = length(dates);
            mod_obj(nrec,1) = inobj;
            hpw = PoolWaitbar(nrec, 'Processing measured timeseries');  %and increment(hpw);
            for i=1:nrec
                anobj = copy(inobj);
                anobj.inpData = params(i,:);
                anobj.inpData.source = 'Wave';
                anobj = getModelSpectrum(anobj);            %compute spectrum for specified conditions
                anobj.Spectrum.date = dates(i);             %date as text string
                if isempty(anobj.Spectrum.SG), continue; end
                anobj.Params = wave_spectrum_params(anobj); %integral properties of spectrum
                mod_obj(i,1) = anobj; 
                increment(hpw);
            end
            delete(hpw)
        end

%%
        function obj = getWaveModel(obj,params)
            %construct a model wave spectrum and return with wave parameters
            obj.inpData = params;
            obj.inpData.source = 'Wave';
            obj.inpData.form = 'Model params table';
            obj = getModelSpectrum(obj);            %compute spectrum for specified conditions
            if isempty(obj.Spectrum.SG), return; end
            obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
        end

%%
        function Sf = getOmniDirSpectrum(obj)
            %extract the omni-directional spectrum from a directional spectrum
            dir = obj.Spectrum.dir;
            SG = obj.Spectrum.SG;
            theta = deg2rad(dir(:));                         %radians and column vector  
            w = trapz_weights_periodic(theta);
            Sf = sum(SG.*w,1); 
        end


%%
        function phi = kit_limit(~,f,ds)
            % Calculate the Kitaigorodskii limit to the spectrum at frquency f
            % f - wave frequency (1/s)
            % ds - water depth at site (m)
            % phi - frquency dependent Kitaigorodskii adjustment to be applied to the 
            % JONSWAP spectrum to take accound of depth limiting effects
            % using Thompson and Vincent, 1983 approximation as given in
            % Hughes, 1984, Eq.(13) and (15)   
            % <duplicated in wave_spectrum.m>
            g = 9.81;
            omega = 2*pi*f.*sqrt(ds/g);
            % if omega<=1 (vectorised) or omega>2
            phi = 0.5*omega.^2.*(omega<=1) + 1.*(omega>2);   
            phi = phi + (1 - 0.5*(2-omega).^2).*(omega>1 & omega<=2);
        end

%% ------------------------------------------------------------------------
% Plotting functions
%--------------------------------------------------------------------------
        function ax = getPlot(obj,ptype,ax)
            %call plot function   
            if nargin<3, ax = []; end

            if strcmp(ptype,'XY')                %single x-y plot
                ax = surfPlot(obj,ax);
            else
                ax = polarPlot(obj,ax);
            end  
        end

%%
        function ax = surfPlot(obj,ax) 
            %plot spectral surface as XY function of frequency and direction
            if nargin<2 || isempty(ax)
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end
            %extract data
            ptxt = obj.Plotxt;
            period = 1./obj.Spectrum.freq;
            dir = obj.Spectrum.dir; 
            
            % dir2dir0 = dir-dir0;
            SG = obj.Spectrum.SG;
            SGpk = max(SG,[],'all');
            %make plot
            surf(ax,period,dir,SG,'Tag','PlotFigSurface');
            view(2);
            shading interp
            axis tight
            %ax.XLim(2) = obj.Params.Tp*2;
            if ~isempty(obj.Params)
                hold on
                plot3(ax,xlim,obj.Params.Dir*[1 1],SGpk*[1,1],'--w','Tag','DirPk')
                plot3(ax,obj.Params.Tp*[1 1],ylim,SGpk*[1,1],'-.w','Tag','TpPk')
                hold off  
            end
            xlabel(ptxt.xtxt)
            ylabel(ptxt.ytxt)

            %add the colorbar and labels   
            cb = colorbar;
            cb.Label.String = ptxt.vtxt;
            cb(1).Tag = 'WaveSpectrum';            

            title(ptxt.ttxt)
            if ~isempty(obj.Params)
                p = obj.Params;
                subtitle(sprintf('Hs=%.2f m; Tp=%.1f s; T_2=%.1f s; Dir=%.3g degTN',...
                                                    p.Hs,p.Tp,p.T2,p.Dir));
            end
        end

%%
        function ax = polarPlot(obj,ax)
            %plot spectral surface as polar function of frequency and direction
            if nargin<2 || isempty(ax)
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end

            %extract data
            ptxt = obj.Plotxt;
            period = 1./obj.Spectrum.freq;
            dir = obj.Spectrum.dir;
            SG = obj.Spectrum.SG;

            %make plot
            wid = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
            maxT = round(obj.Params.Tp*2);
            radrange = [1,maxT];
            %interpolate var(phi,T) onto plot domain defined by tints,rints
            tints = linspace(0,2*pi,360);  
            rints = linspace(1,maxT,30);
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
            p = obj.Params;
            subtitle(sprintf('Hs=%.2fm; Tp=%.1fs; T_2=%.1fs; Dir=%.3gdegTN',...
                                        p.Hs,p.Tp,p.T2,p.Dir));
        end

%%
        function ax = omniSpectrumPlot(obj,ax)
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
        function s = plotObsModel(obj,ptype)
            %plot offshore and inshore spectra  - similar to off_in_plot in SpecralTransfer
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
            getPlot(obj(1),ptype,s(1));
            s(1).Title.String = 'Measured'; 
            s(1).CLim(2) = Smax;

            s(2) = subplot(m,n,2);
            getPlot(obj(2),ptype,s(2));
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
                getPlot(diffobj,ptype,s(3));  
                s(3).Title.String = 'Difference [Meas-Model]';
                s(3).Subtitle.String = [];
                s(3).XLim = s(1).XLim;
            end

            sgtitle(sprintf('Spectrum for: %s',obj(1).Plotxt.ttxt))
        end 

%%
        function parameterPlots(obsobj,modobj)
            %allow user to select from the parameters table and plot
            %measured against model
            obsprop = vertcat(obsobj(:).Params);
            modprop = vertcat(modobj(:).Params);
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
        function getMultiPlot(~,sax)
            %plot a set of axes in a single figure
            hf = figure('Name','SpecTrans','Tag','PlotFig');
            nplot = length(sax);
            t = tiledlayout(hf, nplot, 1); % nplot rows, 1 column
            hfigs = gobjects(0);
            for i=1:nplot
                % Move existing axes into tiles
                hfigs(i) = sax(i).Parent;
                tt = nexttile(t, i);
                delete(tt)           % removes the placeholder axes
                sax(i).Parent = t;   % Reparent to tiledlayout in figure
                sax(i).Layout.Tile = i;
            end
            delete(hfigs)
        end

%%
        function compareProperties(obj)
            %write table of values to the command window
            inp = table2struct(obj(1).inpData.tsprops.DataTable);
            m0 = (inp.Hs/4)^2;
            [~,idx] = max(obj(1).inpData.tsdst.S);        
            inp.Dp = obj(1).inpData.tsdst.Dir(idx);
            inp.Tp = 1/obj(1).inpData.tsdst.Dimensions.freq(idx);
            source = table(inp.Hs,m0,NaN,inp.Smax,inp.Tp,inp.Dp,NaN,NaN,NaN,inp.Tz); %****Smax is now Sp!!!!
            source.Properties.VariableNames = {'Hs','m0','Dir','Sp','Tp','Dp',...
                                                    'Sfdpk','Dfdpk','Tfdpk','T2'};
            out = [source;obj(1).Params;obj(2).Params];
            out.Properties.RowNames = {'Source properties','Source spectrum',...
                                                 'Model spectrum'};
            display(out)
        end
    end

%% ------------------------------------------------------------------------
% Skill functions - NB not wave spectra specific - should not be here???
%--------------------------------------------------------------------------    
    methods
        function skill = getSkillParameters(obj,mobj)
            %extract Skill parameters using MS_RunParams class for input
            x = obj(1).Spectrum.dir;              %direction
            y = obj(1).Spectrum.freq;             %frequency
            robj = MS_RunParams.setInput(mobj);   %default or current values if user cancels
        
            skill.Ro = robj.maxcorr;
            skill.n  = robj.skillexponent;
            skill.Inc = true;                   %flag to include skill score
            skill.W = robj.skillwindow;
            subdomain = robj.skillsubdomain;
            skill.SD = getSubDomain(obj,x,y,subdomain);
            skill.iter = robj.skilliteration;
        end

%%
        function sd = getSubDomain(~,x,y,subdomain)
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

%% ------------------------------------------------------------------------
% Utility functions
%--------------------------------------------------------------------------
        function [dir,freq] = spectrumDimensions(obj)
            %return default direction and frequency intervals
            dir_int = obj.Interp.dir;  %interval used for directions (deg)
            f_int = obj.Interp.freq;   %interval used for frequecy (Hz)
            f_lim = obj.Interp.flim;   %range of frequency bands (Hz)    
            %direction at dir_int degree intervals
            dir = 0:dir_int:360-dir_int;
            % interpolate spectrum to equal frequency intervals       
            freq = f_lim(1):f_int:f_lim(2);
        end

%%
        function [ndir,nfreq] = spectrumNumIntervals(obj)
            %return the number of intervals for direction and frequency
            [dir,freq] = spectrumDimensions(obj);
            ndir = length(dir);
            nfreq = length(freq);
        end

%%
        function [spectra,params] = unpackSpectrum(inobj,offobj)
            %unpack the spectrum property as a set of arrays
            % 1 - input includes inobj and offobj: return spectra and params
            % 2 - input just inobj: returns params
            nrec = length(inobj);
            hpw = PoolWaitbar(nrec, 'Unpacking spectra');
            if nargin==2
                parfor i=1:nrec
                    time(i,1) = offobj(i).Spectrum.date;    
                    swl(i,1) = offobj(i).inpData.tsdst(1).DataTable.swl;
                    Sot(i,:,:) = offobj(i).Spectrum.SG;
                    Sit(i,:,:) = inobj(i).Spectrum.SG;
                    depths(i,1) = inobj(i).Spectrum.depth;  
                    params(i,:) = inobj(i).Params;
                    increment(hpw);
                end   
                spectra = struct('time',time,'swl',swl,'Sot',Sot,'Sit',Sit,'depths',depths);
            else
                parfor i=1:nrec
                    params(i,:) = inobj(i).Params;
                    increment(hpw);
                end
                spectra = [];
            end
            delete(hpw)
        end
    end
end