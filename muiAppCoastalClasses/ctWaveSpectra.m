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
        Spectrum = struct('SG',[],'freq',[],'dir',[],'date',[],'depi',[])          
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
        %struct for interpolation settings (defaults set in constructor)
        Interp         % dir - interval used to interpolate directions (deg)
                       % freq - interval used to interpolate frequency (Hz)
                       %flim - frequency limits
        ModelMovie
    end

    methods
        function obj = ctWaveSpectra
            %class constructor
            obj.Interp.dir = 360/512; %interval used to interpolate directions (deg)
            obj.Interp.freq = 0.005;  %interval used to interpolate frequency (Hz)
            obj.Interp.flim = [0.025,0.58]; %frequency limits
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
                      'Plot a spectrum loaded from file',...
                      'Plot measured against modelled',...
                      'Timeseries animation',...
                      'Model to measured timeseries skill',...
                      'Bimodal model analysis',...
                      };

            selection = listdlg("ListString",listxt,"PromptString",...
                            'Select option:','SelectionMode','single',...
                            'ListSize',[200,120],'Name','Wave spectra');
            if isempty(selection), return; end
            switch selection
                case 1                  %Plot a spectrum using case data
                    ax = ctWaveSpectra.plotCaseSpectrum(mobj);
                case 2                  %Plot a spectrum using a model
                    ax = ctWaveSpectra.plotModelSpectrum();
                case 3                  %Plot a spectrum loaded from file
                    ax = ctWaveSpectra.plotFileSpectrum();
                case 4                  %Compare measured and modelled
                    ax = ctWaveSpectra.plotCFmodelSpectrum(mobj);
                case 5                  %Timeseries animation
                    ctWaveSpectra.animateCaseSpectrum(mobj);
                case 6
                    ctWaveSpectra.fitModel2Measured(mobj);
                case 7
                    ctWaveSpectra.bimodalSpectrum(mobj);
            end
        end

%%
        function ax = plotCaseSpectrum(mobj)
            %Plot a spectrum using case data
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            obj = ctWaveSpectra;            
            [obj,tsdst,plttxt] = getCaseInput(obj,mobj);
            if isempty(tsdst), ax = []; return; end
            tsdst.DataTable = rmmissing(tsdst.DataTable);%remove nans
        
            %select from dates, check for NaNs, get spectrum and plot
            dates = tsdst.DataTable.Properties.RowNames;
            ok = 0; j = 0; ax = gobjects(0);
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                ts = getDSTable(tsdst,irow,[]);         %selected record
                plttxt{4} = sprintf('%s (%s)',tsdst.Description,dates{irow});

                obj.inpData.tsdst = ts;                 %assign data to input

                if strcmp(obj.inpData.source,'Spectrum')
                    obj = getMeasuredSpectrum(obj);     %compute spectrum based on measured form
                else     
                    obj = setInputParams(obj,ts);       %add input needed to construct spectrum
                    obj = getModelSpectrum(obj);        %compute spectrum for specified conditions
                    if isempty(obj.Spectrum.SG), continue; end
                    spm = obj.spModel;
                    spmform = split(spm.form);
                    plttxt{4} = sprintf('%s, gamma=%.2g, and %s, n=%d ',...
                                spmform{1},spm.gamma,spm.spread,spm.nspread);
                    if contains(spm.form,'TMA') && spm.depth>0
                        plttxt{4} = sprintf('%s, d=%.1f',plttxt{4},spm.depth);
                    end
                    inputMessage(obj);                    
                end
                
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                %option to return diagnostics struct to investigate
                %multi-modality
                % [params,diagn] = wave_spectrum_params(obj,true); %integral properties of spectrum
                % obj.Params = params;
                
                %call plot function   
                pax = getPlot(obj,plttxt,ptype);
                ax(1,j+1) = pax;                        %returns array of plot axes
            end
        end
        
%%
        function ax = plotModelSpectrum(plttxt)
            %Plot a spectrum using a model and user defined wave/wind parameters
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');
            
            if nargin<1
                xtxt = 'Wave period (s)';
                ytxt = 'Direction (degTN)';
                plttxt = {'Modelled Spectral Energy (m^2s)',xtxt,ytxt,'model'};   
            end
            inptype = questdlg('Wind or Wave input?',...
                                            'Input','Wind','Wave','Wave');            

            ok = 0; j = 0; ax = gobjects(0);
            while ok<1    
                obj = ctWaveSpectra;
                obj = setSpectrumModel(obj);             %define the model to be used (Jonswap etc)
                if isempty(obj.spModel),return; end

                obj = setForcingConditions(obj,inptype); %UI to input wave/wind conditions
                if isempty(obj), return; end
                %obj = setSprectrumInput(obj,0);          %add input needed for wave_spectrum.m, 0=data     
                obj = getModelSpectrum(obj);             %compute spectrum for specified conditions
                if isempty(obj.Spectrum.SG), continue; end
                obj.Params = wave_spectrum_params(obj);  %integral properties of spectrum
                inp = obj.inpData;
                if strcmp(inptype,'Wave')
                    fprintf('Input parameters: Hs=%.1f m, Tp=%.1f s, Dir=%.1f dTN\n',...
                                                    inp.Hs,inp.Tp,inp.Dir);
                else
                    fprintf('Input parameters: Uw=%.1f m/s, Dir=%.1f dTN, Fetch=%d m\n',...
                                           inp.Uw,inp.Dir,inp.Fetch);
                end
                %call plot function   
                plttxt{4} = getModelInputText(obj,1);
                pax = getPlot(obj,plttxt,ptype);
                ax(1,j+1) = pax;                         %returns array of plot axes
                clear obj
            end
        end
        
%%
        function ax = plotFileSpectrum(plttxt)
            %Plot a spectrum loaded from file
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            obj = ctWaveSpectra;
            if nargin<1
                xtxt = 'Wave period (s)';
                ytxt = 'Direction (degTN)';
                plttxt = {'Measured Spectral Energy (m^2s)',xtxt,ytxt};   
            end
            
            ok = 0; j = 0; ax = gobjects(0);
            while ok<1 
                obj = getSpectrumFile(obj);        %load spectrum data from file
                if isempty(obj),return; end

                obj = getMeasuredSpectrum(obj);    %compute spectrum based on measured form
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                %call plot function
                plttxt{4} = sprintf('File: %s',matlab.lang.makeValidName(obj.inpData.file));
                pax = getPlot(obj,plttxt,ptype);
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
                return;
            else
                tsprops = cobj.Data.Properties;
            end            

            obj.inpData.form = 'Measured';
            obj.inpData.source = 'Spectrum';
            xtxt = 'Wave period (s)';
            ytxt = 'Direction (degTN)';
            plttxt(1,:) = {'Spectral Energy (m^2s)',xtxt,ytxt};
            plttxt(2,:) = {'Spectral Energy (m^2s)',xtxt,ytxt};            

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
                plttxt{1,4} = sprintf('%s (%s)',tsdst.Description,dates{irow});

                obj(1) = getMeasuredSpectrum(obj(1));         %compute spectrum based on measured form
                obj(1).Params = wave_spectrum_params(obj(1)); %integral properties of spectrum
                skill = getSkillParameters(obj,mobj);
                %construct model spectrum
                obj(2) = setSpectrumModel(ctWaveSpectra);     %define the model to be used (Jonswap etc)
                obj(2) = getWaveModel(obj(2),obj(1).Params);
                if isempty(obj(2).Spectrum.SG), return; end

                %plot results
                plttxt{2,4} = getModelInputText(obj(2),1);   
                compareProperties(obj);
                ax(j+1,:) = plotObsModel(obj,plttxt,ptype);
                sax = omniSpectrumPlot(obj,plttxt);
                subtitle(sax,plttxt{2,4})
                
                %plot skill
                metatxt = {plttxt{1,4}, plttxt{2,4}};
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
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            anobj = ctWaveSpectra;     
            [obj,tsdst,plttxt,~] = getCaseInput(anobj,mobj);
            if isempty(tsdst), return; end

            tsdst.DataTable = rmmissing(tsdst.DataTable);   %remove nans
            tsdst = getsampleusingrange(tsdst);            

            if any(contains(tsdst.VariableNames,'Kurt'))
                obj = getMeasuredTS(anobj,tsdst);            %compute spectrum based on measured form
            else  
                dates = tsdst.DataTable.Properties.RowNames; 
                obj = getModelTS(obj,tsdst.DataTable,dates); %compute spectrum for specified conditions                   
            end

            xtxt = 'Wave period (s)';
            ytxt = 'Direction (degTN)';
            plttxt(1,:) = {'Spectral Energy (m^2s)',xtxt,ytxt};
            plttxt{1,4} = sprintf('%s',tsdst.Description);
            wrm_single_animation(mobj,obj,plttxt,ptype);
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
            hf = figure('Tag','PlotFig'); ax = axes(hf);
            spec = vertcat(obsobj(:).Spectrum);
            spgm = vertcat(modobj(:).spModel);
            plot(ax,[spec(:).date],[spgm(:).T_gamma],'x')

            plot_spectrum_model_skill(obsobj,modobj,mobj);
            parameterPlots(obsobj,modobj);
        end

%%
        function bimodalSpectrum(mobj)
            %analyse measured spectrum for bi-modality and explore
            %representing this in a model
            obj = ctWaveSpectra;            
            [obj,tsdst,plttxt] = getCaseInput(obj,mobj);
            if isempty(tsdst), return; end
            tsdst.DataTable = rmmissing(tsdst.DataTable);%remove nans
            dates = tsdst.DataTable.Properties.RowNames;

            ok = 0;
            while ok<1
                irow = listdlg("PromptString",'Select event to plot',...
                    'SelectionMode','single','ListSize',[160,300],...
                    'ListString',dates);
                if isempty(irow), ok = 1; continue; end
                plttxt{4} = sprintf('%s (%s)',tsdst.Description,dates{irow});
                obj(1).inpData.tsdst =  getDSTable(tsdst,irow,[]);  %selected record
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

                obj(2) = setDefaultSpectrumModel(ctWaveSpectra); %define the model to be used (Jonswap etc)
                obj(2).spModel.nspread = obj(2).spModel.nspread(1);
                if isempty(obj(2).spModel),return; end
                obj(2) = getWaveModel(obj(2),w_params);  %compute spectrum and params for specified conditions
                if isempty(obj(2).Spectrum.SG), return; end

                if isempty(idlof)
                    summary = [w_params;obj(2).Params];
                    summary.Properties.RowNames = {'Wind-waves','Model wind'};
                else
                    obj(3) = setDefaultSpectrumModel(ctWaveSpectra); %define the model to be used (Jonswap etc)
                    obj(3).spModel.nspread = obj(3).spModel.nspread(2);
                    if isempty(obj(3).spModel),return; end
                    obj(3) = getWaveModel(obj(3),s_params);  %compute spectrum and params for specified conditions
                    if isempty(obj(3).Spectrum.SG), return; end  
    
                    summary = [w_params;obj(2).Params;s_params;obj(3).Params];
                    summary.Properties.RowNames = {'Wind-waves','Model wind',...
                                               'Swell waves','Model swell'};
                    % obj(4) = copy(obj(2));
                    obj(2).Spectrum.SG = obj(2).Spectrum.SG+obj(3).Spectrum.SG; 
                    obj(3) = [];
                end
               
                display(summary)
                ax = omniSpectrumPlot(obj,plttxt);
                hold(ax,'on')
                plot(ax,[1,1]*w_params.Tp,ylim,'--b')
                if ~isempty(idlof)
                    plot(ax,[1,1]*s_params.Tp,ylim,'--r')
                    plot(ax,[1,1]*1/freq(idmn),ylim,'-.g')                    
                end
                hold(ax,'off')

                %surface plot comparison
                xtxt = 'Wave period (s)';
                ytxt = 'Direction (degTN)';
                txt(1,:) = {'Spectral Energy (m^2s)',xtxt,ytxt,''};
                txt(2,:) = {'Spectral Energy (m^2s)',xtxt,ytxt,''}; 
                txt{1,4} = sprintf('%s (%s)',tsdst.Description,dates{irow});
                txt{2,4} = getModelInputText(obj(2),1);   
                plotObsModel(obj,txt,'XY');
            end

        end
    end

%% ------------------------------------------------------------------------
% Input functions
%--------------------------------------------------------------------------
    methods
        function obj = setSpectrumModel(obj)
            %set spectrum form, data source, and parameters for wave
            %spectrum and directions spreading functions
            %  Defined using varargin as in above function   
            %  source iswind is used to prioritise selection of wind            
            sp = {'JONSWAP fetch limited','TMA shallow water','Bimodal sea (JONSWAP/TMA)'...
                   'Pierson-Moskowitz fully developed','Bretschneider open ocean'};
            src = {'Wave','Wind'};
            if nargin>1 && strcmp(obj.inpData.source,'Wind')
                src = {'Wind','Wave'};
            end
            
            spr = {'SPM cosine function','Donelan secant function'};
            % nsp = string(0:1:10);
            ptxt = sprintf('              Select values to use\nSpread exponent should be 0 if defined in data\nJonswap gamma only used in Jonswap and TMA models\nSet depth to apply saturation in TMA spectrum [0= no saturation]');
            selection = inputgui('FigureTitle','Spectrum',...
                                 'InputFields',{'Wave Spectra','Data type',...
                                    'Spread function','Spread exponent',...
                                    'Jonswap gamma','Water depth'},...
                                 'Style',{'popupmenu','popupmenu','popupmenu',...
                                         'edit','edit','edit'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{sp,src,spr,'2','3.3','0'},...%use nspread=0 if included in wave data
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
                    inp.Hs = str2double(inpt{1});
                    inp.Tp = str2double(inpt{2});
                    inp.Dir = str2double(inpt{3}); 
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
        function [obj,tsdst,plttxt,cobj] = getCaseInput(obj,mobj)
            %get selection and sort input based on data source
            [cobj,tsdst,dsname] = getCaseDataset(obj,mobj); %returns timeseries dst
            if isempty(cobj), tsdst = []; plttxt = {}; return; end
            
            xtxt = 'Wave period (s)';
            ytxt = 'Direction (degTN)';
            %plot labels: colorbar label, xlabel,ylabel,Tag,title
            plttxt = {'Modelled Spectral Energy (m^2s)',xtxt,ytxt};             

            obj.inpData.form = 'Measured'; 
            if strcmp(dsname,'Spectra')            %ts data of spectra
                obj.inpData.source = 'Spectrum';
                plttxt = {'Measured Spectral Energy (m^2s)',xtxt,ytxt}; 
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
                xtsdst = extract_wave_data(tsdst);
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
                inp.Hs = tsdst.Hs;
                inp.Tp = tsdst.Tp;
                inp.Dir = tsdst.Dir;
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
                fprintf('Input parameters: Hs=%.1f m, Tp=%.1f s, Dir=%.1f dTN\n',...
                                                    inp.Hs,inp.Tp,inp.Dir); 
            else
                fprintf('Input parameters: U=%.1f m/s, Dir=%.1f dTN, F=%.0f m\n',...
                                            inp.Uw,inp.Dir,inp.Fetch); 
            end
        end

%%
        function modeltxt = getModelInputText(obj,isone)
            %extract the model input to define a spectrum from spModel
            spm = obj.spModel;
            spmform = split(spm.form);
            if contains(spm.form,{'Pierson-Moskowitz fully developed','Bretschneider open ocean'})
                modeltxt = sprintf('%s and %s, n=%d ',...
                            spmform{1},spm.spread,spm.nspread);
            else
                if isnan(spm.gamma) && isone, spm.gamma = spm.T_gamma; end
                modeltxt = sprintf('%s, gamma=%.2g, and %s, n=%d ',...
                            spmform{1},spm.gamma,spm.spread,spm.nspread);
                if contains(spm.form,'TMA') && spm.depth>0
                    modeltxt = sprintf('%s, d=%.1f',modeltxt,spm.depth);
                end
            end      
        end


%%
        function dst = saveSpectrum(obj,mobj)
            %save an array of spectra to a dstable using the spt format
            %as used to import spectra from a wave buoy
            % input is obj.Spectrum with SG, dir and freq. 
            % saved in 64 frequency intervals
            flim = obj.Interp.flim;
            obsfreq = [flim(1):0.005:0.1,0.11:0.01:flim(2)];
            nvar = length(obj);
            varData = cell(1,nvar);
            for i=1:nvar
                stats = setSpectrum(obj(i),obsfreq);    
                varData{i} = table2cell(stats);
            end

            dsp = wave_cco_spectra('setDSproperties');
            %load the results into a dstable  
            dst.Spectra = dstable(varData{:},'RowNames',myDatetime,'DSproperties',dsp.dspec); 
            dst.Spectra.Dimensions.freq = obsfreq;
            
            %add properties to a dstable
            input = {obj.Params.Hs,obj.Params.T2,obj.Params.Sp,NaN};  %NaN is for SST
            dst.Properties =  dstable(input{:},'RowNames',myDatetime,'DSproperties',dsp.dsprop);

            dst.Spectra.Description = obj.inpData.tsdst.Description;
            dst.Properties.Description = obj.inpData.tsdst.Description;
            %save results
            cobj = ctWaveData;
            cobj.ModelType = 'Spectrum'; 
            setDataSetRecord(cobj,mobj.Cases,dst,'spectrum');
            getdialog(sprintf('Data loaded in class: %s',classname));
        end

%% ------------------------------------------------------------------------
% Spectrum construction functions
%--------------------------------------------------------------------------
        function obj = getModelSpectrum(obj)
            %construct model spectrum from input wind or wave conditions
            sp = obj.spModel;
            % 
            dir_int = obj.Interp.dir; %interval used to interpolate directions (deg)
            f_int = obj.Interp.freq;  %interval used to interpolate periods (s)
            f_lim = obj.Interp.flim;  %range of frequency bands (Hz)
        
            %frequencies at 1/per_int period intervals
            freq = f_lim(1):f_int:f_lim(2);       
            %direction at dir_int degree intervals
            dir = 0:dir_int:360-dir_int;
    
            %directional spreading factor for selected function  
            dir0 = obj.inpData.Dir;                   %mean wave direction 
            G = directional_spreading(dir0,dir,sp.nspread(1),sp.spread);
    
            %spectral energy for selected wave spectrum
            params = obj.inpData;
            if isnan(sp.gamma) && ~istable(obj.inpData)
                sp.gamma = 3.3;
            elseif isnan(sp.gamma)
                %use fit to Tz/Tm==T2/Tp based on Table 1 of MIAS No.4, 1982
                sp.gamma = 70*(obj.inpData.T2/obj.inpData.Tp)^12.23;
                if sp.gamma<1
                    sp.gamma = 0.9999;
                elseif sp.gamma>7.9
                    sp.gamma = 7.9999;
                end
            end
            obj.spModel.T_gamma = sp.gamma;
            S = wave_spectrum(sp.form,freq,params);
            if isempty(S), obj.Spectrum.SG = []; return; end
            obj.Spectrum.SG = G*S;
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir;

            if strcmp(sp.form,'TMA shallow water') && sp.depth>0
                g = 9.81;
                omega = 2*pi*freq.*sqrt(sp.depth/g);
                phi = 0.5*omega.^2.*(omega<=1) + 1.*(omega>2); 
                obj.Spectrum.SG = phi.*obj.Spectrum.SG;
            end
        end

%%
        function obj = getMeasuredSpectrum(obj)
            %construct measured spectrum from input spectrum data
            dir_int = obj.Interp.dir;  %interval used for directions (deg)
            f_int = obj.Interp.freq;   %interval used for frequecy (Hz)
            f_lim = obj.Interp.flim;   %range of frequency bands (Hz)
            
            %direction at dir_int degree intervals
            dir = 0:dir_int:360-dir_int;
            inp =  obj.inpData.tsdst;
            fspec = inp.Dimensions.freq';
            dirinp = {fspec,inp.Dir,inp.Spr,inp.Skew,inp.Kurt};
            G = datawell_directional_spectra(dir,false,dirinp{:}); %isplot=false       
            SG = (G.*inp.S(:)');      %force a row vector for S

            % interpolate spectrum to equal frequency intervals       
            freq = f_lim(1):f_int:f_lim(2);  

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
            %obj = setSprectrumInput(obj,0);         %add input needed for wave_spectrum.m  
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

%% ------------------------------------------------------------------------
% Plotting functions
%--------------------------------------------------------------------------
        function ax = getPlot(obj,plttxt,ptype,ax)
            %call plot function   
            if nargin<4, ax = []; end

            if strcmp(ptype,'XY')                %single x-y plot
                ax = surfPlot(obj,plttxt,ax);
            else
                ax = polarPlot(obj,plttxt,ax);
            end  
        end

%%
        function ax = surfPlot(obj,plttxt,ax) 
            %plot spectral surface as XY function of frequency and direction
            if nargin<3 || isempty(ax)
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end
            %extract data
            period = 1./obj.Spectrum.freq;
            dir = obj.Spectrum.dir; 
            dir0 = obj.Params.Dir;
            
            % dir2dir0 = dir-dir0;
            SG = obj.Spectrum.SG;
            SGpk = max(SG,[],'all');
            %make plot
            surf(ax,period,dir,SG,'Tag','PlotFigSurface');
            view(2);
            shading interp
            axis tight
            %ax.XLim(2) = obj.Params.Tp*2;
            hold on
            plot3(ax,xlim,dir0*[1 1],SGpk*[1,1],'--w','Tag','DirPk')
            plot3(ax,obj.Params.Tp*[1 1],ylim,SGpk*[1,1],'-.w','Tag','TpPk')
            hold off   
            xlabel(plttxt{2})
            ylabel(plttxt{3})

            %add the colorbar and labels   
            cb = colorbar;
            cb.Label.String = plttxt{1};
            cb(1).Tag = 'WaveSpectrum';            

            title(plttxt{4})
            p = obj.Params;
            subtitle(sprintf('Hs=%.2fm; Tp=%.1fs; T_2=%.1fs; Dir=%.3gdegTN',...
                                        p.Hs,p.Tp,p.T2,p.Dir));
        end

%%
        function ax = polarPlot(obj,plttxt,ax)
            %plot spectral surface as polar function of frequency and direction
            if nargin<3 || isempty(ax)
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end

            %extract data
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
            cb.Label.String = plttxt{1};

            text(max(X,[],'all')/2,max(Y,[],'all')*0.95,0,plttxt{2}); %radial label
            cb.Tag = 'WaveSpectrum';
            %improve positioning of colorbar
            cb.Position(3:4) = cb.Position(3:4)*0.6;

            title(plttxt{4})            
            p = obj.Params;
            subtitle(sprintf('Hs=%.2fm; Tp=%.1fs; T_2=%.1fs; Dir=%.3gdegTN',...
                                        p.Hs,p.Tp,p.T2,p.Dir));
        end

%%
        function ax = omniSpectrumPlot(obj,plttxt,ax)
            if nargin<3 || isempty(ax)
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end

            per = 1./obj(1).Spectrum.freq;
            hold(ax,'on')
            for i=1:length(obj)
                Sf = getOmniDirSpectrum(obj(i));
                plot(ax,per,Sf);
            end
            hold(ax,'off')
            xlabel(plttxt{1,2})
            ylabel(plttxt{1,1}) 
            title(ax,plttxt{1,4}) 
        end
%%
        function s = plotObsModel(obj,plttxt,ptype)
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
            getPlot(obj(1),plttxt(1,:),ptype,s(1));
            s(1).Title.String = 'Measured'; 
            s(1).CLim(2) = Smax;

            s(2) = subplot(m,n,2);
            getPlot(obj(2),plttxt(2,:),ptype,s(2));
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
                getPlot(diffobj,plttxt(1,:),ptype,s(3));  
                s(3).Title.String = 'Difference [Meas-Model]';
                s(3).Subtitle.String = [];
                s(3).XLim = s(1).XLim;
            end

            sgtitle(sprintf('Spectrum for: %s',plttxt{1,4}))
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
    end
%% ------------------------------------------------------------------------
% Utility functions
%--------------------------------------------------------------------------
    methods (Static, Access=private)


 %%
        % function [s1,s2] = obs_mod_plot(T,dir,var0,varc,ptype,params,plttxt)
        %     %plot offshore and inshore spectra  - similar to off_in_plot in SpecralTransfer
        % 
        %     hf = figure('Name','SpecTrans','Tag','PlotFig');
        %     ax = axes(hf);
        %     sel = params.modelsel;
        %     obs = params.obstable;
        %     mod = params.modtable;
        % 
        %     if strcmp(ptype,'XY')
        %         s1 = subplot(3,1,1,ax);
        %         ctWaveData.splot(T,dir,var0,plttxt,s1);  %method in SpectralTransfer
        %         s1.Title.String = 'Measured';        
        %         subtitle(s1,sprintf('Hs=%.2f m; Tp=%.1f s; T2=%.1f s; Dir=%.3g degTN',...
        %                                 obs.Hs,obs.Tp,obs.T2,obs.Dir),'Margin',1);
        % 
        %         s2 = subplot(3,1,2);
        %         plttxt{1} = 'Modelled Spectral Energy (m^2s)';
        %         ctWaveData.splot(T,dir,varc,plttxt,s2);
        %         s2.YLim = s1.YLim;
        %         s2.Title.String = sprintf('%s, gamma=%.2g, and %s, n=%d ',sel.form,...
        %                                         sel.gamma,sel.spread,sel.nspread);        
        %         subtitle(s2,sprintf('Hs=%.2f m; Tp=%.1f s; T2=%.1f s; Dir=%.3g degTN',...
        %                                 mod.Hs,mod.Tp,mod.T2,mod.Dir),'Margin',1);
        % 
        %         s3 = subplot(3,1,3);
        %         plttxt{1} = 'Difference in Spectral Energy (m^2s)';  
        %         spdiff = var0-varc;
        %         spmax = max(spdiff,[],'all');
        %         idx = spdiff<spmax/100 & spdiff>-spmax/100;
        %         spdiff(idx) = 0; %remove small differences
        %         ctWaveData.splot(T,dir,spdiff,plttxt,s3);
        %         s3.YLim = s1.YLim;   
        %         s3.Title.String = 'Difference [Meas-Model]';
        % 
        %         sgtitle(sprintf('Spectrum for: %s',plttxt{5}))
        %     else
        %         s1 = subplot(1,2,1,ax);
        %         ctWaveData.polar_plot(T,dir,var0,plttxt,s1);
        %         s1.Title.String = 'Measured';
        %         subtitle(s1,sprintf('Hs=%.2f m; Tp=%.1f s; T2=%.1f s; Dir=%.3g degTN',...
        %                                 obs.Hs,obs.Tp,obs.T2,obs.Dir),'Margin',1);
        % 
        %         s2 = subplot(1,2,2);
        %         plttxt{1} = 'Modelled Spectral Energy (m^2s)';
        %         ctWaveData.polar_plot(T,dir,varc,plttxt,s2);              
        %         s2.YLim = s1.YLim;
        %         s2.Title.String = sprintf('%s, gamma=%.2g, and %s, n=%d ',sel.form,...
        %                                         sel.gamma,sel.spread,sel.nspread);
        %         subtitle(s2,sprintf('Hs=%.2f m; Tp=%.1f s; T2=%.1f s; Dir=%.3g degTN',...
        %                                 mod.Hs,mod.Tp,mod.T2,mod.Dir),'Margin',1);
        % 
        %         sgtitle(sprintf('Spectrum for: %s',plttxt{5}))
        %     end
        % end  
        % function y = adapt_smooth_complex(x, rho, baseLen)
        %     % x: [1 x J] complex series, rho: [1 x J], baseLen: min window
        %     J = numel(x); y = zeros(1,J);
        %     for j = 1:J
        %         k = max(baseLen, round(baseLen + 10*(1 - rho(j))));
        %         h = max(1, floor(k/2));
        %         i0 = max(1, j-h); i1 = min(J, j+h);
        %         w = ones(1, i1-i0+1);
        %         y(j) = sum(x(i0:i1) .* w) / sum(w);
        %     end
        % end



    end
end