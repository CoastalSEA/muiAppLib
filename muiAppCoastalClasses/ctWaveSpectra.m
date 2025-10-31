classdef ctWaveSpectra < handle
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
        Spectrum = struct('SG',[],'freq',[],'dir',[],'depi',[])          
        Params %table for integral properties of the spectrum
        %struct with details of input definition depending on source
        inpData        % Wave - source,Hs,Tp,Dir,swl
                       % Wind  - source,AvSpeed,Dir,Fetch,zw,swl
                       % Spectrum - source,swl,issat,tsdst
        spModel  %struct with specification for spectrum model to use
        Interp   %struct for interpolation settings (defaults set in constructor)
        ModelMovie
    end

    methods
        function obj = ctWaveSpectra
            %class constructor
            obj.Interp.dir = 360/512; %interval used to interpolate directions (deg)
            obj.Interp.freq = 0.0005;  %interval used to interpolate frequency (Hz)
        end
    end
%%
    methods (Static)
        function ax = getPlotOption(mobj)
            %user selected plotting options
            listxt = {'Plot a spectrum using Case data',...
                      'Plot a spectrum using a model',...
                      'Plot a spectrum loaded from file',...
                      'Compare measured and modelled',...
                      'Timeseries animation',...
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

            end
        end

%%
        function ax = plotCaseSpectrum(mobj)
            %Plot a spectrum using case data
            obj = ctWaveSpectra;
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');
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
                ts = getDSTable(tsdst,irow,[]);    %selected record
                plttxt{4} = sprintf('%s (%s)',tsdst.Description,dates{irow});

                obj.inpData.tsdst = ts;            %assign data to input

                if strcmp(obj.inpData.source,'Spectrum')
                    obj = getMeasuredSpectrum(obj);%compute spectrum based on measured form
                else     
                    obj = setInputParams(obj,ts);  %add input needed to construct spectrum
                    obj = getModelSpectrum(obj);   %compute spectrum for specified conditions
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

                %call plot function   
                pax = getPlot(obj,plttxt,ptype);
                ax(1,j+1) = pax;                   %returns array of plot axes
            end
        end
        
%%
        function ax = plotModelSpectrum(plttxt)
            %Plot a spectrum using a model and user defined wave/wind parameters
            obj = ctWaveSpectra;
            if nargin<1
                xtxt = 'Wave period (s)';
                ytxt = 'Direction (degTN)';
                plttxt = {'Modelled Spectral Energy (m^2s)',xtxt,ytxt,'model'};   
            end
            inptype = questdlg('Wind or Wave input?',...
                                            'Input','Wind','Wave','Wave');
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            ok = 0; j = 0; ax = gobjects(0);
            while ok<1          
                obj = setSpectrumModel(obj);       %define the model to be used (Jonswap etc)
                if isempty(obj.spModel),return; end

                obj = setForcingConditions(obj,inptype); %UI to input wave/wind conditions
                if isempty(obj), return; end
                obj = setSprectrumInput(obj,0);    %add input needed for wave_spectrum.m, 0=data     
                obj = getModelSpectrum(obj);       %compute spectrum for specified conditions
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                inp = obj.inpData;
                fprintf('Input parameters: Hs=%.1f, Tp=%.1f, Dir=%.1f\n',...
                                                    inp.Hs,inp.Tp,inp.Dir);
                %call plot function    
                spm = obj.spModel;
                spmform = split(spm.form);
                plttxt{4} = sprintf('%s, gamma=%.2g, and %s, n=%d ',...
                                spmform{1},spm.gamma,spm.spread,spm.nspread);
                if contains(spm.form,'TMA') && spm.depth>0
                    plttxt{4} = sprintf('%s, d=%.1f',plttxt{4},spm.depth);
                end
                pax = getPlot(obj,plttxt,ptype);
                ax(1,j+1) = pax;                   %returns array of plot axes
            end
        end
        
%%
        function ax = plotFileSpectrum(plttxt)
            %Plot a spectrum loaded from file
            obj = ctWaveSpectra;
            if nargin<1
                xtxt = 'Wave period (s)';
                ytxt = 'Direction (degTN)';
                plttxt = {'Measured Spectral Energy (m^2s)',xtxt,ytxt};   
            end

            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');
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
            obj = ctWaveSpectra;
            classops = {'ctWaveData'};
            [~,tsdst,dsname] = getCaseDataset(obj,mobj,classops,1);
            if ~strcmp(dsname,'Spectra')
                warndlg('Spectral data required for this option')
                return;
            end
            casedesc = tsdst.Description;

            obj.inpData.form = 'Measured';
            obj.inpData.source = 'Spectrum';
            xtxt = 'Wave period (s)';
            ytxt = 'Direction (degTN)';
            plttxt(1,:) = {'Measured Spectral Energy (m^2s)',xtxt,ytxt};
            plttxt(2,:) = {'Modelled Spectral Energy (m^2s)',xtxt,ytxt};

            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');

            dates = tsdst.DataTable.Properties.RowNames;
            ok = 0; j = 0; ax = gobjects(0);
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                obj(1).inpData.tsdst = getDSTable(tsdst,irow,[]);    %selected record
                plttxt{1,4} = sprintf('%s (%s)',casedesc,dates{irow});

                obj(1) = getMeasuredSpectrum(obj(1));  %compute spectrum based on measured form
                obj(1).Params = wave_spectrum_params(obj(1)); %integral properties of spectrum
                
                obj(2) = ctWaveSpectra;
                obj(2) = setSpectrumModel(obj(2));     %define the model to be used (Jonswap etc)
                obj(2).inpData = obj(1).Params;
                obj(2).inpData.source = 'Wave';
                obj(2) = setSprectrumInput(obj(2),0);  %add input needed for wave_spectrum.m  
                obj(2) = getModelSpectrum(obj(2));    %compute spectrum for specified conditions
                obj(2).Params = wave_spectrum_params(obj(2)); %integral properties of spectrum
                spm = obj(2).spModel;    
                spmform = split(spm.form);
                plttxt{2,4} = sprintf('%s, gamma=%.2g, and %s, n=%d ',...
                                spmform{1},spm.gamma,spm.spread,spm.nspread);
                if contains(spm.form,'TMA') && spm.depth>0
                    plttxt{2,4} = sprintf('%s, d=%.1f',plttxt{2,4},spm.depth);
                end   
                
                ax(j+1,:) = plotObsModel(obj,plttxt,ptype);
            end
        end

        %%
        function animateCaseSpectrum(mobj)
            %
            obj = ctWaveSpectra;
            ptype = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');
            [obj,tsdst,plttxt,cobj] = getCaseInput(obj,mobj);
            if isempty(tsdst), return; end
            tsdst.DataTable = rmmissing(tsdst.DataTable);%remove nans

            while height(tsdst)>5000
                promptxt = sprintf('Times series contains %d records\nThis could take a while to run and genearte large file\nUse time sub-selection to extract shorter time period',...
                                                            height(tsdst));
                hw = warndlg(promptxt);
                waitfor(hw)
                tsdst = getsampleusingrange(tsdst);
            end

            %select from dates, check for NaNs, get spectrum and plot
            dates = tsdst.DataTable.Properties.RowNames;
            %ax = gobjects(0);
            hf = figure('Name','Animation', 'Units','normalized', ...
                    'Resize','on','HandleVisibility','on','Visible','on',...
                    'Position',[0.38,0.42,0.30,0.42],'Tag','PlotFig');

            
            %create an instance of muiPlots and populate the properties that are   
            %needed for the newAnimation method
            if isa(mobj.mUI.Plots,'muiPlots')
                pobj = mobj.mUI.Plots;    %get existing instance          
                clearPreviousPlotData(pobj);
            else
                pobj = muiPlots.get_muiPlots();                   %create new instance
            end       
            pobj.Plot.CurrentFig = hfig;
            pobj.Plot.FigNum = hfig.Number;
            pobj.ModelMovie = [];
            pobj.Title = sprintf('Case: %s',tsdst.Description);

            %extract the timeseries data and dimensions for plot
            Dims = obj.inpData.tsdst.Dimensions;
            pobj.Data.X = 1./Dims.freq;
            pobj.Data.Y = Dims.dir;
            pobj.Data.Z = {SGo,SGi};
            pobj.Data.T = tsdst.RowNames; 


            ax = axes(hf);
            ax.Position = [0.13,0.58,0.70,0.34]; %make space for slider bar
            nrec = length(dates);
            Mframes(nrec) = struct('cdata',[],'colormap',[]);
            for irow=1:nrec
                obj.inpData.tsdst = getDSTable(tsdst,irow,[]);    %selected record
                plttxt{4} = sprintf('%s (%s)',tsdst.Description,dates{irow});

                if strcmp(obj.inpData.source,'Spectrum')
                    obj = getMeasuredSpectrum(obj);%compute spectrum based on measured form  
                end
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                ax = getPlot(obj,plttxt,ptype,ax);
                Mframes(irow) = getframe(gcf);

                %ax(1,irow) = pax;                   %returns array of plot axes
            end


            obj.ModelMovie = Mframes;   %save movie to class property





        end
    end

%% ------------------------------------------------------------------------
% Input functions
%--------------------------------------------------------------------------
    methods
        function obj = setSpectrumModel(obj)
            %get spectrum form, data source, and parameters for wave
            %spectrum and directions spreading functions
            %  Defined using varargin as in above function   
            %  source iswind is used to prioritise selection of wind            
            sp = {'JONSWAP fetch limited','TMA shallow water',...
                   'Pierson-Moskowitz fully developed','Bretschneider open ocean'};
            src = {'Wave','Wind'};
            if nargin>1 && strcmp(obj.inpData.source,'Wind')
                src = {'Wind','Wave'};
            end
            
            spr = {'SPM cosine function','Donelan secant function'};
            nsp = string(0:1:10);
            ptxt = sprintf('              Select values to use\nSpread exponent should be 0 if defined in data\nJonswap gamma only used in Jonswap and TMA models\nSet depth to apply saturation in TMA spectrum [0= no saturation]');
            selection = inputgui('FigureTitle','Spectrum',...
                                 'InputFields',{'Wave Spectra','Data type',...
                                    'Spread function','Spread exponent',...
                                    'Jonswap gamma','Water depth'},...
                                 'Style',{'popupmenu','popupmenu','popupmenu',...
                                         'popupmenu','edit','edit'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{sp,src,spr,nsp,'3.3','0'},...%use nspread=0 if included in wave data
                                 'PromptText',ptxt,'Position',[0.25,0.56,0.25,0.1]);
            if isempty(selection)
                spectrum = [];
            else
                spectrum.form = sp{selection{1}};
                spectrum.source = src{selection{2}};
                spectrum.spread = spr{selection{3}};
                spectrum.nspread = str2double(nsp{selection{4}});
                spectrum.gamma = str2double(selection{5});
                spectrum.depth = str2double(selection{6});
                % satxt = 'excluding depth saturation';
                % if strcmp(spectrum.form,'TMA shallow water')
                %     satxt = 'including depth saturation';
                % end
                % spectrum.inptxt = sprintf('%s using %s (%s)\nSpreading function is %s with an exponent of %s. Jonswap gamma = %s',...
                %                  sp{selection{1}},src{selection{2}},satxt,...
                %                  spr{selection{3}},nsp{selection{4}},selection{5});
            end 
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
                    %inp.swl = str2double(inpt{4});  
                    inp.form = 'Model';

                case 'Wind'           %define spectrum using wind parameters                     
                    promptxt = {'Wind Speed (m/s)','Wind Direction (degTN)',...
                                'Height above msl (m)','Fetch Length (m)'};                                 
                    defaults = {'20.0','185','10.0','20000'};
                    inpt = inputdlg(promptxt,'Input conditions',1,defaults);
                    if isempty(inpt), return; end  %user cancelled
                    inp.AvSpeed = str2double(inpt{1});
                    inp.Dir = str2double(inpt{2});
                    inp.zW = str2double(inpt{3});
                    inp.Fetch = str2double(inpt{4});
                    %inp.swl = str2double(inpt{5});
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

            % promptxt = {'Still water level','Include depth saturation (1=Yes,0=No)'};
            % defaults = {'0.0','1'};
            % inpt = inputdlg(promptxt,'Input conditions',1,defaults);
            % if isempty(inpt), return; end  %user cancelled
            % inputs.swl = str2double(inpt{1});
            % inputs.issat = logical(str2double(inpt{2}));  
             
            [filename,path,~] = getfiles('MultiSelect','off',...
                'FileType',{'*.spt; *.txt'},'PromptText','Select file:');
            if filename==0, obj = []; return; end  %user cancelled
            varlist = {'',[path,filename]};
            specdst = wave_cco_spectra('getData',varlist{:});
            tsdst = horzcat(specdst.Spectra,specdst.Properties);
            tsdst = activatedynamicprops(tsdst,specdst.Properties.VariableNames);
            % inputs.swl = 0;
            % tsdst = addvars(tsdst,inputs.swl,'NewVariableNames','swl');
            inputs.form = 'Measured';
            inputs.source = 'Spectrum';
            inputs.file = filename;
            inputs.tsdst = tsdst;               
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
            %setup input parameters mapping variables from input dataset

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
            %
            if strcmp(obj.inpData.source,'Wave')
                obj.inpData.Hs = tsdst.Hs;
                obj.inpData.Tp = tsdst.Tp;
                obj.inpData.Dir = tsdst.Dir;
            elseif strcmp(obj.inpData.source,'Wind')
                obj.inpData.AvSpeed = tsdst.AvSpeed;
                obj.inpData.Dir = tsdst.Dir;
                obj.inpData.zW = tsdst.MetaData.zW;
                obj.inpData.Fetch = tsdst.MetaData.Fetch;
            end
            obj = setSprectrumInput(obj,1);  %1=timeseries input
        end

%%
        function obj = setSprectrumInput(obj,ists)
            %match input used in getForcingConditions when the input is a
            %measured spectrum and is called for same purpose in runWaves
            % inputs = []; 
            % promptxt = {'Include depth saturation (1=Yes,0=No)'};
            % defaults = {'1'};
            % inpt = inputdlg(promptxt,'Input conditions',1,defaults);
            % if isempty(inpt), obj = []; return; end  %user cancelled
            % inputs.issat = logical(str2double(inpt{1}));    
            % inputs.form = 'Measured';      
            % inputs.source = 'Spectrum'; 
            % inputs.tsdst = tsdst;
            % if inputs.logical(str2double(inpt{1}))
            %     inputs.txt = 'Measured spectrum including depth saturation';
            % else
            %     inputs.txt = 'Measured spectrum excluding depth saturation';
            % end
            % obj.inpData = inputs;
            % obj.spModel.inptxt = [];

            % if ists
            %     inp = obj.inpData.tsdst;
            % else
                
            % end
            inp = obj.inpData;
            sptype = obj.inpData.source;
            if strcmp(sptype,'Wave')
                obj.inpData.wvsp_inputs = {sptype,inp.Hs,inp.Tp};
            elseif strcmp(sptype,'Wind')
                obj.inpData.wvsp_inputs = {sptype,inp.AvSpeed,inp.zW,inp.Fetch};
            end        
        end
%%
        function inputMessage(obj)
            %write details of the input conditions to the command window
            inp = obj.inpData;
            if strcmp(inp.source,'Wave')
                fprintf('Input parameters: Hs=%.1fm, Tp=%.1fs, Dir=%.1fdTN\n',...
                                                    inp.Hs,inp.Tp,inp.Dir); 
            else
                fprintf('Input parameters: U=%.1fm/s, Dir=%.1fdTN, F=%.0fm\n',...
                                            inp.AvSpeed,inp.Dir,inp.Fetch); 
            end
        end
%% ------------------------------------------------------------------------
% Spectrum construction functions
%--------------------------------------------------------------------------
        function obj = getModelSpectrum(obj)
            %construct model spectrum from input wind or wave conditions
            sp = obj.spModel;
            % 
            dir_int = obj.Interp.dir;  %interval used to interpolate directions (deg)
            % radint = deg2rad(dir_int);
            f_int = obj.Interp.freq;  %interval used to interpolate periods (s)
            f_lower = 0.025;          %lower bound of frequency range  
            f_upper = 1;            %upper bound of frequency range  
        
            %frequencies at 1/per_int period intervals
            freq = f_lower:f_int:f_upper;         
            %offshore direction at dir_int degree intervals
            beta = 0:dir_int:360-dir_int;
    
            %directional spreading factor for selected function  
            dir0 = obj.inpData.Dir;                       %mean wave direction 
            diro = dir0-90:dir_int:dir0+90;
            G = directional_spreading(dir0,beta,sp.nspread,sp.spread);
            %G = interp1(diro,G,beta,'linear',0);
    
            %spectral energy for selected wave spectrum
            params = obj.inpData.wvsp_inputs;
            S = wave_spectrum(sp.form,freq,params{:},sp.gamma);
            obj.Spectrum.SG = G'*S;
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = beta;

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
            dir_int = obj.Interp.dir;  %interval used to interpolate directions (deg)
            % radint = deg2rad(dir_int);
            per_int = obj.Interp.per;  %interval used to interpolate periods (s)
            per_lower = 1;             %lower bound of period range  
            per_upper = 40;            %upper bound of period range  
        
            %frequencies at 1/per_int period intervals
            freq = 1./(per_lower:per_int:per_upper);             
            %direction at dir_int degree intervals for interpolation
            beta = 0:dir_int:360-dir_int;
            
            inp =  obj.inpData.tsdst;
            fspectra = inp.Dimensions.freq;
            dirinp = {fspectra,inp.Dir,inp.Spr,inp.Skew,inp.Kurt};
            %use range of directional means to define direction range which can
            %be greater than 180 if bidrectional sea-state
            % Dmnmx = minmax(inp.Dir,'omitnan');
            % dspectra = Dmnmx(1)-90:dir_int:Dmnmx(2)+90;
            G = datawell_directional_spectra(beta,false,dirinp{:}); %isplot=false

            obj.Spectrum.SG = (G.*inp.S);
            obj.Spectrum.freq = fspectra;
            obj.Spectrum.dir = beta; 
            obj.Params = wave_spectrum_params(obj); %integral properties of spectrum

            [XoFo,XiFi] = ctWaveSpectra.transferDims(beta,fspectra,beta,freq); 
            Gint = interpn(XoFo{:},G,XiFi{:},'linear',0);  
            Sint = interp1(fspectra,inp.S,freq,'linear',0);  

            %spectral energy for selected wave spectrum
            obj.Spectrum.SG = real(Gint.*Sint); 
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = beta;    
        end 

%%
        % function dst = setSpectrum(obj)
        %     %reduce a detailed spectrum to the format defined by the Datawell buoy
        %     %spt file format          
        %     SG = obj.Spectrum.SG;                 %spectrum to be saved
        %     dir = obj.Spectrum.dir;
        %     freq = obj.Spectrum.freq;
        % 
        %     obsfreq = obj.inpData.tsdst.Dimensions.freq;
        %     raddir = deg2rad(dir);
        %     %output: f, S(f)/Smax, Dir(f), Spr(f), Skew(f), Kurt(f) 
        %     %interpolate spectral density to each fequency
        %     S = trapz(raddir,SG,1);
        %     % figure; plot(freq,S);
        %     %NB spt file contains Sf/Smax but Sf is saved as S in wave_cco_spectra
        %     Sf = interp1(freq,S,obsfreq','pchip');
        %     Sf(Sf<0) = 0; 
        % 
        %     %extract the other output parameters at each frequency  
        %     m0 = trapz(raddir,abs(trapz(freq,SG,2)));        %zero moment
        %     nfreq = length(obsfreq);
        %     Dir = zeros(1,nfreq); Spr = Dir; Skew = Dir; Kurt = Dir;
        %     [F,D] = meshgrid(freq,dir);
        %     G = 1\S.*SG;
        %     for i=1:nfreq
        %         Sdf = interp2(F,D,G,obsfreq(i)*ones(size(dir)),dir);
        %         if all(isnan(Sdf)), continue; end
        %         % Normalize weights
        %         Sdf = Sdf/sum(Sdf);
        %         % figure; plot(dir,Sdf);
        % 
        %         % First trigonometric moments
        %         C = sum(Sdf.*cos(raddir));
        %         S = sum(Sdf.*sin(raddir));
        %         R = sqrt(C^2 + S^2);
        % 
        %         %mean direction
        %         meanDir_rad = atan2(S,C);
        %         Dir(1,i) = mod(rad2deg(meanDir_rad),360);
        % 
        %         % Directional spread (in degrees)
        %         spread_rad = sqrt(-2*log(R));
        %         Spr(1,i) = rad2deg(spread_rad);
        % 
        %         % Skewness and kurtosis (circular)
        %         % Skew(1,i) = (sum(Sdf.*sin(2*(raddir-meanDir_rad))));
        %         % Kurt(1,i) = (sum(Sdf.*cos(2*(raddir-meanDir_rad))));
        % 
        %         Skew(1,i) = skewness(Sdf);
        %         Kurt(1,i) = kurtosis(Sdf);
        %         %fprintf('%.3f, %.3e, %.1f, %.1f, %.1f, %.1f\n',obsfreq(i),Sf(i),Dir(i),Spr(i),Skew(i),Kurt(i));
        %     end
        % 
        %     myDatetime = obj.inpData.tsdst.RowNames;
        %     varData = {Sf,Dir,Spr,Skew,Kurt};
        % 
        %     opt = wave_cco_spectra('setDSproperties');
        %     %load the results into a dstable  
        %     dst.Spectra = dstable(varData{:},'RowNames',myDatetime,'DSproperties',opt.dspec); 
        %     dst.Spectra.Dimensions.freq = obsfreq;
        %     dst.Spectra.Description = obj.inpData.tsdst.Description;
        % 
        %     %add properties to table
        %     input = {obj.Params.Hs,obj.Params.T2,obj.Params.Spk,NaN};  %NaN is for SST
        %     dst.Properties =  dstable(input{:},'RowNames',myDatetime,'DSproperties',opt.dsprop);
        %     dst.Properties.Description = obj.inpData.tsdst.Description;
        % end

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
            SG = obj.Spectrum.SG;
            %make plot
            surf(ax,period,dir,SG,'Tag','PlotFigSurface');
            view(2);
            shading interp
            axis tight
            ax.XLim(2) = obj.Params.Tp*2;
            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = plttxt{1};
            xlabel(plttxt{2}); 
            ylabel(plttxt{3}); 
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
            radrange = [0,maxT];
            %interpolate var(phi,T) onto plot domain defined by tints,rints
            tints = linspace(0,2*pi,360);  
            rints = linspace(1,maxT,30);
            [Tq,Rq] = meshgrid(tints,rints); 
            warning('off',wid)
            vq = griddata(deg2rad(dir),period,SG',Tq,Rq);
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
        function s = plotObsModel(obj,plttxt,ptype)
            %plot offshore and inshore spectra  - similar to off_in_plot in SpecralTransfer
            hf = figure('Name','SpecTrans','Tag','PlotFig');
            ax = axes(hf);

            if strcmp(ptype,'XY')
                m = 3; n = 1;  %vertical plots of xy
            else
                m = 1; n = 2;  %horizontal plots if polar
            end
            casedesc = plttxt{1,4};

            s(1) = subplot(m,n,1,ax);
            getPlot(obj(1),plttxt(1,:),ptype,s(1));
            s(1).Title.String = 'Measured'; 

            s(2) = subplot(m,n,2);
            getPlot(obj(2),plttxt(2,:),ptype,s(2));

            if strcmp(ptype,'XY')
                s(3) = subplot(m,n,3);
                spdiff = obj(1).Spectrum.SG-obj(2).Spectrum.SG;
                spmax = max(spdiff,[],'all');
                idx = spdiff<spmax/100 & spdiff>-spmax/100;
                spdiff(idx) = 0; %remove small differences
                obj(1).Spectrum.SG = spdiff;
                getPlot(obj(1),plttxt(1,:),ptype,s(3));  
                s(3).Title.String = 'Difference [Meas-Model]';
                s(3).Subtitle.String = [];
            end

            sgtitle(sprintf('Spectrum for: %s',casedesc))

        end 
    end

%% ------------------------------------------------------------------------
% Utility functions
%--------------------------------------------------------------------------
    methods (Static, Access=private)
        function [PFW,XoFoWo,fw,fowo] = transferDims(InDir,fray,beta,freq)
            %replicate grid vectors to produce grids for interpn using
            %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
                [P,F] = ndgrid(InDir,fray);
                PFW = {P,F};  
                fw = {fray};
                %frequency, direction and water level arrays to interpolate to 
                %var(xso,fro)
                [Xso,Fro] = ndgrid(beta,freq);  
                XoFoWo = {Xso,Fro};  
                fowo = {freq};
        end
    end
end