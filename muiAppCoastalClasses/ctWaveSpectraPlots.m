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
        % Interp         % rlim - radial limits for polar plot
        %                % tlim - theta limits for polar plot
        %MetaData       % additional variable selection details 
        hFig = gobjects(0)         % figure handle
        %ModelMovie
    end

    methods
        function obj = ctWaveSpectraPlots
            %class constructor
            %Interp property
            % obj.Interp.rlim = {1,20,20};      %radial limits for polar plot
            % obj.Interp.tlim = {0,2*pi,360};   %theta limits for polar plot
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
                    obj = ctWaveSpectraPlots.plotCaseSpectrum(mobj);
                case 2                  %Plot a spectrum using a model
                    obj = ctWaveSpectraPlots.plotModelSpectrum();
                case 3                 %Plot a spectrum loaded from file
                    obj = ctWaveSpectraPlots.plotFileSpectrum();
                case 4                  %Compare case spectra
                    obj = ctWaveSpectraPlots.cfCaseSpectra(mobj);
                case 5                  %Compare measured and modelled
                    obj = ctWaveSpectraPlots.cfSpectrum2Model(mobj);
                case 6                  %Timeseries animation
                    ctWaveSpectraPlots.animateCaseSpectrum(mobj);
                case 7                  %Fit a model to a timeseries of measured spectr
                    ctWaveSpectraPlots.fitModel2Measured(mobj);
                case 8                  %decompose measured spectrum to bimodal form
                    ctWaveSpectraPlots.bimodalSpectrum(mobj);
            end
        end

%%
        function runPlotOption(src,~,mobj)
            %main menu callback as alternative to getPlotOption (used in
            %WaveRayModel)
            switch src.Text
                case 'Case'                  %Plot a spectrum using case data
                    ctWaveSpectraPlots.plotCaseSpectrum(mobj);
                case 'Model'                 %Plot a spectrum using a model
                    ctWaveSpectraPlots.plotModelSpectrum();
                case 'SPT File'              %Plot a spectrum loaded from file
                    ctWaveSpectraPlots.plotFileSpectrum();
                case 'cf Cases'              %Compare case spectra
                    ctWaveSpectraPlots.cfCaseSpectra(mobj);
                case 'SPT v Model'           %Compare measured and modelled
                    ctWaveSpectraPlots.cfSpectrum2Model(mobj);
                case 'Animation'             %Timeseries animation
                    ctWaveSpectraPlots.animateCaseSpectrum(mobj);
                case 'ModelvMeasured skill'  %Fit a model to a timeseries of measured spectr
                    ctWaveSpectraPlots.fitModel2Measured(mobj);
                case 'Bimodal analysis'      %decompose measured spectrum to bimodal form
                    ctWaveSpectraPlots.bimodalSpectrum(mobj);
            end
        end

%% ------------------------------------------------------------------------
% Static plotting methods
%--------------------------------------------------------------------------
        function obj = plotCaseSpectrum(mobj,vis,~)
            %Plot a spectrum using case data
            if nargin<2, vis = 'on'; end
            ptype = ctWaveSpectraPlots.plotType();

            %get Case dataset to be used
            [~,tsdst,~] = ctWaveSpectraPlots.getCaseInputParams(mobj);
            if isempty(tsdst), return; end
        
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
                obj.spModel.selection = spselect;
                obj = getSpectrumObject(obj,tsdst,irow);
                if isempty(obj.spModel)
                    return
                else
                    spselect = obj.spModel.selection;
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
            [cobj,tsdst,~] = ctWaveSpectraPlots.getCaseInputParams(mobj,1);  %1 limits to ctWaveData class
            if isempty(tsdst), return; end

            %select record from dateset get spectrum and plot
            dates = tsdst.DataTable.Properties.RowNames;
            ok = 0; isfirst = true;
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                %construct measured spectrum
                obj = ctWaveSpectraPlots;             %initialise class object
                obj = getSpectrumObject(obj,tsdst,irow);
                %add spectrum input properties
                obj.inpData.properties = cobj.Data.sptProperties.DataTable(irow,:);
                 %construct model spectrum
                obj(2) = setSpectrumModel(ctWaveSpectraPlots);     %define the model to be used (Jonswap etc)
                obj(2) = getWaveModel(obj(2),obj(1).Params);
                if isempty(obj(2).Spectrum.SG), return; end

                %plot results
                obj(2) = getModelInputText(obj(2));   
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
            [~,tsdst,~] = ctWaveSpectraPlots.getCaseInputParams(mobj);
            if isempty(tsdst), return; end

            tsdst(1) = getsampleusingrange(tsdst(1));  
            if isempty(tsdst(1)), return; end         %invalid selection

            anobj = ctWaveSpectraPlots();
            if any(contains(tsdst(1).VariableNames,'Kurt'))
                obj = getMeasuredTS(anobj,tsdst);    %compute spectrum based on measured form
            else  
                obj = getModelTS(anobj,tsdst); %compute spectrum for specified conditions                   
            end

            obj(1).Plotxt.ttxt = sprintf('%s',tsdst.Description);
            wrm_single_animation(obj,mobj,ptype);
        end

%%
        function fitModel2Measured(mobj)
            %create a measured and model timeseries of spectra to determine
            %the best model fit parameters for a given sea location
            [~,tsdst,~] = ctWaveSpectraPlots.getCaseInputParams(mobj,1);  %1 limits to ctWaveData class
            if isempty(tsdst), return; end

            tsdst = getsampleusingrange(tsdst);  
            if isempty(tsdst), return; end       %invalid selection

            anobj = ctWaveSpectraPlots;
            %get timeseries of measured spectra
            obsobj = getMeasuredTS(anobj,tsdst); 
 
            spectra = [obsobj(:).Spectrum]; 
            dates = [spectra(:).date];
            intable = vertcat(obsobj(:).Params);
            intable = dstable(intable,'RowNames',dates);
            intable.MetaData.sptype = 'Wave';
            modobj = getModelTS(anobj,intable);
            if isempty(modobj), return; end   %user cancelled
            %check on gamma values
            if isnan(anobj.spModel.gamma)
                hf = figure('Tag','PlotFig'); ax = axes(hf);
                spec = vertcat(obsobj(:).Spectrum);
                spgm = vertcat(modobj(:).spModel);
                plot(ax,[spec(:).date],[spgm(:).T_gamma],'x')
            end

            %plot model skill and allow user to examine individual parameters
            plotSpectrumModelSkill(obsobj,modobj,mobj)
            parameterPlots(obsobj,modobj);
        end
%%
        function bimodalSpectrum(mobj)
            %analyse measured spectrum for bi-modality and explore
            %representing this in a model
            %get the measured wave spectrum to be modelled
            [~,tsdst,~] = ctWaveSpectraPlots.getCaseInputParams(mobj,1);  %1 limits to ctWaveData class
            if isempty(tsdst), return; end

            %select record from dateset get spectrum and plot
            dates = tsdst.DataTable.Properties.RowNames;
            obj = ctWaveSpectraPlots;             %initialise class object
            ok = 0;
            while ok<1
                irow = listdlg("PromptString",'Select event to plot',...
                    'SelectionMode','single','ListSize',[160,300],...
                    'ListString',dates);
                if isempty(irow), ok = 1; continue; end
                obj.Plotxt.ttxt = sprintf('%s (%s)',tsdst.Description,dates{irow});
                tsdstrow = ctWaveSpectrum.getDatasetRow(tsdst,irow);
                %set input parameters for selected record
                obj = setInputParams(obj,tsdstrow);
                obj = getMeasuredSpectrum(obj);  %compute spectrum based on measured form
                if isempty(obj.Spectrum.SG); return; end
                freq = obj.Spectrum.freq;
                dir = obj.Spectrum.dir;
                SG = obj.Spectrum.SG;
                Sf = getOmniDirSpectrum(obj);
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
                anobj = setSpectrumModel(ctWaveSpectraPlots); 
                obj(2) = copy(anobj);
                obj(2).spModel.nspread = anobj.spModel.nspread(1);
                obj(2).spModel.gamma = anobj.spModel.gamma(1);
                if isempty(obj(2).spModel),return; end
                obj(2) = getWaveModel(obj(2),w_params);  %compute spectrum and params for specified conditions
                if isempty(obj(2).Spectrum.SG), return; end

                if isempty(idlof)
                    summary = [w_params;obj(2).Params];
                    summary.Properties.RowNames = {'Wind-waves','Model W-W'};
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
                    summary.Properties.RowNames = {'Wind-waves','Model W-W',...
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
                obj(2) = getModelInputText(obj(2));  
                %obj(2).Plotxt.vtxt = 'Modelled Spectral Energy (m^2s)';  
                hf = plotObsModel(obj,'XY');
                %add button to access wave parameters
                addDataButton(obj,hf,summary); 

                obj = ctWaveSpectraPlots; %reset to blank instance
            end
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
    methods (Access=private)
        function hf = getMultiPlot(obj)
            %plot a set of axes in a single figure
            hf = figure('Name','SpecTrans','Tag','PlotFig');
            nplot = numel(obj);
            t = tiledlayout(hf, nplot, 1); % nplot rows, 1 column
            hfigs = [obj(:).hFig];
            for i=1:nplot
                % Move existing axes into tiles
                sax = findobj(hfigs(i).Children,'Type','axes');
                tt = nexttile(t, i);
                delete(tt)           % removes the placeholder axes
                sax.Parent = t;   % Reparent to tiledlayout in figure
                sax.Layout.Tile = i;
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
            %extract Skill parameters using MS_RunParams class for input
            x = obj(1).Spectrum.dir;              %direction
            y = obj(1).Spectrum.freq;             %frequency
            robj = MS_RunParams.setInput(mobj);   %default or current values if user cancels

            skill.Ro = robj.maxcorr;
            skill.n  = robj.skillexponent;
            skill.Inc = true;                   %flag to include skill score
            skill.W = robj.skillwindow;
            if isempty(skill.W), skill.Inc = false; end
            subdomain = robj.skillsubdomain;
            skill.SD = ctWaveSpectraPlots.getSubDomain(x,y,subdomain);
            skill.iter = robj.skilliteration;
        end

%%
        function plotSpectrumModelSkill(obsobj,modobj,mobj)
            %compute the skill of model v measured spectrum data and produce Taylor
            %plot of timeseries results
            %calls MS_RunParams class and taylor_plot function    
            skill = getSkillParameters(obsobj(1),mobj);    %get the parameters for skill model
        
            %get the statistics for the Taylor plot
            stats = get_skill_stats(obsobj,modobj,skill);
        
            %add the timeseries results to the Taylor plot
            ndteststd = [stats(:).teststd]./[stats(:).refstd]; %normalised std
            rLim = ceil(max(ndteststd));                       %radial limit for the plot
            ax = taylor_plot_figure(rLim);    
            metatxt = {'Measured','Model'};
            ax = taylor_plot_ts(ax,stats,skill,metatxt);
            % subtitletxt = getModelInputText(modobj(1));
            ax.Title.String = obsobj(1).inpData.source;
            ax.Subtitle.String = modobj(1).Plotxt.stxt;    
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
    end

%% ------------------------------------------------------------------------
% Static plotting utilities
%--------------------------------------------------------------------------  
    methods (Static, Access=private)
        function ptype = plotType()
            %define 
            ptype = questdlg('What type of plot','XY Polar','XY','Polar','XY'); 
        end

%%
        function [cobj,tsdst,dsname] = getCaseInputParams(mobj,~)
            %get the case dataset input parameters
            if nargin==2
                [cobj,tsdst,dsname] = ctWaveSpectraPlots.getCaseDataset(mobj,{'ctWaveData'},1);
            else
                [cobj,tsdst,dsname] = ctWaveSpectraPlots.getCaseDataset(mobj);
            end
            if isempty(tsdst), return; end
            tsdst = ctWaveSpectraPlots.getInputParams(cobj,tsdst,dsname);  %extract required variables
            if isempty(tsdst), return; end
            tsdst(1).DataTable = rmmissing(tsdst(1).DataTable);%remove nans
        end

%%
        function [cobj,tsdst,dsname] = getCaseDataset(mobj,classops,idd)
            %get selection and load case. option to limit classopt in call
            if nargin<2
                idd = [];
                classops = {'ctWaveData','ctWindData','WRM_WaveModel','muiUserModel'};
            elseif nargin<3
                idd = [];                
            end
            promptxt = 'Select Case to use:';
            [cobj,~,datasets,idd] = selectCaseDataset(mobj.Cases,...
                                          [],classops,promptxt,idd);
            if isempty(cobj), tsdst = []; dsname = []; return; end
            dsname = datasets{idd};
            tsdst = cobj.Data.(dsname);
        end

%%        
        function xtsdst = getInputParams(cobj,tsdst,dsname)
            %check for valid variable names when timeseries wave or wind
            %data are used to define conditions
            xtsdst = [];
            if strcmp(dsname,'sptProperties')     %ts data of spectra properties
                warndlg('Select ''sptSpectra'' dataset when using spt input')
                return;
            elseif strcmp(dsname,'sptSpectra')         %ts data of spectra 
                intype = 'Spectrum'; meta = [];
                xtsdst = tsdst;
            else
                if isa(cobj,'ctWaveData')
                    intype = 'Wave';
                    [xtsdst,meta] = extract_wave_data(tsdst); %returns 1xN array if multi-modal
                elseif isa(cobj,'WRM_WaveModel') 
                    intype = 'Wave';
                    [xtsdst,meta] = extract_wave_data(tsdst); %returns 1xN array if multi-modal
                elseif isa(cobj,'ctWindData')
                    intype = 'Wind';
                    [xtsdst,meta] = extract_wind_data(tsdst,1); %isfetch=true
                else
                    warndlg('muiUserModel not yet handled in ctWaveSpectraPlots.getInputParams')
                    return
                end 

            end
            if isempty(xtsdst), return; end

            newmeta = struct('dstmeta',xtsdst(1).MetaData,'sptype',intype,...
                                                     'sptmeta',meta);
            xtsdst(1).MetaData = newmeta;            
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
    end
end