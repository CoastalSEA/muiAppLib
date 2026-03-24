classdef ctWaveSpectrum < matlab.mixin.Copyable
%
%-------class help------------------------------------------------------
% NAME
%   ctWaveSpectrum.m
% PURPOSE
%   Class that creates and holds a wave spectrum in terms of spectral density 
%   as a function of direction and frequency
% NOTES
%   Spectral data are imported to the ctWaveData class using functions
%   such as wave_cco_spectrum. Spectrum records hold two datasets named as
%   sptSpectrum and sptProperties. Wave data can be simple unimodal or 
%   multimodal descriptions of the sea state.
% SEE ALSO
%   see ct_costal_plots and SpectralTransfer for examples of use. used by
%   ctWaveSpectraPlots class as a super class
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2025
%--------------------------------------------------------------------------
%    
    properties 
        %struct for Spectral density matrix, S(dir,f) and dimensions (freq and dir)
        Spectrum = struct('SG',[],'freq',[],'dir',[],'date',[],'depth',[])          
        Params %table for integral properties of the spectrum
        %struct with details of input definition depending on source
        inpData        % Wave - Hs,Tp,Dir,ds
                       % Wind  - Uw,zW,Dir,Fetch,df,ds
                       % Spectrum - S,Dir,Spr,Skew,Kurt,f,ds,issat
                       % depth at site - ds, includes any variation in swl
                       % gamma - wave record specific gamma value
                       % - all include fields for input, output and source
        %struct with specification for spectrum model to use               
        spModel        % form - wave spectrum model
                       % source - wind or wave data
                       % spread - relationship to define directional spreading
                       % nspread - exponent for directional spreading
                       % gamma - model peakiness exponent in JONSWAP spectrum
                       % depth - site depth for saturation limit
                       % inptxt - summary text for spectrum model
                       % selection - last model selections used
        %struct for interpolation settings (defaults set in constructor)
        Interp         % dir - interval used to interpolate directions (deg)
                       % freq - interval used to interpolate frequency (Hz)
                       % flim - frequency limits
                       % rlim - radial limits for polar plot
                       % tlim - theta limits for polar plot
        %struct for plotting text of wave spectra
        Plotxt         % xtxt - label for x-axis
                       % ytxt - label for y-axis
                       % vtxt - label for variable in surface plot
                       % ttxt - label for plot title   
                       % stxt - label for plot subtitle
    end

    methods
        function obj = ctWaveSpectrum
            %class constructor
            obj.Interp.dir = 360/512;       %interval used to interpolate directions (deg)
            obj.Interp.freq = 0.005;        %interval used to interpolate frequency (Hz)
            obj.Interp.flim = [0.025,0.58]; %frequency limits
            obj.Interp.rlim = {1,20,20};    %radial limits for polar plot
            obj.Interp.tlim = {0,2*pi,360}; %theta limits for polar plot            
            %default plot labels
            obj.Plotxt.xtxt = 'Wave period (s)';
            obj.Plotxt.ytxt = 'Direction (degTN)';
            obj.Plotxt.vtxt = 'Spectral Energy (m^2s)';  
            obj.Plotxt.ttxt = '';
            obj.Plotxt.stxt = '';
        end
    end
%% ------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
    methods (Static)
        function tsdstrow = getDatasetRow(tsdst,irow)
            %extract single row
            dates = tsdst(1).DataTable.Properties.RowNames;
            w_dst = getDSTable(tsdst(1),irow,[]);     %selected record
            tsdstrow(1) = w_dst;                      %assign data to input
            for j=2:numel(tsdst)
                %tsdst(1) defines wind-wave. get same record for swell(s)
                idx = find(tsdst(j).RowNames==dates{irow});
                s_dst = getDSTable(tsdst(j),idx,[]);  %selected record
                tsdstrow(j) = s_dst; %#ok<AGROW>
            end
        end
    end

    methods
%% ------------------------------------------------------------------------
% set functions
%--------------------------------------------------------------------------
function obj = setSpectrumModel(obj)
            %set spectrum form, data source, and parameters for wave
            %spectrum and directions spreading functions           
            sp = {'JONSWAP fetch limited','TMA shallow water',...
                   'Pierson-Moskowitz fully developed','Bretschneider open ocean'};

            if isempty(obj.spModel) || ...
               isfield(obj.spModel,'selection') && isempty(obj.spModel.selection)
                def1 = {1,1,'','',''};
                def2 = {'2','0','0'};
            else
                def1 = [obj.spModel.selection(1:2),{'','',''}];
                def2 = obj.spModel.selection(3:end);
            end

            spr = {'SPM cosine function','Donelan secant function'};
            % nsp = string(0:1:10);
            ptxt = inputHelpTxt();
            selection = inputgui('FigureTitle','Spectrum',...
                                 'InputFields',{'Wave Spectra',...
                                    'Spread function','Spread exponent',...
                                    'Jonswap gamma','Water depth'},...
                                 'Style',{'popupmenu','popupmenu',...
                                         'edit','edit','edit'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',[{sp},{spr},def2(:)'],...%use nspread=0 if included in wave data
                                 'UserData',def1,...  %initial selection for popupmenus
                                 'PromptText',ptxt,'Position',[0.25,0.56,0.25,0.1]);
            if isempty(selection)
                spectrum = [];
            else
                spectrum.form = sp{selection{1}};
                spectrum.spread = spr{selection{2}};
                spectrum.nspread = str2num(selection{3}); %#ok<ST2NM> handles vectors
                spectrum.gamma = str2num(selection{4});   %#ok<ST2NM>
                spectrum.depth = str2double(selection{5});
                satxt = 'excluding depth saturation';
                if strcmp(spectrum.form,'TMA shallow water') && spectrum.depth>0
                    satxt = 'including depth saturation';
                end
                spectrum.inptxt = sprintf('%s (%s)\nSpreading function is %s with an exponent of %s. Jonswap gamma = %s',...
                                 sp{selection{1}},satxt,...
                                 spr{selection{2}},selection{3},selection{4});
                spectrum.selection = selection;
            end 
            obj.spModel = spectrum;

            %nested function-----------------------------------------------
            function txt = inputHelpTxt()
                %txt1 = '              Select values to use';
                txt2 = 'Selection: Spread exponent should be 0 if defined in data';
                txt3 = 'Jonswap gamma only used in Jonswap and TMA models';
                txt4 = 'Setting gamma>0 overrides built-in relationships';
                txt5 = 'If input defines sea and swell, gamma and spread can accept multiple values';
                txt6 = 'Set depth to apply saturation in TMA spectrum [0 = no saturation]';

                txt = sprintf('%s\n%s\n%s\n%s\n%s',txt2,txt3,txt4,txt5,txt6);
            end
        end

%%
        function obj = setForcingConditions(obj)
            %set the user input for wave or wind conditions. If intype is
            %known call setWaveForcing or setWindForcing directly
            intype = questdlg('Wind or Wave input?','Input',...
                                              'Wind','Wave','Quit','Wave');
            switch intype
                case 'Quit'
                    obj.inpData = []; return;
                case 'Wave'           %define spectrum using wave parameters
                    obj = setWaveForcing(obj);
                case 'Wind'           %define spectrum using wind parameters                     
                    obj = setWindForcing(obj);
                otherwise
                    getdialog('Unknown type of source in ctWaveSpectrum.getForcingConditions')
                    obj.inpData = []; return; 
            end
        end

%%
        function obj = setWindForcing(obj)
            %set the user input for wind conditions
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
            inp.input = 'Wind';
            inp.output = 'Modelled';
            inp.source = 'User';
            inp.date = 'now';
            obj.inpData = inp;
        end

%%
        function obj = setWaveForcing(obj)
            %set the user input for wave conditions
            ok = 0;
            inpt = {'1.1','8.2','185'};
            while ok<1
                promptxt = {'Wave height (m)','Peak period (s)',...
                                 'Wave direction (degTN)'};
                inpt = inputdlg(promptxt,'Input conditions',1,inpt);
                if isempty(inpt), return; end  %user cancelled
                inp.Hs = str2num(inpt{1}); %#ok<ST2NM>
                inp.Tp = str2num(inpt{2}); %#ok<ST2NM>
                inp.Dir = str2num(inpt{3}); %#ok<ST2NM>
                %check inputs are the same length
                if numel(inp.Hs)==numel(inp.Tp) && numel(inp.Hs)==numel(inp.Dir) 
                    ok = 1;
                end
            end
            inp.input = 'Wave';
            inp.output = 'Modelled';
            inp.source = 'User';
            inp.date = 'now';
            obj.inpData = inp;
        end

%%
        function obj = setSpectrumFile(obj)
            %load spectrum from a file
            [filename,path,~] = getfiles('MultiSelect','off',...
                'FileType',{'*.spt; *.txt'},'PromptText','Select file:');
            if filename==0, obj = []; return; end  %user cancelled
            varlist = {'',[path,filename]};
            specdst = wave_cco_spectra('getData',varlist{:});
            inputs.spectrum = specdst.sptSpectrum;
            inputs.properties = specdst.sptProperties;            
            inputs.input = 'Spectrum';
            inputs.output = 'Measured';
            inputs.source = specdst.sptSpectrum.Description;
            inputs.date = specdst.sptSpectrum.RowNames;
            obj.inpData = inputs;
            obj = setPlotText(obj);               
        end

%%          
function obj = setInputParams(obj,tsdst,inptype)
            %extract input data from dstable and assign in format needed for
            %spectrum function
            if strcmp(inptype,'Spectrum')
                inp.spectrum = tsdst; 
                inp.input = 'Spectrum';
                inp.output = 'Measured';
                inp.date = inp.spectrum.RowNames;   

            elseif strcmp(inptype,'Wave')

                for i=1:numel(tsdst)
                    inp.Hs(i) = tsdst(i).Hs;
                    inp.Tp(i) = tsdst(i).Tp;
                    inp.Dir(i) = tsdst(i).Dir;
                end
                inp.gamma = obj.spModel.gamma;
                inp.input = 'Wave';
                inp.output = 'Modelled';
                inp.date = tsdst(1).RowNames; 

            elseif strcmp(inptype,'Wind')
                inp.Uw = tsdst.AvSpeed;
                inp.Dir = tsdst.Dir;
                inp.zW = tsdst.MetaData.zW;
                inp.Fetch = tsdst.MetaData.Fetch;
                inp.gamma = obj.spModel.gamma;
                inp.input = 'Wind';
                inp.output = 'Modelled';
                inp.date = tsdst.RowNames; 

            end
            inp.source = tsdst.Description; 
            inp.ds = [];
            
            obj.inpData = inp;           
        end

%%
        function gamma = setGamma(obj,seastate)
            %replicate gamma if single valued and set values if <=0
            ncomp = numel(obj.inpData.Hs);   %number of sea state components

            if obj.spModel.gamma(1)<=0
                %estimate the value of gamma for the total sea state
                istma = contains(obj.spModel.form,'TMA');
                [~,f] = spectrumDimensions(obj);                
                gamma0 = obj.spModel.gamma(1);
                gamma = estimate_jonswap_gamma(seastate,f,gamma0,istma);
                if ncomp>1
                    gamma = repmat(gamma,1,ncomp);
                end                
            else
                gamma = obj.spModel.gamma;
            end
        end
%% ------------------------------------------------------------------------
% get functions
%--------------------------------------------------------------------------

        function obj = getSpectrumObject(obj,meta,tsdst,irow)
            %get the spectrum and wave properties for selected input
            tsdstrow = ctWaveSpectrum.getDatasetRow(tsdst,irow);
     
            if ~strcmp(meta.inptype,'Spectrum')        
                obj = setSpectrumModel(obj);            %define model if not measured spectrum
                if isempty(obj.spModel), return; end    %user cancelled
            end

            %set input parameters for selected record
            obj = setInputParams(obj,tsdstrow,meta.inptype);
            %get the spectrum data for selected record
            if isempty(meta.variables) || isempty(meta.variables.seastate)
                seastate = [];
            else
                seastate = meta.variables.seastate.DataTable(irow,:);
            end
            obj = getSpectrum(obj,seastate);  
            if isempty(obj.Spectrum.SG), return; end
            inputMessage(obj);

            %option to return diagnostics struct to investigate
            %multi-modality
            % [params,diagn] = wave_spectrum_params(obj,true); %integral properties of spectrum
            % obj.Params = params;
        end

%%
        function obj = getSpectrum(obj,meta)
           %return measured or model spectrum depending on data inputs in
           %calling instance of ctWaveSpectrum.        
            if strcmp(obj.inpData.input,'Spectrum') 
                obj = getMeasuredSpectrum(obj);         %compute spectrum based on measured form
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                obj.Params.Properties.RowNames = {'Properties'};

            elseif strcmp(obj.inpData.input,'Wind') 
                obj = getModelSpectrum(obj);            %compute spectrum for specified conditions
                obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
                obj.Params.Properties.RowNames = {'Properties'};

            else
                obj = getMultiModalSpectrum(obj,meta);    %compute spectrum for specified conditions                  
            end
            obj = setPlotText(obj);   
        end

%%
        function obj = getMeasuredSpectrum(obj)
            %construct measured spectrum from input spectrum data
            [dir,freq] = spectrumDimensions(obj);
            % dir_int = obj.Interp.dir;  %interval used for directions (deg)
            % f_int = obj.Interp.freq;   %interval used for frequecy (Hz)
            % f_lim = obj.Interp.flim;   %range of frequency bands (Hz)

            inp =  obj.inpData.spectrum(1); %waveModels.getInputParams loads sptSpectrum as first tsdst
            if ~(any(strcmp(inp.VariableNames,'S')))
                warndlg('Invaldid input data type'); return; 
            end
            fspec = inp.Dimensions.freq';
            dirinp = {fspec,inp.Dir,inp.Spr,inp.Skew,inp.Kurt};
            G = datawell_directional_spectra(dir,false,dirinp{:}); %isplot=false       
            SG = (G.*inp.S(:)');      %force a row vector for S

            SGint = zeros(size(SG,1),length(freq));
            for i = 1:size(SG,1)
                SGint(i,:) = interp1(fspec, SG(i,:), freq, 'pchip', 0);
            end 

            obj.Spectrum.SG = SGint;
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir; 
            obj.Spectrum.date = inp.RowNames;
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
            [S,gamma] = wave_spectrum(sp.form,freq,params);
            if isempty(S), obj.Spectrum.SG = []; return; end

            obj.inpData.gamma = gamma;  %update gamma if modified in wave_spectrum
            obj.Spectrum.SG = G*S;
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir;
            %handle depth saturation if TMA spectrum used            
            if strcmp(sp.form,'TMA shallow water') && sp.depth>0
                obj.Spectrum.depth = sp.depth;
            end
        end

%%
        function obj = getMultiModalSpectrum(obj,seastate)
            %construct spectrum from input wave conditions 
            ncomp = numel(obj.inpData.Hs);   %number of sea state components
            %pad input variables if required - assumes either scalar or correct number           
            if numel(obj.spModel.nspread)~=ncomp
                obj.spModel.nspread = repmat(obj.spModel.nspread,1,ncomp);
            end
            
            if ~isempty(seastate)
                gamma = setGamma(obj,seastate);
            else
                %if negative input gamma and no seastate set as +ve value
                if obj.spModel.gamma(1)<0, fact = -1; else, fact = 1; end
                gamma = repmat(obj.spModel.gamma(1),1,ncomp)*fact;
            end

            %loop over each component to get component spectrum
            wavecomp(1,ncomp) = copy(obj);
            idx = 1;
            for i=1:ncomp 
                wavecomp(i) = copy(obj);
                wavecomp(i).spModel.nspread = obj.spModel.nspread(i);
                % wavecomp(i).spModel.gamma = gamma(i);
                compparams.Hs = obj.inpData.Hs(i);
                compparams.Tp = obj.inpData.Tp(i);
                compparams.Dir = obj.inpData.Dir(i);
                compparams.gamma = gamma(i);

                if compparams.Hs==0 || compparams.Tp==0 
                    %create zero spectrum if no waves (can be wind or swell)
                    wavecomp(i) = zeroSpectrum(wavecomp(i),compparams);
                else
                    wavecomp(i) =  getWaveModel(wavecomp(i),compparams);  %compute spectrum and params for specified conditions
                end 

                if ~isempty(wavecomp(i).Params), idx = [idx,i+1]; end %#ok<AGROW>                
            end
            spect = [wavecomp(:).Spectrum];
            obj.Spectrum.freq = wavecomp(1).Spectrum.freq;
            obj.Spectrum.dir = wavecomp(1).Spectrum.dir;
            nf = numel(obj.Spectrum.freq);  
            nd = numel(obj.Spectrum.dir); 
            obj.Spectrum.SG = sum(reshape([spect.SG],nd,nf,[]),3);            
            inp = [wavecomp(:).inpData];
            obj.inpData.gamma = [inp(:).gamma];% inp(1).gamma;   % for Jonswap swell [inp(:).gamma];
               
            rownames = {'Combined','Wind-waves','Primary swell waves','Secondary swell waves'};

            combined = wave_spectrum_params(obj);
            
            waveparams = vertcat(wavecomp.Params);
            if height(waveparams)>1
                summary = [combined;waveparams];                
            else
                summary = waveparams; idx = 2; %only wind-waves in results
            end

            summary.Properties.RowNames = rownames(idx);
            obj.Params = summary;
        end

%%
        function meas_obj = getMeasuredTS(obj,tsdst)
            %use a timeseries of wave spectra input data to create spectra timeseries
            nrec = height(tsdst);
            meas_obj(nrec,1) = obj;      
            hpw = PoolWaitbar(nrec, 'Processing measured timeseries');  %and increment(hpw);
            parfor i=1:nrec                                   %parfor loop
                anobj = copy(obj);                
                itsdst = getDSTable(tsdst,i,[]);              %selected record
                anobj = setInputParams(anobj,itsdst,'Spectrum');
                anobj = getMeasuredSpectrum(anobj);           %compute spectrum based on measured form
                anobj.Params = wave_spectrum_params(anobj);   %integral properties of spectrum
                meas_obj(i,1) = anobj;
                increment(hpw);
            end
            delete(hpw)
        end

%%
function mod_obj = getModelTS(obj,tsdst,meta)
            %use a timeseries of wave data to create a timeseries of spectra  
            obj = setSpectrumModel(obj);
            if isempty(obj.spModel), mod_obj = []; return; end
  
            dates = tsdst(1).DataTable.Properties.RowNames; 
            nrec = length(dates);
            mod_obj(nrec,1) = obj;
            hpw = PoolWaitbar(nrec, 'Processing measured timeseries');  %and increment(hpw);
            parfor i=1:nrec                                 %parfor loop
                anobj = copy(obj);
                itsdst = ctWaveSpectrum.getDatasetRow(tsdst,i); %selected record
                anobj = setInputParams(anobj,itsdst,meta.inptype);
                if strcmp(meta.inptype,'Wind')
                    anobj = getModelSpectrum(anobj);        %compute spectrum for specified conditions
                else
                    if isempty(meta.variables) || isempty(meta.variables.seastate)
                        seastate = [];
                    else
                        seastate = meta.variables.seastate.DataTable(i,:);
                    end
                    anobj = getMultiModalSpectrum(anobj,seastate);%compute spectrum for specified conditions
                end
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
            %only used in getMultiModalSpectrum
            obj.inpData = params;
            obj.inpData.ds = 0;
            obj.inpData.input = 'Wave';
            obj.inpData.output = 'Modelled';            
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
            % Calculate the Kitaigorodskii limit to the spectrum at frequency f
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

%%
        function [idpks,idmn] = spectrumPeaks(obj)
            %find the peaks in a unimodal spectrum and the minima between
            %the peaks
                Sf = getOmniDirSpectrum(obj);
                [maxSf,~] = max(Sf);        %maximum spectral density

                minpeakdist = 10;        %minimum no of points separating peaks **
                minpeakht = maxSf*0.2;   %minimum height of peaks **
                [mxlocs,~] = peakseek(Sf,minpeakdist,minpeakht); 
                [mnlocs,~] = peakseek(-Sf,minpeakdist,-minpeakht); 

                idmns = zeros(1,numel(mxlocs)-1);
                for i=1:numel(mxlocs)-1
                    idx = find(mnlocs>mxlocs(i)  & mnlocs<mxlocs(i+1)); 
                    if ~isempty(idx)
                        [~,idl] = min(Sf(mnlocs(idx)));
                        idmns(i) = mnlocs(idx(idl));
                    end
                end

                if sum(idmns)==0
                     idpks = mxlocs(1); idmn = [];
                else
                    idmns = idmns(idmns>0);
                    [~,idd] = min(Sf(idmns));
                    idmn = idmns(idd);
                    if isempty(idmn)
                        idpks = mxlocs(1);  %only a single separated peak
                    else
                        idpks(1) = mxlocs(find(mxlocs<idmn,1));
                        idpks(2) = mxlocs(find(mxlocs>idmn,1));
                    end
                end
        end
%% ------------------------------------------------------------------------
% Save spectrum
%--------------------------------------------------------------------------
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
                stats = setspectrum(obj(i),obsfreq);    
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
% Plotting metadata
%--------------------------------------------------------------------------
%%
        function obj = setPlotText(obj)
            %assign additional metadata and plotting text
            obj.Plotxt.ttxt = sprintf('%s: %s',obj.inpData.source,...
                                                 string(obj.inpData.date));
            if strcmp(obj.inpData.input,'Spectrum')     %ts data of spectra
                obj.Plotxt.vtxt = 'Measured Spectral Energy (m^2s)';
                obj.Plotxt.ttxt = sprintf('%s: %s',obj.inpData.source,...
                                                 string(obj.inpData.date));
            else
                obj =  setModelInputText(obj); 
                obj.Plotxt.vtxt = 'Modelled Spectral Energy (m^2s)';
            end
        end

%%
        function obj = setModelInputText(obj)
            %extract the model input to define a spectrum from spModel property
            spm = obj.spModel;
            if isempty(spm), return; end   %no Spectrum model defined
            spmform = split(spm.form);
            if contains(spm.form,{'Pierson-Moskowitz fully developed','Bretschneider open ocean'})
                ttxt = sprintf('%s and %s, spread=%d ',...
                                        spmform{1},spm.spread,spm.nspread);
            else
                if spm.gamma(1)<=0 && ~isempty(obj.inpData)
                    spm.gamma = obj.inpData.gamma;
                end
                ttxt = sprintf('%s, gamma=%.2g, and %s, spread=%d ',...
                        spmform{1},spm.gamma(1),spm.spread,spm.nspread(1));
            end
        
            %handle TMA depth limit
            if contains(spm.form,'TMA') && spm.depth>0
                ttxt = sprintf('%s, d=%.1f',ttxt,spm.depth);
            end
        
            stxt = '';            
            if strcmp(obj.inpData.input,'Wave')
                ninp = numel(obj.inpData.Hs);
                for i=2:ninp
                % for i=1
                    if i>2, stxt = sprintf('%s; ',stxt); end
                    stxt = sprintf('%sswell-%d: gamma=%.2g; spread=%d',...
                                    stxt,i-1,spm.gamma(i),spm.nspread(i));
                    % stxt = sprintf('%sswell-%d: spread=%d',...
                    %                 stxt,i-1,spm.nspread(i));
                end
            end
        
            if isempty(stxt)
                obj.Plotxt.stxt = ttxt;
            else
                obj.Plotxt.stxt = sprintf('%s\n%s',ttxt,stxt);
            end
        end

%% ------------------------------------------------------------------------
% Utility functions
%--------------------------------------------------------------------------
        function plotobj = Spectrum2SpectralPlots(obj)
            %create ctWaveSpectralPlots instance from ctWaveSpectrum
            nrec = numel(obj);
            plotobj(nrec) = ctWaveSpectraPlots;
            for i=1:nrec
                plotobj(i) = ctWaveSpectraPlots;
                plotobj(i).Spectrum = obj(i).Spectrum;     
                plotobj(i).Params = obj(i).Params;
                plotobj(i).inpData = obj(i).inpData;     
                plotobj(i).spModel = obj(i).spModel;
                plotobj(i).Interp = obj(i).Interp;
                plotobj(i).Plotxt = obj(i).Plotxt; 
            end
        end

%%
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
        function obj = zeroSpectrum(obj,params)
            %spectrum object with zero energy density (eg for Hs=0)
            [dir,freq] = spectrumDimensions(obj);
            obj.Spectrum.freq = freq;
            obj.Spectrum.dir = dir;
            obj.Spectrum.SG = zeros(numel(dir),numel(freq));
            obj.inpData = params;
            obj.inpData.ds = [];
            obj.inpData.input = 'Wave';
            obj.inpData.output = 'Modelled';
            obj.inpData.gamma = 0.0;
            
            obj.Params = wave_spectrum_params(obj); %integral properties of spectrum
        end

%%
        function [spectra,params] = unpackSpectrum(inobj,offobj)
            %unpack the spectrum property as a set of arrays
            % 1 - input includes inobj and offobj: return spectra and params
            % 2 - input just inobj: returns params
            nrec = length(inobj);
            hpw = PoolWaitbar(nrec, 'Unpacking spectra');
            if nargin==2
                parfor i=1:nrec                               %parfor loop
                    time(i,1) = offobj(i).Spectrum.date;    
                    swl(i,1) = offobj(i).inpData.swl;
                    Sot(i,:,:) = offobj(i).Spectrum.SG;
                    Sit(i,:,:) = inobj(i).Spectrum.SG;
                    depths(i,1) = inobj(i).Spectrum.depth;  
                    params(i,:) = inobj(i).Params;
                    increment(hpw);
                end   
                spectra = struct('time',time,'swl',swl,'Sot',Sot,'Sit',Sit,'depths',depths);
            else
                parfor i=1:nrec                               %parfor loop
                    params(i,:) = inobj(i).Params;
                    increment(hpw);
                end
                spectra = [];
            end
            delete(hpw)
        end

%%
        function inputMessage(obj)
            %write details of the input conditions to the command window
            inp = obj.inpData;
            if strcmp(inp.input,'Wave')
                if ~isfield(inp,'Hs') && isfield(inp,'tsdst')
                    inp = inp.tsdst;
                end

                for i=1:length(inp.Hs)
                    fprintf('Input parameters: Hs=%.1f m, Tp=%.1f s, Dir=%.1f dTN\n',...
                                           inp.Hs(i),inp.Tp(i),inp.Dir(i)); 
                end
            elseif strcmp(inp.input,'Wind')
                fprintf('Input parameters: U=%.1f m/s, Dir=%.1f dTN, F=%.0f m\n',...
                                            inp.Uw,inp.Dir,inp.Fetch); 
            else
                fprintf('Wave spectrum for %s\n',string(inp.date));
            end
        end
    end
end