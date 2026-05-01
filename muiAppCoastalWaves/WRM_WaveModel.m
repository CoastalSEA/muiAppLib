classdef WRM_WaveModel < waveModels
%
%-------class help------------------------------------------------------
% NAME
%   WRM_WaveModel.m
% PURPOSE
%    Class for wave refraction using backward ray transfer function.
%    Constructs inshore time series from an offshore timeseries.Inherits 
%    abstract class waveModels which is a subclass of muiDataSet.
% SEE ALSO
%   muiDataSet, waveModels, WaveRayModel, RayTracks, SpectralTransfer
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
    end
    
    properties (Hidden)
        ModelType        %model used for the particular instance
    end                  %abstract property required by waveModels
    
    methods %(Access={?muiDataSet,?muiStats,?ctWaveData,?ctWaveSpectra,?subsample_dst})
        function obj = WRM_WaveModel()                    
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run the wave refraction spectral transfer model.
            obj = WRM_WaveModel;                            
            [dspec,dsprop] = setDSproperties(obj);
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            %get the timeseries input data and site parameters
            [tsdst,meta] = obj.getInputData(mobj);            
            if isempty(tsdst), return; end   %user cancelled data selection

            if strcmp(meta.inptype,'Spectrum')
                ModelType = 'transfer spectra';
            else
                ModelType = 'transfer timeseries';
            end

            %select a spectral transfer case to use
            [sptobj,sptmeta] = SpectralTransfer.getSTcase(mobj);
            if isempty(sptobj), return; end
            setRunParam(obj,mobj,meta.caserecs{:})     %assign run parameters
            %add spectral transfer selection to meta data
            inputxt = sprintf('%s, %s',meta.inptxt,sptmeta.inptxt);

            spobj = [];
            if ~strcmp(meta.inptype,'Spectrum')       
                spobj = setSpectrumModel(ctWaveSpectrum); 
            end

            %check whether wave spectra should also be saved
            nrec = height(tsdst(1));
            meta.issave = false;
            if nrec<=50000
                %only save results if arrays size is manageable
                answer = questdlg('Save the spectra?','Wave model','Yes','No','No');
                if strcmp(answer,'Yes') 
                    meta.issave = true; 
                    questxt = 'Save the full wave spectraor spt format';
                    savetype = questdlg(questxt,'Wave model','Full','SPT format','SPT format');
                end
            else
                getdialog(sprintf('Timeseries is too large to save spectra\nN = %d and current limit in WRM_WaveModel.runModel is 50000)',nrec))
            end
            
            tic
            [vartime,results,spectra] = runWaves(sptobj,tsdst,meta,spobj);
            elapsedTime = toc/60;  % Stop timer
            fprintf('Run time for %d steps: %.2f minutes\n',nrec,elapsedTime);
            if isempty(results), obj = []; return; end
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------  
            dst.Properties = dstable(results,'RowNames',vartime,'DSproperties',dsprop);                      
            %assign metadata about model            
            dst.Properties.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                         ModelType);
            dst.Properties.MetaData = inputxt;  
            % %add depths of inshore point for which there are backward rays
            % dst.Properties.UserData = depths;  
            
            if meta.issave
                specobj = ctWaveSpectrum;
                source = sprintf('Class %s, using %s',metaclass(obj).Name,ModelType); 
                if strcmp(savetype,'SPT format')
                    %saveSPTspectra outputs Spectrum and Properties but
                    %only the Spectrum dataset is saved
                    specobj.inpData.source = sprintf('Offshore spectrum from %s',tsdst.Description);
                    dst.OffshoreSpectra = saveSPTspectra(specobj,vartime,spectra,1);
                    dst.OffshoreSpectra.Source = source;
                    dst.OffshoreSpectra.MetaData = inputxt; 
                    specobj.inpData.source = sprintf('Inshore spectrum from %s',tsdst.Description);
                    dst.InshoreSpectra = saveSPTspectra(specobj,vartime,spectra,2);
                    dst.InshoreSpectra.Source = source;
                    dst.InshoreSpectra.MetaData = inputxt;                
                else
                    [dir,freq] = spectrumDimensions(specobj); 
                    nrec = numel(spectra);
                    hpw = PoolWaitbar(nrec,'Saving spectra');
                    parfor i=1:nrec
                        Sot(i,:,:) = spectra(i).Sot;
                        Sit(i,:,:) = spectra(i).Sit;
                        increment(hpw)
                    end
                    dst.oiSpectra = dstable(Sot,Sit,'RowNames',vartime,...
                                                     'DSproperties',dspec); 
                    dst.oiSpectra.Dimensions.dir = dir;    %NB order is X,Y and must
                    dst.oiSpectra.Dimensions.freq = freq;  %match variable dimensions  
                    delete(hpw)
                    %assign metadata about model
                    dst.oiSpectra.Source = source;
                    dst.oiSpectra.MetaData = inputxt;  
                end
            end
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------  
            %save results
            obj.ModelType = 'Inwave_model';           %added to match ctWaveModel
            setDataSetRecord(obj,mobj.Cases,dst,'model');
            getdialog('Run complete');
        end

%% --------------------------------------------------------------------------
% Other utilities
%--------------------------------------------------------------------------
        function runBatchMode(mobj)
            %run spectral transfer for multiple points and load a set of
            %dstables as datasets in a single case
            obj = WRM_WaveModel;                            
            [~,dsprop] = setDSproperties(obj);
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            %get the timeseries input data and site parameters
            muicat = mobj.Cases;
            [tsdst,meta] = obj.getInputData(mobj);            
            if isempty(tsdst), return; end   %user cancelled data selection

            %select a spectral transfer case to use
            [sptrecs,sptmeta] = SpectralTransfer.getSTcase(mobj,true);
            caserecs = [meta.caserecs,sptmeta.caserecs];
            setRunParam(obj,mobj,caserecs{:})     %assign run parameters

             %define the model to be used (Jonswap etc)
            inptype = meta.inptype; 
            if strcmp(inptype,'Spectrum') 
                spmodelobj = [];
            else
                spmodelobj = setSpectrumModel(ctWaveSpectrum); 
            end

            meta.issave = false;  %do not save spectra for a batch run
            ModelType = 'transfer timeseries';

            hw = waitbar(0,'Processing point 0');
            npnts = length(sptrecs);
            tic
            for i=1:npnts      %NOT parfor because used in runWaves               
                waitbar(i/npnts,hw,sprintf('Processing point %d',i));
                sptobj = getCase(muicat,sptrecs(i));
                [vartime,results,~] = runWaves(sptobj,tsdst,meta,spmodelobj);

                adst = dstable(results,'RowNames',vartime,'DSproperties',dsprop);
                %assign metadata about model            
                adst.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                             ModelType);
                adst.MetaData = sprintf('%s, %s used for spectral transfer',meta.inptxt,...
                                              sptobj.Data.Inshore.Description);
                %add depths of inshore point for which there are backward rays
                adst.UserData = sptobj.Data.Inshore.UserData.Depths;
                pname{i} = sprintf('Point%d',i);
                pdst(i) = adst; clear results
            end
            elapsedTime = toc/60;  % Stop timer
            fprintf('Run time for %d points: %.2f minutes\n',npnts,elapsedTime);
            
            for j = 1:npnts
                dst.(pname{j}) = pdst(j);
            end
            delete(hw)
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------  
            %save results
            obj.ModelType = 'Inwave_model';           %added to match ctWaveModel
            setDataSetRecord(obj,mobj.Cases,dst,'model');
            getdialog('Run complete');          
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab            
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 

%%
    methods (Access=private)
        function results = callRefraction(~,tsdst,site_params)
            %parse inputs for call to refraction function
            z0 = site_params.OffshoreBedLevel;
            theta0 = site_params.OffshoreAngle;
            Kf = site_params.FrictionCoefficient;
            
            prompt = {'Angle of bed contour at refraction point(s) (degTN)',...
                      'Bed level at refraction point(s) (mOD)'};
            default = {num2str(theta0),num2str(-100)};
            answer = inputdlg(prompt,'Deepwater refraction',1,default);
            if isempty(answer), return; end            
            thetaT = str2double(answer{1});
            zT = str2double(answer{2});
            dep0 = tsdst.swl-z0;   %source water depth from swl to bed level 
            depT = tsdst.swl-zT;   %target water depth from swl to bed level 
            
            %this is a general refraction model so 'isshore' is false
            [Hso,Diro] = refraction(tsdst.Hs,tsdst.Tp,tsdst.Dir,[dep0,depT],...
                                                [theta0,thetaT],Kf,false); 
            results = {Hso,tsdst.Tp,Diro};
        end
%%
        function [dspec,dsprop] = setDSproperties(~)
            %define the variables in the dataset
            %define the metadata properties for the demo data set
            dspec = struct('Variables',[],'Row',[],'Dimensions',[]); dsprop = dspec;    
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dspec.Variables = struct(...
                'Name',{'So','Si'},...
                'Description',{'Offshore spectral density','Inshore spectral density'},...
                'Unit',{'m^2/Hz','m^2/Hz'},...
                'Label',{'Spectral density (m^2/Hz)','Spectral density (m^2/Hz)'},...
                'QCflag',repmat({'raw'},1,2)); 
            dspec.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dspec.Dimensions = struct(...    
                'Name',{'dir','freq'},...
                'Description',{'Direction','Frequency'},...
                'Unit',{'degTN','Hz'},...
                'Label',{'Direction (degTN)','Frequency (Hz)'},...
                'Format',{'',''});    

            dsprop.Variables = struct(...   
                'Name',{'Hs','m0','Dir','Sp','Tp','Dp',...
                        'Sfdpk','Tfdpk','Dfdpk','T1','T2','T10'...
                        'kw','kt2','ktp','kd','swl','depi'},...
                'Description',{'Inshore wave height',...
                               'Inshore zero moment',...
                               'Inshore wave direction',...
                               'Inshore peak spectral density'...
                               'Inshore peak period',...
                               'Inshore peak direction',...
                               'Inshore f-d spectral density peak',...
                               'Inshore f-d peak period',...
                               'Inshore f-d peak direction',...
                               'Inshore mean period (T1)',...
                               'Inshore zero-upcross period (T2)',...
                               'Inshore energy period (T-10)',...
                               'Wave transfer coefficient',...
                               'Mean period coefficient',...
                               'Peak period coefficient',...
                               'Mean direction shift',...
                               'Still water level',...
                               'Inshore depth'},...                               
                'Unit',{'m','m^2','degTN','m^2/Hz','s','degTN','m^2/Hz','s',...
                            'degTN','s','s','s','-','-','-','deg','mOD','m'},...
                'Label',{'Wave height (m)','Zero moment (m^2)',...
                         'Wave direction (degTN)','Spectral density (m^2/Hz)',...                
                         'Wave period (s)','Wave direction (degTN)',... 
                         'Spectral density (m^2/Hz)','Wave period (s)',... 
                         'Wave direction (degTN)','Wave period (s)',...
                         'Wave period (s)','Wave period (s)',...
                         'Transfer coefficient, kw','Transfer coefficient, kt2',...
                         'Transfer coefficient, ktp','Direction shift (deg)',...
                         'Water level (mOD)','Water depth (m)'},...
                'QCflag',repmat({'model'},1,18)); 
            dsprop.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dsprop.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''}); 
        end   
    end  
%%
    methods (Static, Access=private)
        function off = addWaveConditions(SGo,Dims,off)
            %when using wind input define the offshore wave conditions
            g = 9.81;
            off.Tp = 0.54*g^-0.77*off.AvSpeed.^0.54.*off.Fetch.^0.23;

            dir_int = 0.5;     %interval used to interpolate directions (deg)
            radint = deg2rad(dir_int);
            So = trapz(radint,abs(trapz(Dims.freq,SGo,2))); %integral of offshore spectrum
            off.Hs = 4*sqrt(So);
        end
    end
end    