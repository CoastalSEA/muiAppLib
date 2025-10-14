function ct_coastal_plots(mobj)                       
%
%-------function help------------------------------------------------------
% NAME
%   ct_coastal_plots.m
% PURPOSE
%   functions to do provide additional bespoke plot options using the 
%   wave and beach process data in coastal tools
% USAGE
%   ct_coastal_plots(mobj)
% INPUTS
%   mobj - ModelUI instance
% OUTPUT
%   plots for littoral drift (run directly from a wave time series),
%   scatter plot of 2 or 3 variables.
% NOTES
%    called as part of CoastalTools App.
% SEE ALSO
%   
%
% Author: Ian Townend
% CoastalSEA (c) May 2025
%--------------------------------------------------------------------------
%
    listxt = {'Littoral Drift','Generic Scatter Plot','Wave Scatter Plot',...
              'Frequency Analysis','Pos-Neg Change Plot'};
    % ok = 1;
    % while ok>0
        selection = listdlg("ListString",listxt,"PromptString",...
                            'Select option:','SelectionMode','single',...
                            'ListSize',[150,200],'Name','CoastalPlots');
        %if isempty(selection), ok = 0; continue; end
        if isempty(selection), return; end

        switch listxt{selection}
            case 'Littoral Drift'
                drift_plot(mobj);
            case 'Generic Scatter Plot'
                scatter_plot(mobj);    %muitoolbox function              
            case 'Wave Scatter Plot'
                wave_scatter_plot(mobj);
            case 'Frequency Analysis'
                frequency_analysis(mobj);
            case 'Pos-Neg Change Plot'
               positive_negative_plot(mobj);
        end
    % end
end

%%
function drift_plot(mobj)
    %get the inshore wave data, compute the littoral drift and plot as
    %monthly and annual values. Option to save results.
    promptxt = 'Select wave data set:'; 
    inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt); 
    if isempty(inwave), warndlg('No suitable data available'); return; end
    %retrieve an inshore wave data set
    [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
    if isempty(wv), return; end
    mtime = wv.RowNames;
    site = mobj.Inputs.ctWaveParameters; 

    %user selects which model to use
    tlist = {'CERC formula (SPM, 1994)',...
             'Dynamics of Marine Sands, Soulsby',...
             'Kamphuis formula',...
             'Damgaard & Soulsby (shingle)'};
    [h_dlg,ok] = listdlg('Name','Plot profile', ...
                         'PromptString','Select formula', ...
                         'ListSize',[200,80], ...
                         'SelectionMode','single', ...
                         'ListString',tlist);
    if ok==0, return; end 

    % input angle of shore, D50 and SPM drift coefficient
    [theta,d50,Kc] = getDriftSettings(site);

    %properties for bed slope within surf zone (half depth of inshore wave point)  
    ubs = site.UpperBeachSlope;
    z1km = site.BedLevelat1km;            
            
    %call drift model and add longhshore drift to wave time series
    g = mobj.Constants.Gravity;
    rhw = mobj.Constants.WaterDensity;
    rhs = mobj.Constants.SedimentDensity;
    vsc = mobj.Constants.KinematicViscosity;
    
    % dset = 'Dataset';
    % if isfield(waves.Data,'Properties'), dset = 'Properties'; end
    % wv = waves.Data.(dset);
    bs = profileslope(wv.depi/2,wv.swl,z1km,ubs); %first argument is depth
    Qall = littoraldrift(wv.Hsi,wv.Tp,wv.Diri,wv.depi,...
                                            theta,bs,d50,Kc,g,rhs,rhw,vsc);
    qs = Qall(:,h_dlg);
    dst = littoraldriftstats(qs,mtime,'month');

    if isa(dst,'struct')  %save results        
        dsp = driftDSproperties();    %add the full timeseries
        dst.Drift = dstable(qs,'RowNames',mtime,'DSproperties',dsp);
                       
        %assign metadata about model
        obj = CT_WaveModels;   %ModelType 1 for Littoral Drift
        dst.Drift.Source = sprintf('Class %s, using %s model','CT_WaveModels',...
                                                        obj.ModelName{1});
        formula = tlist{h_dlg};
        bs = mean(bs,'omitnan');
        zi = mean((wv.swl-wv.depi),'omitnan');
        mtxt1 = sprintf('Drift using %s; Theta=%g; d50=%g; Kc=%g; Beach slope=1:%.1f; Zi=%g',...
                                    formula,theta,d50,Kc,bs,zi);    
        mtxt2 = sprintf('Using %s case for wave input',wv.Description);
        output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);
        dst.Drift.MetaData = output.metatxt;
        %save results
        setDataSetRecord(obj,mobj.Cases,dst,'model');
        getdialog('Run complete');          
    else 
        %post result if single valued or text
        msgbox(dst,'Drift Results');
    end
end

%%
function [theta,d50,Kc] = getDriftSettings(site)
    %check that the current sit parameters settings are correct
    %modifications used to update RunParams stored with Case
    theta = site.ShorelineAngle;
    d50 = site.GrainSize;
    Kc = site.DriftCoefficient;
    promptxt = {'Shoreline angle (degTN)','Sediment grain size (m)',...
                'CERC Drift coefficient, Kc'};
    defaults = {num2str(theta),num2str(d50),num2str(Kc)};
    data = inputdlg(promptxt,'Drift settings',1,defaults);
    if isempty(data), return; end  %no change to default settings
    theta = str2double(data{1});
    d50 = str2double(data{2});
    Kc = str2double(data{3}); 
end

%%
function dsp = driftDSproperties() 
    %define a dsproperties struct and add the model metadata
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
    dsp.Variables = struct(...                       
        'Name',{'Qs'},...
        'Description',{'Alongshore drift rate potential'},...
        'Unit',{'m^3/s'},...
        'Label',{'Transport rate (m^3/s)'},...
        'QCflag',{'model'});
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

%%
function wave_scatter_plot(mobj)
    %plot wave scatter diagrams using depth or direction as 3rd variable
    promptxt = 'Select wave height:';    
    [H,indH] = get_variable(mobj,promptxt,'XYZmxvar',1);
    if isempty(H), return; end
    promptxt = 'Select wave period:';    
    [T,indT] = get_variable(mobj,promptxt,'XYZmxvar',1);
    if isempty(T), return; end  
    isvalid = checkdimensions(H.data,T.data);
    if ~isvalid, return; end
   
    answer = questdlg('Scatter plot using depth or direction?','Scatter plot',...
                      'Depth','Direction','Depth');
    if strcmp(answer,'Depth')
        promptxt = 'Select inshore depth:'; 
        [d,indd] = get_variable(mobj,promptxt,'XYZmxvar',1);
        if isempty(d), return; end  
        %input is a struct of dstables to be compatible with muiUserModel usage
        dst(1).data = getDSTable(indH.cobj.Data.Dataset,'VariableNames',H.name);
        dst(2).data = getDSTable(indT.cobj.Data.Dataset,'VariableNames',T.name);
        dst(3).data = getDSTable(indd.cobj.Data.Dataset,'VariableNames',d.name);
        wave_scatter(dst);
    else
        promptxt = 'Select wave direction:'; 
        [D,~] = get_variable(mobj,promptxt,'XYZmxvar',1);
        if isempty(D), return; end
        wave_scatter_3d(H.data,T.data,D.data);
    end
end

%%
function frequency_analysis(mobj)
    %time series analysis to generate a range of plots of frequency, 
    %probability of exceedance and duration of exceedance
    [cobj,~,dsname,ivar] = selectCaseDatasetVariable(mobj.Cases);
    if isempty(cobj), return; end
    dst = cobj.Data.(dsname);
    var = dst.DataTable{:,ivar};    
    t = dst.RowNames;
    vardesc = dst.VariableDescriptions{ivar};
    if contains(vardesc,'Water Level','IgnoreCase',true)
        res = waterlevelfreqplots(var,t);
    else
        res = frequencyanalysis(var,t,vardesc);
    end
    getdialog(res);
end

%%
function positive_negative_plot(mobj)
    %compute the rate of change of variable and plot histograms for positive
    %and negative components    
    [cobj,~,dsname,ivar] = selectCaseDatasetVariable(mobj.Cases);
    if isempty(cobj), return; end
    dst = cobj.Data.(dsname);
    var = dst.DataTable{:,ivar};
    t = dst.RowNames;
    
    vardesc = dst.VariableDescriptions{ivar};
    res = posneg_dv_stats(var,t,{vardesc,dst.Description});
    msgbox(res,'Pos-Neg plot');
end
