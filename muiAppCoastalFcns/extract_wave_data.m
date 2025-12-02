function wvdst = extract_wave_data(inwvdst)
%
%-------function help------------------------------------------------------
% NAME
%   extract_wave_data.m
% PURPOSE
%   Extract Hs, Tp and Dir from a dataset that does not use default naming
%   convention (e.g. Copernicus re-analysis data)
% USAGE
%   wvdst = extract_wave_data(inwvdst)
% INPUTS
%   inwvdst - dstable of data to use for data selection
% OUTPUT
%   wvdst - dstable with default naming of Hs, Tp and Dir
% SEE ALSO
%   used in ctWaveModel and WRM_WaveModel
%   
% Author: Ian Townend
% CoastalSEA (c) April 2025
%--------------------------------------------------------------------------
%  
    varnames = inwvdst.VariableNames;
    if sum(ismatch(varnames,{'Hs','Tp','Dir'}))==3 
        wvdst = copy(inwvdst);
        return; 
    end 

    %handle unimodal and bimodal input data options
    Ntype = 1; t_txt = {'Wave data'};
    answer = questdlg('Type of data','Wave selection','unimodal','bimodal','unimodal');
    if strcmp(answer,'bimodal')
        Ntype = 2; t_txt = {'Wind-wave','Swell wave'}; 
    end
    vardesc = inwvdst.VariableDescriptions;

    for i=1:Ntype
        %call UI to select all required fields
        sel = getInputUI(vardesc,t_txt{i});
        if isempty(sel), wvdst = []; return; end

        factor = 1;
        if contains(vardesc{sel{2}},'mean')
            factor = 1.2; %scale mean period to peak period. this value            
        end               %is for Jonswap with gamma=3.3
    
        inwv = inwvdst.DataTable;
        indata = {inwv{:,sel{1}},inwv{:,sel{2}}*factor,inwv{:,sel{3}}};
        
        wvtime = inwvdst.RowNames;
        dsp = setDSproperties();
        wvdst(i) = dstable(indata{:},'RowNames',wvtime,'DSproperties',dsp);
    end
end


%%
function selection = getInputUI(vardesc,titletxt)
    %define inputgui for the selection of variables
    % to see field defintions use >>help inputgui
    inp.fields = {'Sig. wave height','Peak period','Mean direction'};                              
    inp.style = {'popupmenu','popupmenu','popupmenu'};
    inp.defaults = {vardesc,vardesc,vardesc};

    selection = inputgui('FigureTitle',titletxt,...
                         'InputFields',inp.fields,...
                         'Style',inp.style,...
                         'ActionButtons', {'Select','Cancel'},...
                         'DefaultInputs',inp.defaults,...
                         'PromptText','Select variables to use');
end


%%
function dsp = setDSproperties()
    %define the metadata properties for the data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'Hs','Tp','Dir'},...
        'Description',{'Significant wave height',...
                'Peak period','Wave direction'},...
        'Unit',{'m','s','deg'},...
        'Label',{'Wave height (m)','Wave period (s)','Wave direction (deg)'},...
        'QCflag',repmat({'raw'},1,3)); 
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