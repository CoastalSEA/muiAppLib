function [wvdst,meta] = extract_wave_data(inwvdst,nvar)
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
%   wvdst - array of dstables with default naming of Hs, Tp and Dir for
%           each component selected (N=1-3)
%   meta -  array of variable names for each component defined (N=1-3)
% SEE ALSO
%   used in ctWaveModel and WRM_WaveModel
%   
% Author: Ian Townend
% CoastalSEA (c) April 2025
%--------------------------------------------------------------------------
%  
    varnames = inwvdst.VariableNames;
    vardesc = inwvdst.VariableDescriptions;

    if sum(ismatch(varnames,{'Hs','Tp','Dir'}))==3 
        wvdst = copy(inwvdst);
        idel = ~ismatch(varnames,{'Hs','Tp','Dir'});
        wvdst = removevars(wvdst,varnames(idel));
        meta.inputs(1,:) = vardesc(ismatch(varnames,{'Hs','Tp','Dir'}));
        return; 
    end 

    %handle unimodal and multimodal input data options 
    if nargin<2
        nvar = inputdlg({'Number of components, 1-3'},'Multimodal',1,{'1'});
    end
    getdialog('Default selection for Copernicus data is the combined sea state');
       
    if isempty(nvar)
        wvdst = copy(inwvdst); 
        meta = msgdlg('No selection made'); return
    else
        Ntype = str2double(nvar{1});
    end
    %
    switch Ntype
        case 1
            t_txt = {'Wave data'};
        case 2
            t_txt = {'Wind-wave','Swell wave'};
        case 3
            t_txt = {'Wind-wave','Primary swell wave','Secondary swell wave'};
        otherwise
            meta = msgdlg('Invalid number of components'); return
    end

    for i=1:Ntype
        %call UI to select all required fields
        sel = getInputUI(vardesc,t_txt{i});
        if isempty(sel), wvdst = []; return; end

        factor = 1;
        if contains(vardesc{sel{2}},'mean')
            factor = 1.2; %scale mean period to peak period. this value            
        end               %is for Jonswap with gamma=3.3
        %extract data selected
        inwv = inwvdst.DataTable;
        indata = {inwv{:,sel{1}},inwv{:,sel{2}}*factor,inwv{:,sel{3}}};
        wvtime = inwvdst.RowNames;
        dsp = setDSproperties();
        wvdst(i) = dstable(indata{:},'RowNames',wvtime,'DSproperties',dsp); %#ok<AGROW>
        wvdst(i).Description = inwvdst.Description;
        %assign metadata of selection
        meta.selection(1,:) = [sel{:}];
        meta.inputs(i,:) = vardesc([sel{:}]);
    end
end


%%
function selection = getInputUI(vardesc,titletxt)
    %define inputgui for the selection of variables
    % to see field defintions use >>help inputgui
    inp.fields = {'Sig. wave height','Peak period','Mean direction'};                              
    inp.style = {'popupmenu','popupmenu','popupmenu'};
    inp.defaults = {vardesc,vardesc,vardesc};
    defsel = {1,2,3};
    if numel(vardesc)==17,defsel = {1,17,5}; end %Copernicus combined sea state parameters
        
    selection = inputgui('FigureTitle',titletxt,...
                         'InputFields',inp.fields,...
                         'Style',inp.style,...
                         'ActionButtons', {'Select','Cancel'},...
                         'DefaultInputs',inp.defaults,...
                         'UserData',defsel,...
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