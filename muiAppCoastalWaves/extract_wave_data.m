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
        meta.selection = [];
        meta.seastate = [];
        return; 
    end 

    %handle unimodal and multimodal input data options 
    if nargin<2
        nvar = inputdlg({'Number of components, 1-3'},'Multimodal',1,{'1'});
    end
    getdialog('Default selection for Copernicus data is the combined sea state',...
                                                                     [],1);
    if isempty(nvar)
        wvdst = copy(inwvdst); 
        meta = msgbox('No selection made'); return
    else
        nvar = str2double(nvar{1});
    end
    %
    switch nvar
        case 1
            t_txt = {'Wave data'};
        case 2
            t_txt = {'Wind-wave','Swell wave'};
        case 3
            t_txt = {'Wind-wave','Primary swell wave','Secondary swell wave'};
        otherwise
            meta = msgdlg('Invalid number of components'); return
    end


    inwv = inwvdst.DataTable;
    for i=1:nvar
        %call UI to select all required fields
        sel = getInputUI(vardesc,t_txt,nvar,i);
        if isempty(sel), wvdst = []; return; end
        
        %factor = 1;
        % if contains(vardesc{sel{2}},'mean period') 
        %     %indata = setWavePeriods(inwvdst,sel{2});
        %     %Tp=var.no.17
        %     indata = {inwv{:,sel{1}},inwv{:,17},inwv{:,sel{3}},...
        %         inwv{:,sel{2}},inwv{:,15},inwv{:,16}};
        %     ismean = true;
            
        % elseif sel{2}==17 && nvar==1
        %     T1 = getWavePeriod(inwvdst);
        %     indata = {inwv{:,sel{1}},inwv{:,17},inwv{:,sel{3}},T1,inwv{:,15},inwv{:,16}};
        %     ismean = true;
        % else   
            %extract data selected            
            indata = {inwv{:,sel{1}},inwv{:,sel{2}},inwv{:,sel{3}}};
            ismean = false;
        % end
        wvtime = inwvdst.RowNames;
        dsp = setDSproperties(ismean);
        wvdst(i) = dstable(indata{:},'RowNames',wvtime,'DSproperties',dsp); %#ok<AGROW>
        wvdst(i).Description = inwvdst.Description; %#ok<AGROW>
        %assign metadata of selection
        meta.selection(i,:) = [sel{:}];
        meta.inputs(i,:) = vardesc([sel{:}]);        
    end
    
    if strcmp(varnames{1},'VHM0') &&strcmp(varnames{17},'VTPK') %test for Copernicus data!
        dsp = setDSproperties(true);
        T1 = getWavePeriod(inwvdst);
        %table of total sea state parameters Hs,Tp,Dir,T1(est),T2,T10
        indata = {inwv{:,1},inwv{:,17},inwv{:,5},T1,inwv{:,15},inwv{:,16}};                
        meta.seastate = dstable(indata{:},'RowNames',wvtime,'DSproperties',dsp);
    else
        meta.seastate = [];
    end
end

%%
function selection = getInputUI(vardesc,titletxt,nvar,idx)
    %define inputgui for the selection of variables
    % to see field defintions use >>help inputgui
    inp.fields = {'Sig. wave height','Peak period','Mean direction'};                              
    inp.style = {'popupmenu','popupmenu','popupmenu'};
    inp.defaults = {vardesc,vardesc,vardesc};
    defsel = {1,2,3};
    if numel(vardesc)==17
        defsel = getCopComponents(nvar,idx);
    end     
    selection = inputgui('FigureTitle',titletxt{idx},...
                         'InputFields',inp.fields,...
                         'Style',inp.style,...
                         'ActionButtons', {'Select','Cancel'},...
                         'DefaultInputs',inp.defaults,...
                         'UserData',defsel,...
                         'PromptText','Select variables to use');
end

%%
function defsel = getCopComponents(nvar,idx)
    %selection indices for Copernicus 1 to 3 components
    % 1  'Sea surface wave significant height',...
    % 2  'Sea surface primary swell wave significant height',...
    % 3  'Sea surface secondary swell wave significant height',...
    % 4  'Sea surface wind wave significant height',...
    % 5  'Sea surface wave from direction',...
    % 6  'Sea surface primary swell wave from direction',...
    % 7  'Sea surface secondary swell wave from direction',...
    % 8  'Sea surface wind wave from direction',...
    % 9  'Sea surface wave from direction at variance spectral density maximum',...
    % 10 'Sea surface wave stokes drift x velocity',...
    % 11 'Sea surface wave stokes drift y velocity',...
    % 12 'Sea surface primary swell wave mean period',...
    % 13 'Sea surface secondary swell wave mean period',...
    % 14 'Sea surface wind wave mean period',...
    % 15 'Sea surface wave mean period from variance spectral density second frequency moment',...
    % 16 'Sea surface wave mean period from variance spectral density inverse frequency moment',...
    % 17 'Sea surface wave period at variance spectral density maximum'
    switch nvar
        case 1 %Copernicus combined sea state parameters
            defsel = {1,17,5}; 
        case 2 %Copernicus wind and primary swell sea state parameters
            if idx==1
                defsel = {4,14,8}; 
            else
                defsel = {2,12,6}; 
            end
        case 3 %Copernicus wind and both swell sea state parameters
            if idx==1
                defsel = {4,14,8}; 
            elseif idx==2
                defsel = {2,12,6}; 
            else
                defsel = {3,13,7}; 
            end
    end
end

%%
function T1_tot = getWavePeriod(inwvdst)
    %use energy‑weighted harmonic mean of the sea state components to estimate 
    %the value of T1 for the totla sea state - bespoke for Copernicus wave data
    inwv = inwvdst.DataTable;

    T1(:,1) = inwv{:,12};  %primary swell mean period
    T1(:,2) = inwv{:,13};  %primary swell mean period
    T1(:,3) = inwv{:,14};  %wind wave mean period

    Hs(:,1) = inwv{:,2};
    Hs(:,2) = inwv{:,3};
    Hs(:,3)= inwv{:,4};

    % Valid components: both Hs>0 and T1>0
    valid = (Hs > 0) & (T1 > 0);

    % Preallocate numerator and denominator
    num = sum(Hs.^2 .* valid, 2);   % Σ Hs_i^2 over valid components

    % Compute denominator safely: only divide where valid
    den = zeros(size(num));
    for j = 1:size(Hs,2)
        idx = valid(:,j);           % rows where component j is valid
        den(idx) = den(idx) + (Hs(idx,j).^2 ./ T1(idx,j));
    end

    % Final total T1
    T1_tot = num ./ den;

    % Handle rows where no components are valid
    T1_tot(den == 0) = 0;           % or NaN if preferred
end

%%
function dsp = setDSproperties(ismean)
    %define the metadata properties for the data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
    %struct entries are cell arrays and can be column or row vectors
    if ismean
        dsp.Variables = struct(...
            'Name',{'Hs','Tp','Dir','T1','T2','T10'},...
            'Description',{'Significant wave height',...
                    'Peak period','Wave direction','Mean period',...
                    'Second moment period','Inverse moment period'},...
            'Unit',{'m','s','deg','s','s','s'},...
            'Label',{'Wave height (m)','Wave period (s)',...
                     'Wave direction (deg)','Wave period (s)',...
                     'Wave period (s)','Wave period (s)'},...
            'QCflag',repmat({'raw'},1,6)); 
    else
        dsp.Variables = struct(...
            'Name',{'Hs','Tp','Dir'},...
            'Description',{'Significant wave height',...
                    'Peak period','Wave direction'},...
            'Unit',{'m','s','deg'},...
            'Label',{'Wave height (m)','Wave period (s)','Wave direction (deg)'},...
            'QCflag',repmat({'raw'},1,3)); 
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