function [wndst,meta] = extract_wind_data(inwndst,isfetch)
%
%-------function help------------------------------------------------------
% NAME
%   extract_wind_data.m
% PURPOSE
%   Extract AvSpeed,MaxSpeed,Dir from a dataset that does not use default naming
%   convention
% USAGE
%   [wndst,meta] = extract_ind_data(inwndst)
% INPUTS
%   inwndst - dstable of data to use for data selection
%   isfetch - true if fetch is also required
% OUTPUT
%   wndst - dstable with default naming of AvSpeed, MaxSpeed, Dir
%   meta -  array of variable names for slected variables
% SEE ALSO
%   used in ctWaveSpectra
%   
% Author: Ian Townend
% CoastalSEA (c) Oct 2025
%--------------------------------------------------------------------------
%  
    wndst = []; meta = [];
    vardesc = inwndst.VariableDescriptions;
    if isfield(inwndst.MetaData,'zw')
        zw = num2str(inwndst.MetaData.zw);
    else
        zw = '10';
    end

    %call UI to select all required fields
    sel = getInputUI(vardesc,zw,isfetch);
    if isempty(sel), return; end

    inwn = inwndst.DataTable;
    indata = {inwn{:,sel{1}},inwn{:,sel{2}},inwn{:,sel{3}}};
    wntime = inwndst.RowNames;
    dsp = setDSproperties();
    wndst = dstable(indata{:},'RowNames',wntime,'DSproperties',dsp);
    wndst.Description = inwndst.Description;
    wndst.MetaData.zw = str2double(sel{4});
    wndst.MetaData.Fetch = str2double(sel{5});
    %assign metadata of selection
    meta.selection(1,:) = [sel{1:3}];
    meta.inputs(1,:) = vardesc([sel{1:3}]);
end


%%
function selection = getInputUI(vardesc,zw,isfetch)
    %define inputgui for the selection of variables
    % to see field defintions use >>help inputgui
    inp.fields = {'Av. Speed','Direction','Max. Speed','Elevation aMSL'};                              
    inp.style = {'popupmenu','popupmenu','popupmenu','edit'};
    inp.defaults = {vardesc,vardesc,vardesc,zw};
    defsel = {1,2,3,[],[]};

    if isfetch
        inp.fields = [inp.fields,{'Fetch (m)'}];
        inp.style = [inp.style,{'edit'}];
        inp.defaults = [inp.defaults,{'10000'}];
    end

    selection = inputgui('FigureTitle','Wind data',...
                         'InputFields',inp.fields,...
                         'Style',inp.style,...
                         'ActionButtons', {'Select','Cancel'},...
                         'DefaultInputs',inp.defaults,...
                         'UserData',defsel,...
                         'PromptText','Select variables to use');
end

%%
%--------------------------------------------------------------------------
% dataDSproperties
%--------------------------------------------------------------------------
function dsp = setDSproperties()
    %define the variables in the dataset
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'AvSpeed','Dir','MaxSpeed'},...              
        'Description',{'Mean wind speed','Mean wind direction','Maximum wind speed'},...
        'Unit',{'m/s','deg','m/s'},...
        'Label',{'Wind speed (m/s)','Wind direction (deg)','Wind speed (m/s)'},...
        'QCflag',repmat({'raw'},1,3)); 
    dsp.Row = struct(...
        'Name',{'Time'},...
        'Description',{'Time'},...
        'Unit',{'h'},...
        'Label',{'Time'},...
        'Format',{'dd-MM-yyyy HH:mm:ss'});        
    dsp.Dimensions = struct(...    
        'Name',{'Position'},...
        'Description',{'Latitude and Longitude'},...
        'Unit',{'deg'},...
        'Label',{'Latitude and Longitude'},...
        'Format',{''});         
end