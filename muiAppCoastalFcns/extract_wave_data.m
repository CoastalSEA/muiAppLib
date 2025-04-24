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
    wvdst = [];

    vardesc = inwvdst.VariableDescriptions;
    %select wave height
    select_Hs = inputSelection(vardesc,'Wave height');
    if isempty(select_Hs), return; end
    %select wave period
    select_Tp = inputSelection(vardesc,'Wave period');
    if isempty(select_Tp), return; end
    %select wave direction
    select_Dir = inputSelection(vardesc,'Wave direction');
    if isempty(select_Dir), return; end

    inwv = inwvdst.DataTable;
    indata = {inwv{:,select_Hs},inwv{:,select_Tp},inwv{:,select_Dir}};
    
    wvtime = inwvdst.RowNames;
    dsp = setDSproperties();
    wvdst = dstable(indata{:},'RowNames',wvtime,'DSproperties',dsp);
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

%%
function selection = inputSelection(vardesc,varname)
    promptxt = sprintf('Select variable for %s',varname);
    selection = listdlg("PromptString",promptxt,"ListSize",[320,160],...
                        "ListString",vardesc,"SelectionMode","single");
end