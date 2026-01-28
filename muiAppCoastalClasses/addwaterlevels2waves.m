function [tsdst,meta] = addwaterlevels2waves(wvdst,mobj,meta)
    %
%-------class function help------------------------------------------------------
% NAME
%   addwaterlevels2waves.m
% PURPOSE
%   add water levels to a selected wave dataset. Used by waveModels and
%   ctWaveData when being used for derivative models such as runup in
%   CT_WaveModels.
% USAGE
%   [tsdst,meta] = addwaterlevels2waves(wvdst,mobj,meta);
% INPUTS
%   wvdst - wave dataset dstable
%   mobj - model instance
%   meta - meta data struct
% RESULTS
%   tsdst - updated dstable with water levels added if available
%   meta - updated meta data with details of water level record used
% NOTES
%   used by waveModels and ctWaveData
%
% Author: Ian Townend
% CoastalSEA (c)Jan 2026
%--------------------------------------------------------------------------
%
    %add water levels to a wave dataset - can be model or imported data
    muicat = mobj.Cases;
    wvtime = wvdst.RowNames;      %timesteps in wave record

    wlclassprops = {'ctWaterLevelData','ctTidalAnalysis','muiUserModel'};                                              
    promptxt = 'Select input water level data set (Cancel to use SWL=0):';           
    [wl_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                     'CaseClass',wlclassprops,'ListSize',[300,100]);                                 
    swl = zeros(size(wvtime));               
    if ok<1 || isempty(wl_crec)
        getdialog('Using SWL=0');
        inputxt = sprintf('%s, 0mOD used for water level',meta.inptxt);
        wl_crec = 0;     %assign a null value if no water level data available
    else
        wldst = getDataset(muicat,wl_crec,1); 
        %check that there is water level data for period of interest
        [idst,idnd] = ts2_endpoints_in_ts1(wvdst(1),wldst);
        if isempty(idst)
            getdialog('Data do not overlap. Using SWL=0');
            inputxt = sprintf('%s, 0mOD used for water level',meta.inptxt);
            wl_crec = 0; %assign a null value if no water level data available
        else 
            %select a variable from the water level dataset
            varnames = wldst.VariableNames;
            idx = 1;
            if length(varnames)>1
                [idx,ok] = listdlg('Name','WL options', ...
                    'PromptString','Select a variable:', ...
                    'SelectionMode','single','ListSize',[200,100],...
                    'ListString',varnames);
                if ok<1, idx = 1; end
            end
            wldata = wldst.(varnames{idx});
            wltime = wldst.RowNames;
            swltime = wvtime(idst:idnd);                                        
            %now interpolate water levels onto wave height times
            swl(idst:idnd,1) = interp1(wltime,wldata,swltime,'linear','extrap');
            swl(isnan(swl)) = 0;

            inputxt = sprintf('%s, %s used for water levels',...
                                        meta.inptxt,wldst.Description);
        end
    end

    tsdst = wvdst;
    tsdst(1) = addvars(wvdst(1),swl,'NewVariableNames','swl');
    
    meta.inptxt = inputxt;

    %assign the run parameters to the model instance
    if wl_crec~=0
        meta.caserecs = [meta.caserecs,{wl_crec}];    
    end
end