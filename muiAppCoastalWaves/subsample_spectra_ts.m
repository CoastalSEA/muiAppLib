function newdst = subsample_spectra_ts(dst,mobj,method,tol)
%-------function help------------------------------------------------------
% NAME
%   subsample_spectra_ts.m
% PURPOSE
%   create a timeseries by interpolating a spectrum timeseries using times
%   from another timeseries
% USAGE
%   newdst = subsample_spectra_ts(dst,mobj,method,tol)
% INPUT
%   dst - dstable holding spectral data to the CCO/spt format
%   mobj - handle to muiModelUI instance to allow access to data
%   method - interpolation method used in interp1 (optional, default = linear)
%            when method defined as 'none', function only selects data with a
%            date match.
%   tol - only required if method ='none' and the datetimes are to be
%         matched using a tolerance. tol in seconds
% OUTPUT
%   newdst - new dstable is a copy of the source dstable with the datasets
%            interpolated or subsampled to the new time intervals
%
% Author: Ian Townend
% CoastalSEA (c)March 2026
%--------------------------------------------------------------------------
%
    muicat = mobj.Cases;
    if nargin<3
        method = 'linear';
        tol = [];
    elseif nargin<4
        tol = [];
    end  

    %get the dataset used for the subsample time intervals
    promptxt = 'Select dataset to define sub-sample time intervals';
    [caserec,isok] = selectRecord(muicat,'PromptText',promptxt,...
                                                    'ListSize',[300,100]);    
    if isok<1, newdst = []; return; end %user cancelled
    cobj = getCase(muicat,caserec);
    if isempty(cobj), return; end

    dnames = fields(cobj.Data);
    if numel(dnames)>1
        [~,idd] = selectDataset(muicat,cobj);
    else
        idd = 1;
    end
    timedst =  cobj.Data.(dnames{idd});          
    

    [newtime,idt] = getTimes2Use(dst,    timedst,method,tol);
    clear cobj timedst
    
    adst = copy(dst);
    for j=1:2
        subtable = dst(j).DataTable(idt,:);  
        adst(j).DataTable = subtable;
        adst(j).RowNames = newtime; %update times in case there is an offset
    end
    newdst.sptSpectrum = adst(1);
    newdst.sptProperties = adst(2);
end

%%
function [newtime,idt] = getTimes2Use(dst,timedst,method,tol)
    %get the times to use for sampling
    vartime = dst.RowNames;            %dataset being sampled
    newtime = timedst.RowNames;        %dataset to define new times
    
    if strcmp(method,'none')
        if isempty(tol)
            inp = inputdlg({'Use a tolerance (s)? [0 for exact match]'},'Subsample',1,{'0'});
            if isempty(inp) || strcmp(inp{1},'0')
                tol = [];
            else
                tol = seconds(str2double(inp{1}));
            end
        else
            tol = seconds(tol);
        end
        %
        if isempty(tol)
            [idn,idv] = ismember(newtime, vartime);
            idt = find(idv(idv>0));%??????
            newtime = newtime(idn);
        else
            D = abs(newtime - vartime');      % duration matrix
            [minDiff, idx] = min(D, [], 2);
            tf = minDiff <= tol;
            idt = idx(tf);
            newtime = newtime(tf);
        end

    else        
        idt = [];
    end
end
