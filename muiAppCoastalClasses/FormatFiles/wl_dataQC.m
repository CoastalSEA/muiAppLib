function dst = wl_dataQC(dst)
%
%-------function help------------------------------------------------------
% NAME
%   wl_dataQC.m
% PURPOSE
%   generic quality control for a water level data timeseries 
% USAGE
%   dst = wl_dataQC(dst)
% INPUTS
%   dst - the obj passed to the function is ta dstable
% OUTPUT
%   dst - dst having been passed though quality checks
% NOTES
%   Applies user defined minimum and maximum limits to water levels
% SEE ALSO
%   wl_cco_format.m and wl_bodc_format.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%
    WLOD = dst.WLOD;
    WLCD = dst.WLCD;

    varnames = dst.VariableNames;
    if any(strcmp(varnames,'WLFlag'))
        WLflag = dst.WLFlag;
        %BODC NSTLF uses M which is converted to 9 when data is loaded
        idx = WLflag==9;  
        %CCO uses the following flags for water level elevation
        %1 (correct value)
        %2 (interpolated)
        %3 (doubtful)
        %4 (isolated spike or wrong value)%
        %5 (correct but extreme)
        %6 (reference change detected)
        %7 (constant value)
        %8 (out of range)
        %9 (missing)
        idx = idx | (WLflag>2 & WLflag<5);
        idx = idx | WLflag==8;
    else
        idx = false(size(WLOD));
    end            
    WLOD(idx) = NaN;

    figure('Name','QC Plot','Tag','PlotFig');
    plot(WLOD);                  
    maxWL = ceil(max(WLOD));
    minWL = floor(min(WLOD));
    %get input parameters from user
    prompt = {'Maximum water level (mOD):','Minimum water level (mOD):'};
    title = 'Define limiting water levels';
    numlines = 1;
    defaultvalues = {num2str(maxWL),num2str(minWL)};
    useInp=inputdlg(prompt,title,numlines,defaultvalues);
    if isempty(useInp), return; end %user cancelled
    maxWL = str2double(useInp{1});
    minWL = str2double(useInp{2});    
    idx = idx | (WLOD>maxWL);
    idx = idx | (WLOD<minWL);

    WLOD(idx) = NaN;
    WLCD(idx) = NaN;

    hw = waitbar(0, 'Loading data. Please wait');            

    dst.WLOD = WLOD;             
    dst.WLCD = WLCD;            
    dst.VariableQCflags(1:2) = repmat({'qc'},1,2);
    
    waitbar(1); 
    close(hw);
end