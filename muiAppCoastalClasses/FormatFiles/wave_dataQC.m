function dst = wave_dataQC(dst)
%
%-------function help------------------------------------------------------
% NAME
%   wave_dataQC.m
% PURPOSE
%   generic quality control for a wave data timeseries 
% USAGE
%   dst = wave_dataQC(dst)
% INPUTS
%   dst - the obj passed to the function is a dstable
% OUTPUT
%   dst - dst having been passed though quality checks
% NOTES
%   Applies user defined limit to wave period and 1/7 steepness limit to waves
% SEE ALSO
%   wave_cco_format.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%
    Hs = dst.Hs;
    Tp = dst.Tp;
    Dir = dst.Dir;

    maxTp = 25;
    %get input parameters from user
    prompt = {'Maximum value of Wave Peak Period (s):'};
    title = 'Define limiting peak period';
    numlines = 1;
    defaultvalues{1} = num2str(maxTp);
    useInp=inputdlg(prompt,title,numlines,defaultvalues);
    if isempty(useInp), return; end %user cancelled
    maxTp = str2double(useInp{1});  

    idx = false(size(Hs));
    idx = idx | (Tp>maxTp);
    %originally used Hs as 1/18 and Hmax as 1/16 but changed to
    %single wave limit of 1/7 = 0.14 for both
    Steep = Hs./(1.56*Tp.^2);  %denominator is gT^2/2pi
    idx = idx | (Steep>0.14);         
    idx = idx | (Tp<=0);   %Negative wave periods

    Dir(Dir<0) = Dir(Dir<0)+360;
    Dir(Dir>360) = Dir(Dir>360)-360;

    Hs(idx) = NaN;    
    Tp(idx) = NaN;
    Dir(idx) = NaN;

    hw = waitbar(0, 'Loading data. Please wait');

    dst.Hs = Hs;  
    dst.Tp = Tp;  
    dst.Dir = Dir;
    varname = {'Hs','Tp','Dir'};
    for i=1:3
        idx = strcmp(dst.VariableNames,varname{i});
        dst.VariableQCflags{idx} = 'qc';
    end

    waitbar(1); 
    close(hw);
end