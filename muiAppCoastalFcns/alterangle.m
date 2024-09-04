function outdir = alterangle(indir)
%
%-------function help------------------------------------------------------
% NAME
%   alterangle.m
% PURPOSE
%   Function to adjust the direction angle of a wind or wave dataset 
% USAGE
%   outdir = alterangle(indir)
% INPUTS
%   indir - vector of directions  
% OUTPUT
%   outdir - vector of modified directions    
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    promptxt = {'Maximum direction shift (deg)',...
                'Angle for zero shift (degTN)',...
                'Angle for maximum shift (degTN)'...
                'Angle for limit of adjustment (degTN)',...
                'Direction of shift (+1 or -1)'};
    defaultxt = {'58','100','200','270','-1'};
    answer = inputdlg(promptxt,'Angle input',1,defaultxt);
    
    alp0 = str2double(answer{1});  %Maximum direction shift (deg)
    Dir0 = str2double(answer{2});  %Angle for zero shift (degTN)
    Dir1 = str2double(answer{3});  %Angle for maximum shift (degTN)
    Dir2 = str2double(answer{4});  %Angle for limit of adjustment (degTN)
    asgn = str2double(answer{5});  %Direction of shift (+1 or -1)
    alp = alp0*(indir-Dir0)/(Dir1-Dir0); %linear scaling of shift betweeb Dir0 and Dir1
    alp(indir>Dir1) = alp0;        %maximum shift beyond Dir1
    alp(indir>Dir2) = alp0;        %maximum shift beyond Dir2
    outdir = indir+asgn*alp;       %adjust direction by +/- alp
end