function times = get_profile_times(mobj,caserec)
%
%-------function help------------------------------------------------------
% NAME
%   get_profile_times.m
% PURPOSE
%   %get the composite time intervals for all profiles
% USAGE
%   times = get_profile_times(mobj,caserec)
% INPUTS
%   mobj - handle to model UI
%   caserec - vector of cse record ids of the profiles to be sorted
%   (numeric or logical)  
% OUTPUT
%   times - vector of unique record times used in profiles defined by caserec 
% SEE ALSO
%   used in CT_BeachAnalysis.m and CT_Plots.m
%
% Author: Ian Townend
% CoastalSEA (c)July 2021
%--------------------------------------------------------------------------
%
    npro = length(caserec);            
    times = [];           %time intervals in any of the surveys
    for j=1:npro
       [dst,~] = getDataset(mobj.Cases,caserec(j),1);  %idset=1 
       newtime = dst.RowNames;
       times = unique(vertcat(times,newtime));
    end
end