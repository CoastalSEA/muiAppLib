function [sortedE,sortedN,idd] = shore_profile_order(mobj,caserec)
%
%-------function help------------------------------------------------------
% NAME
%   shore_profile_order.m
% PURPOSE
%   sort the profile order based on the E,N of the base point - min(Chainage)
% USAGE
%   [sortedE,sortedN,idd] = shore_profile_order(mobj,caserec)
% INPUTS
%   mobj - handle to model UI
%   caserec - vector of cse record ids of the profiles to be sorted
%   (numeric or logical)
% OUTPUT
%   sortedE - vector of sorted x-coordinates
%   sortedN - vector of sorted y-coordinates
%   idd - indices of sorting
% NOTES
%   sorting is based on distance from the point with the minimum Easting
%   if data set has no Eastings and Northings returns existing order
%
% Author: Ian Townend
% CoastalSEA (c)July 2021
%--------------------------------------------------------------------------
%     
    sortedE = []; sortedN = [];
    npro = length(caserec);   %number of profiles
    Es = NaN(1,npro); Ns = Es;
    for k=1:npro
        [dst,~] = getDataset(mobj.Cases,caserec(k),1);  %idset=1
        %check whether data set has Eastings and Northings
        if ~isprop(dst,'Eastings')
            idd = 1:npro;
            return; 
        end
        
        E = dst.Eastings;
        N = dst.Northings; 
        Ch = dst.Chainage;
        [~,idminCh] = min(Ch(:,1),[],'omitnan'); %minimum in first column
        Es(1,k) = E(idminCh,1);
        Ns(1,k) = N(idminCh,1); 
    end
    %sort profiles into alongshore order
    [sortedE,sortedN,idd] = sortENdata2line(Es,Ns);
end