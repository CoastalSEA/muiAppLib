function gd_colormap(zlimits)
%
%-------function help------------------------------------------------------
% NAME
%   gd_colormap.m
% PURPOSE
%   check if Mapping toolbox is installed to use land/sea colormap, or call
%   cmap_selection if not available
% USAGE
%    gd_colormap(zlimits)
% INPUTS
%   zlimits - minimum and maximum elevations [minz,maxz]
% OUTPUT
%   sets the colormap based on land/sea map, or user selection
% SEE ALSO
%   cmap_selection.m
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    if license('test','MAP_Toolbox')
        cmapsea = [0,0,0.2;  0,0,1;  0,0.45,0.74;  0.30,0.75,0.93; 0.1,1,1];
        cmapland = [0.95,0.95,0.0;  0.1,0.7,0.2; 0,0.4,0.2; 0.8,0.9,0.7;  0.4,0.2,0];
        demcmap(zlimits,128,cmapsea,cmapland) %mapping toolbox
    else
        cmap = cmap_selection;
        if isempty(cmap), cmap = 'parula'; end
        colormap(cmap)        
    end
end