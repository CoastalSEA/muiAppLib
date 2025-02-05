function varargout = gd_colormap(zlimits)
%
%-------function help------------------------------------------------------
% NAME
%   gd_colormap.m
% PURPOSE
%   check if Mapping toolbox is installed to use land/sea colormap, or call
%   cmap_selection if not available
% USAGE
%    [cmap,clim] = gd_colormap(zlimits);
% INPUTS
%   zlimits - minimum and maximum elevations [minz,maxz]
% OUTPUT
%   sets the colormap based on land/sea map, or user selection and also 
%   returns the colormap matrix and color axis limit vector
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
        [cmap,lims] = demcmap(zlimits,128,cmapsea,cmapland); %mapping toolbox
    else
        cmap = cmap_selection;
        if isempty(cmap), cmap = 'parula'; end           
        try
            lims = clim;  %v2022a
        catch
            lims = caxis; %#ok<CAXIS> 
        end        
    end
    clim(lims);
    colormap(cmap);

    %return output if request
    output = {cmap,lims}; %handles to cross and and circle
    for k=1:nargout
        varargout{k} = output{k}; %#ok<AGROW> 
    end
end