function varargout = gd_colormap(zlimits,isselect)
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
%   isselect - alllow user to select a color map - optional default is false
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
     if nargin<1, isselect = false; end
     iscaxis = isMATLABReleaseOlderThan("R2022a");
     lims = [];

    if license('test','MAP_Toolbox')
        cmapsea = [0,0,0.2;  0,0,1;  0,0.45,0.74;  0.30,0.75,0.93; 0.1,1,1];
        cmapland = [0.95,0.95,0.0;  0.1,0.7,0.2; 0,0.4,0.2; 0.8,0.9,0.7;  0.4,0.2,0];
        [cmap,lims] = demcmap(zlimits,128,cmapsea,cmapland); %mapping toolbox
    elseif isselect
        cmap = cmap_selection;
    else
        cmap = 'parula';      
    end

    if isempty(lims) && iscaxis
        lims = caxis; %#ok<CAXIS>
    elseif isempty(lims)
        lims = clim;  %v2022a
    end

    if iscaxis
        caxis(lims) %#ok<CAXIS> 
    else
        clim(lims);
    end

    colormap(cmap);

    %return output if request
    output = {cmap,lims}; %handles to cross and and circle
    for k=1:nargout
        varargout{k} = output{k}; %#ok<AGROW> 
    end
end