function shp = gd_readshapefile(path,filename)
%
%-------function help------------------------------------------------------
% NAME
%    gd_readshapefile.m
% PURPOSE
%   read the x and y coordinates from a shape file. Lines are concatenated
%   and separated by NaNs in single x and y vectors. Suitable for reading
%   boundaries or sections into a single array.
% USAGE
%   shp gd_readshapefile(path,filename);
% INPUTS
%   path - full path to shapefile to be loaded
%   filename - name of shaperfile to be loaded
% OUTPUT
%   shp - struct with x and y fields from the XY data in the shape file.
% NOTES
%   checks whether Mapping toolbox is available and if not checks whether
%   m_shaperead from the M-Map1.4 toolbox is available. Returns empty shp
%   if neither toolboxes are available. Mapping toolbox return 2xN row
%   vectors and M-Map returns Nx2 column vectors.
% SEE ALSO
%   used in edb_surfacearea_table and loading of bourndaries and
%   cross-sections in EstuaryDB
%   
% Author: Ian Townend
% CoastalSEA (c) Jan 2025
%--------------------------------------------------------------------------
%
    istoolbox = license('test','MAP_Toolbox');   %toolbox is licensed to use
    if istoolbox
        addons = matlab.addons.installedAddons;
        istoolbox = any(matches(addons.Name,'Mapping Toolbox')); %toolbox is installed
    end  

    [~,shapename,~] = fileparts(filename);
    shp.x = []; shp.y = [];
    if istoolbox
        Shp = shaperead([path,shapename]);   %requires Mapping toolbox
        for i=1:length(Shp)
            shp.x = [shp.x,Shp(i).X];    %lines separated by NaNs
            shp.y = [shp.y,Shp(i).Y];    %
        end     
    elseif isfile(which('m_shaperead.m'))    %isfile only works for specified or current path        
        Shp = m_shaperead([path,shapename]); %use M-Map function instead
        for i=1:length(Shp.ncst)    
            if isnan(Shp.ncst{i,1}(end,1))
                shp.x = [shp.x;Shp.ncst{i,1}(:,1)]; 
                shp.y = [shp.y;Shp.ncst{i,1}(:,2)];                
            else
                shp.x = [shp.x;Shp.ncst{i,1}(:,1);NaN];    %the added NaN fixes a problem when
                shp.y = [shp.y;Shp.ncst{i,1}(:,2);NaN];    %using the polygon in insidepoly
            end
        end
    else
        shp = [];
        errordlg('Unable to load Shape file\nCheck Mapping toolbox or M-Map is installed'); 
        return;
    end
end