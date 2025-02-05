function [subgrid,ixo,iyo] = gd_subgrid(grid,subdomain)
%
%-------function help------------------------------------------------------
% NAME
%   gd_subgrid.m
% PURPOSE
%   extract a subdomain from a grid (xi,yi,zi) and return the extracted
%   grid and the source grid indices of the bounding rectangle 
% USAGE
%   [subgrid,ixo,iyo] = gd_subgrid(grid,subdomain)
% INPUTS
%   grid - x,y,z struct of input grid as x,y vectors and z array
%   subdomain - subdomain to be extracted, defined as [x0,xN,y0,yN]
% OUTPUT
%   subgrid - x,y,z struct of ouput grid for subdomain as x,y vectors and z array
%   ixo,iyo - indices of bounding subdomain in the input xi,yi grid 
% SEE ALSO
%   used in ModelSkill App getUserTools function.
%   previously called getsubgrid. moved to muiuAppGrdFcns in Jan 2025
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%
%%
    xi = grid.x; yi = grid.y; zi = grid.z';
    %
    if isempty(subdomain)
        subdomain = [min(xi),max(xi),min(yi),max(yi)];
    end
    gd_dir = gd_ax_dir(grid);
    %find x indices for bounding rectangle of subgrid
    if gd_dir.x==3 || gd_dir.x==4 %x-axis is reversed
        ix0 = find(xi>=subdomain(2),1,'last');
        ixN = find(xi<=subdomain(1),1,'first');        
    else
        ix0 = find(xi<=subdomain(1),1,'last');
        ixN = find(xi>=subdomain(2),1,'first');        
    end
    %find y indices of bounding rectangle of subgrid
    if gd_dir.y==3    %y-axis is reversed
        iy0 = find(yi>=subdomain(4),1,'last');
        iyN = find(yi<=subdomain(3),1,'first');        
    else
        iy0 = find(yi<=subdomain(3),1,'last');
        iyN = find(yi>=subdomain(4),1,'first');        
    end
    ixo = [ix0,ix0,ixN,ixN,ix0];
    iyo = [iyN,iy0,iy0,iyN,iyN];
    %extract grid values within subdomain    
    subgrid.x = xi(min(ixo):max(ixo));
    subgrid.y = yi(min(iyo):max(iyo));
    subgrid.z = zi(min(iyo):max(iyo),min(ixo):max(ixo))';
end