function Wz = gd_plan_form(grid,wl)
%
%-------function help------------------------------------------------------
% NAME
%   gd_plan_form.m
% PURPOSE
%   compute planform variation along the x-axis at specified planar levels 
%   for impored grids.
% USAGE
%   Wz = gd_plan_form(grid,zwl)
% INPUTS
%   grid - struct of x,y,z,t values that define grid 
%   wl - struct of levels for zhw, zmt, zlw
% OUTPUTS
%   Wz - table of widths for each level, Whw, Wmt, Wlw (to match model
%        output, e.g. cf_exp_models.m)
% NOTES
%   function is used to generate plan form data set for grids that have 
%   been imported. Grids from models output yz with the same output format.
% SEE ALSO
%   called by cf_valley_model as part of ChannelForm model and
%   addFormProperties in GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    nx = length(grid.x);
    ich = gd_basin_indices(grid);  %account for offset to mouthand head (if defined)
    if isscalar(wl.zhw)
        zwl(:,1) = ones(nx,1)*wl.zhw;          
        zwl(:,3) = ones(nx,1)*wl.zlw;
        zwl(:,2) = (zwl(:,1)+zwl(:,3))/2;           %mean tide level(m)  
    else 
        zwl(:,1) = wl.zhw;
        zwl(:,2) = wl.zmt;
        zwl(:,3) = wl.zlw;
    end
    nz = 3;
    Wz = zeros(nx,nz);
    dely = abs(grid.y(2)-grid.y(1)); %grid interval

    zi = grid.z;
    zmax = max(zi,[],'all');
    if any(zmax<zwl)
        %grid does not extend to highest water level
        warndlg('Grid does not does not extend to high water. Plan form not set');
        %yz = num2cell(yz',2)';      %formatted to load into dstable
        Wz = [];
        return;
    end

    for jz=1:nz
        for ix=ich
            zij = zi(ix,:);
            zij(zij>zwl(ix,jz)) = NaN;         %set values above zwl to NaN
            Wz(ix,jz) = sum(~isnan(zij)*dely); %width at zwl
        end
    end
    
    Wz = num2cell(Wz',2)';           %formatted to load into dstable
    Wz = table(Wz{:},'VariableNames',{'Whw','Wmt','Wlw'});
%     dsp = setDSproperties();
%     Wzdst = dstable(Wz{:},'RowNames',grid.t,'DSproperties',dsp);
%     Wzdst.Dimensions.X = grid.x;  %NB uses full grid
end
%%
% function dsp = setDSproperties()
%     %define the metadata properties for the demo data set
%     dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
%     %define each variable to be included in the data table and any
%     %information about the dimensions. dstable Row and Dimensions can
%     %accept most data types but the values in each vector must be unique
%     % *note: these would be better in the gd_property functions so
%     % *that the defintion is in the same file as the assignment
% 
%     %struct entries are cell arrays and can be column or row vectors
%     dsp.Variables = struct(...
%         'Name',{'Whw','Wmt','Wlw'},...                  
%         'Description',{'HW Width','MT Width','LW Width'},...
%         'Unit',{'m','m','m'},...
%         'Label',{'Width (m)','Width (m)','Width (m)'},...
%         'QCflag',repmat({'model'},[1,3]));  
%     dsp.Row = struct(...
%                     'Name',{'Time'},...
%                     'Description',{'Time'},...
%                     'Unit',{'y'},...
%                     'Label',{'Time (yr)'},...
%                     'Format',{'y'});         
%     dsp.Dimensions = struct(...    
%                     'Name',{'X'},...
%                     'Description',{'Distance'},...
%                     'Unit',{'m'},...
%                     'Label',{'Distance (m)'},...
%                     'Format',{'-'});  
% end  