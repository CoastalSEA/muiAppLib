function pdst = gd_section_properties(grid,wl,histdst)
%
%-------function help------------------------------------------------------
% NAME
%   gd_section_properties.m
% PURPOSE
%   compute the width, cross-sectional area and prism along channel
% USAGE
%   pdst = gd_section_properties(grid,wl,histdst)
% INPUTS
%   grid - struct of x,y,z,t values that define grid 
%   wl - struct of zhw, zmt and zlw with scalar values at the mouth, or a
%        vector with the same dimension as the x vector, or CF_HydroData instance
%   histdst - elevation histogram, SArea, as a dstable with dimension of x and z 
% OUTPUT
%   pdst - Along channel/x-axis properties as a dstable, including
%            widths at HW, MT and LW 
%            cross-sectional area at HW, MT and LW
%   	     cumulative tidal prism and csa/prism ratio using CSA
%            cumulative surface areas from head at HW, MT and LW 
%            cumulative volumes from head at HW, MT and LW 
%            cumulative tidal prism from head using volumes
%            Dronkers Gamma
%            cumulative storage volumes from head
%            cumulative channel volumes from head
%            hydraulic depth of upstream channel
%            tidal amplitude along channel
% SEE ALSO
%   used in GDinterface, along with gd_basin_hypsometry and
%   gd_gross_properties
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    delx = abs(grid.x(2)-grid.x(1));
    xi = grid.x;
    
    ich = gd_basin_indices(grid);  %account for offset to mouthand head (if defined)
    if isscalar(wl.zhw)
        zhw = ones(size(xi))*wl.zhw;          
        zlw = ones(size(xi))*wl.zlw;
        zmt = (zhw+zlw)/2;           %mean tide level(m)  
    else 
        zhw = wl.zhw;
        zmt = wl.zmt;
        zlw = wl.zlw;
    end
    [vprops,~] = gd_basin_properties(grid,wl,histdst);    

    %initialise arrays
    nint = length(xi); 
    Whw = zeros(1,nint); Wmt = Whw; Wlw = Whw;
    Ahw = zeros(1,nint); Amt = Ahw; Alw = Ahw; pr = Ahw;
    Dhw = zeros(1,nint); Dmt = Dhw; Dlw = Dhw;

    for ix = ich            
        [Whw(1,ix),Ahw(1,ix),Dhw(1,ix)] = interplevel(grid,zhw,ix);
        [Wmt(1,ix),Amt(1,ix),Dmt(1,ix)] = interplevel(grid,zmt,ix);
        [Wlw(1,ix),Alw(1,ix),Dlw(1,ix)] = interplevel(grid,zlw,ix);
        pr(1,ix) = (Ahw(1,ix)-Alw(1,ix));
    end     
    
    PrA =fliplr(cumtrapz(delx,fliplr(pr))); %trapezoidal summation handles ends
    %assign results
    aprops = table(Whw,Wmt,Wlw,Ahw,Amt,Alw,Dhw,Dmt,Dlw,PrA); %formated for dstable input
    props = horzcat(aprops,vprops);

    pdsp = setDSproperties();
    pdst = dstable(props,'RowNames',grid.t,'DSproperties',pdsp);
    pdst.Dimensions.X = grid.x;  %NB uses full grid
    
    % nested function ---------------------------------------------
    function [w_zwl,a_zwl,d_zwl] = interplevel(grid,zwl,idx)
        %get the width and cross-sectional area at distance x and elevation z0
        % grid - struct with grid definition as used in GDinterface
        % zwl - level to use (scalar or vector with size of grid.x)
        % idx - index of section to extract along x-axis
        yi = grid.y;
        dely = abs(yi(2)-yi(1));  %grid interval
        
        zi = grid.z;
        zi = zi(idx,:);
        if ~isscalar(zwl)
            zwl = zwl(idx);
        end         
        
        zi(zi>zwl) = NaN;               %set values above zwl to NaN
        w_zwl = sum(~isnan(zi)*dely);   %width at zwl
        dzwl = zwl-zi;                  %depths below zwl
        d_zwl = max(dzwl,[],'omitnan'); %deepest depth in section to zwl
        idd = ~isnan(dzwl);             %index of valid depths
        if sum(idd)<=2
            a_zwl = 0;
        else
            a_zwl = trapz(yi(idd),dzwl(idd)); %csa below zwl
        end
    end
end      
%%
function dsp = setDSproperties()
    %define the metadata properties for the demo data set
    dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
    %define each variable to be included in the data table and any
    %information about the dimensions. dstable Row and Dimensions can
    %accept most data types but the values in each vector must be unique
    % *note: these would be better in the gd_property functions so
    % *that the defintion is in the same file as the assignment

    %struct entries are cell arrays and can be column or row vectors
    dsp.Variables = struct(...
        'Name',{'Whw','Wmt','Wlw','CSAhw','CSAmt','CSAlw',...
                'Dhw','Dmt','Dlw','PrA',...
                'Shw','Smt','Slw','Vhw','Vmt','Vlw',...
                'PrV','Gamma','Vs','Vc','hyd','amp'},...                  
        'Description',{'HW Width','MT Width','LW Width',...
                       'HW Cross-sectional Area',...
                       'MT Cross-sectional Area',...
                       'LW Cross-sectional Area',...
                       'HW Depth','MT Depth','LW Depth',...
                       'Tidal Prism (CSA)',...
                       'HW Surface area',...
                       'MT Surface area',...
                       'LW Surface area',...
                       'HW Volume','MT Volume','LW Volume',...
                       'Tidal Prism (Hypsometry)',...
                       'Dronkers Gamma',...
                       'Storage volume','Channel volume',...
                       'Hydraulic depth at MTL','Tidal amplitude'},...
        'Unit',{'m','m','m','m2','m2','m2','m3','m2','m2','m2',...
                'm','m','m','m3','m3','m3','m3','-','m3','m3','m','m'},...
        'Label',{'Width (m)','Width (m)','Width (m)',...
                 'Cross-sectioinal area (m^2)',...
                 'Cross-sectioinal area (m^2)',...
                 'Cross-sectioinal area (m^2)',...
                 'Depth (m)','Depth (m)','Depth (m)',...
                 'Tidal prism (m^3)',...
                 'Surface area (m^2)',...
                 'Surface area (m^2)',...
                 'Surface area (m^2)',...
                 'Volume (m^3)','Volume (m^3)','Volume (m^3)',...
                 'Tidal prism (m^3)',...
                 'Dronkers Gamma','Volume (m^3)','Volume (m^3)',...
                 'Hydraulic depth (m)','Amplitude (m)'},...
        'QCflag',repmat({'model'},[1,22]));  
    dsp.Row = struct(...
                    'Name',{'Time'},...
                    'Description',{'Time'},...
                    'Unit',{'y'},...
                    'Label',{'Time (yr)'},...
                    'Format',{'y'});         
    dsp.Dimensions = struct(...    
                    'Name',{'X'},...
                    'Description',{'Distance'},...
                    'Unit',{'m'},...
                    'Label',{'Distance (m)'},...
                    'Format',{'-'});  
end  