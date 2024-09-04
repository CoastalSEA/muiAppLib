function gdst = gd_gross_properties(grid,wl,pdst) 
%
%-------function help------------------------------------------------------
% NAME
%   gd_gross_properties.m
% PURPOSE
%   compute the gross properties of a gridded bathymetry
% USAGE
%   gdst = gd_gross_properties(grid,wl,pdst)
% INPUTS
%   grid - struct of x,y,z values that define grid 
%   wl - table of zhw, zmt and zlw with scalar values at the mouth, or a
%        vector with the same dimension as the x vector, or CF_HydroData instance
%   pdst - dstable from gd_section_properties (contains along-channel
%           properties such as width, CSA, Surface area and Volume
% OUTPUT
%   gdst - dstable of the gross properties at mouth, including:
%                Shw, Slw, Vhw, Vlw, PrA, PrV, Gamma, Vs, Vc, Wm, Am, Dm,
%                Wr, Ar, amp, hyd, aoh, VsoVc, PrvAm, SflShw , Lw, La
% NOTES
%   uses the output, pdst, from gd_section_properties for most values
% SEE ALSO
%   used in GDinterface, along with gd_basin_hypsometry and
%   gd_section_properties
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    xi = grid.x;
    [ich,ixM,ixH] = gd_basin_indices(grid); %account for offset to mouth 
                                            %and head (if defined)
    Shw = pdst.Shw(ixM);      %surface area at high water
    Slw = pdst.Slw(ixM);      %surface area at low water
    Vhw = pdst.Vhw(ixM);      %volume at high water
    Vlw = pdst.Vlw(ixM);      %volume at low water
    PrA = pdst.PrA(ixM);      %volume of tidal prism from area of cross-sections
    PrV = pdst.PrV(ixM);      %volume of tidal prism from hypsometry volumes
    Gamma = pdst.Gamma(ixM);  %Dronkers gamma (~1)
    Vs = pdst.Vs(ixM);        %storage volume over intertidal
    Vc = pdst.Vc(ixM);        %channel volume to mean tide level
    VsoVc = Vs/Vc;            %ratio fo storage to channel volumes
    amp = pdst.amp(ixM);      %tidal amplitude
    hyd = pdst.hyd(ixM);      %hydraulic depth
    aoh = amp/hyd;            %tidal amplitude/hydraulic depth   
    Wm = pdst.Wmt(ixM);       %width at mouth 
    Am = pdst.CSAmt(ixM);     %csa at mouth
    Dm = pdst.Dmt(ixM);       %deepest depth in mouth cross-section
    PrvAm = PrV/Am;           %Prism to CSA ratio
    SflShw = (Shw-Slw)/Shw;   %ratio of intertidal area to basin area
    
    if isfield(grid,'Rv') && ~isempty(grid.Rv)
        %if river has been defined for grid use these values
        Wr = grid.Rv.Wr;
        Ar = grid.Rv.Ar;
    else
        %use the values at the tidal limit Lt
        Wr = pdst.Whw(ixH);       %width at head
        Ar = pdst.CSAhw(ixH);     %csa at head    
    end
    if isempty(Wr), Wr = 0; Ar = 0; end

    %convergence length from mouth accounting for any offset
    ewidth = pdst.Wmt-Wr; ewidth(ewidth<0) = 0;
    Lw = -getconvergencelength(xi(ich),ewidth(ich)); %width convergence length at mean tide level
    ecsa = pdst.CSAmt-Ar; ecsa(ecsa<0) = 0;
    La = -getconvergencelength(xi(ich),ecsa(ich)); %csa convergence length at mean tide level

    %check width and CSA at mouth
    % z0 = wl.zmt;              %mean tide level(m)
    % [Wo,Ao] = getCSAat_z0_X(grid,z0,ixM); %width and csa at mouth (cross-check only)
    % fprintf(sprintf('Width: Wm %g, Wo %g;\n CSA: Am %g, Ao %g\n',Wm,Wo,Am,Ao))
    
    %assign output to table
    grossprops = table(Shw,Slw,Vhw,Vlw,PrA,PrV,Gamma,Vs,Vc,Wm,Am,Dm,...
                              Wr,Ar,amp,hyd,aoh,VsoVc,PrvAm,SflShw,Lw,La);                          
    gdsp = setDSproperties();
    gdst = dstable(grossprops,'RowNames',grid.t,'DSproperties',gdsp);                                 
end
%%
function [w,csa] = getCSAat_z0_X(grid,z0,ixM) %#ok<DEFNU>
    %get thw width and cross-sectional area at distance x and elevation z0
    %assumes constant z0 defines a water level surface
    if ~isscalar(z0)
        z0 = z0(ixM); %x increases from mouth - use value at mouth
    end
    yi = grid.y;
    dely = abs(yi(2)-yi(1));  %grid interval
    x = grid.x;
    idx = [1,length(x)];
    %find the width and CSA at both ends of grid and use the maximum values
    W = zeros(1,2); CSA = W;
    for i=1:2
        zi = grid.z(idx(i),:);
        zi(zi>z0) = NaN;
        W(i) = sum(~isnan(zi)*dely);
        dzwl = z0 - zi;
        idd = ~isnan(dzwl);           %index of valid depths
        if sum(idd)<=2
            CSA(i) = 0;
        else
            CSA(i) = trapz(yi(idd),dzwl(idd)); %csa below zwl
        end
    end
    w = max(W);
    csa = max(CSA);
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
                'Name',{'Shw','Slw','Vhw','Vlw','PrA','PrV',...
                         'Gamma','Vs','Vc','Wm','Am','Dm','Wr',...
                         'Ar','amp','hyd','aoh',...
                         'VsoVc','PrvAm','SflShw','Lw','La'},...                                                   
                'Description',{'HW Surface Area','LW Surface Area',...
                         'HW Volume','LW Volume',...
                         'Tidal Prism using CSA','Tidal Prism using hypsomety',...
                         'Gamma','Storage Volume','Channel Volume',...
                         'Width at mouth to MTL','CSA at mouth to MTL',...
                         'Depth at mouth to MTL',...
                         'Width at head','CSA at head',...                                 
                         'Tidal amplitude at mouth','Hydraulic depth to MTL',...
                         'Amplitude to Depth ratio',...
                         'Storage to Channel Volume ratio',...
                         'Prism to CSA ratio at mouth',...
                         'Intertidal to Basin Area ratio',...
                         'Mean tide Width convergence length',...
                         'Mean tide CSA convergence length'},...
                'Unit',{'m2','m2','m3','m3','m3','m3','-','m3','m3',...
                        'm','m2','m','m','m2','m','m','-','-','m','-','m','m'},...
                'Label',{'Surface area (m^2)','Surface area (m^2)',...
                         'Volume (m^3)','Volume (m^3)',...
                         'Prism (m^3)','Prism (m^3)',...
                         'Gamma','Volume (m^3)','Volume (m^3)',...
                         'Width (m)','CSA (m^2)','Depth (m)',...
                         'Width (m)','CSA (m^2)',... 
                         'Tidal amplitude (m)','Hydraulic depth (m)',...
                         'Amplitude to Depth ratio (-)',...                                 
                         'Storage to Channel Volume ratio (-)',...
                         'Prism to CSA ratio (m)',...
                         'Intertidal to Basin Area ratio (-)',...
                         'Convergence lengh (m)','Convergence lengh (m)'},...                                 
                'QCflag',repmat({'model'},[1,22]));  
            dsp.Row = struct(...
                            'Name',{'Time'},...
                            'Description',{'Time'},...
                            'Unit',{'y'},...
                            'Label',{'Time (yr)'},...
                            'Format',{'y'});       
            dsp.Dimensions = struct('Name',{''},'Description',{''},...
                                    'Unit',{''},'Label',{''},'Format',{''});                    
end
