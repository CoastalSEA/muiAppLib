function [Hsi,Diri,depS,bS] = hs_surf(dst,inp)
%
%-------function help------------------------------------------------------
% NAME
%   hs_surf.m
% PURPOSE
%   Calculate the inshore wave height at the edge of the surf zone
% USAGE
%   [Hsi,Diri,depS,bS] = hs_surf(inp)
% INPUTS 
%   dst - dstable with variables:
%       Hs  - offshoresignificant wave height (m)
%       Tp   - peak wave period (s)
%       Dir - wave direction (degrees TN)
%       swl  - still water level (mOD)
%   inp - struct derived from ctWaveParameters instance with properties:
%       z0   - offshore bed level (mOD)
%       zi   - inshore bed level (mOD) [if used instead of surf zone depth]
%       offtheta - angle of contours at offshore wave data location (degrees TN)
%       intheta- angle of contours at shoreline (degrees TN) 
%       Kf   - friction coefficient (default=1)
%       zBC  - beach crest level (mOD)
%       ubs  - upper beach slope (1:bs)
%       z1km - bed level 1km from shore (mOD)
%       hboption - flag for wave breaking model in hb_break,1=SPM breaking on a slope
%       hsbflag - 0=no breaking; 1=SPM Hb; 2=La Roux Hb; 3=SPM+Holmes/TCP Hsb used in hs_surf.getHsb
%       g    - acceleration due to gravity (m/s-2)
% OUTPUTS
%   Hsi  - inshore wave height at edge of surf zone (m)
%   Diri - wave direction at edge of surf zone (degrees TN)
%   depS - depth at edge of surf zone, or the defined inshore point (m)
%   bS   - slope at the edge of surf zone, or the defined inshore point (1:bS)
% NOTES
%   Iterates using the depth at which the depth limited breaker height and 
%   the refracted wave height converge. See Le Roux J P, 2007, A simple 
%   method to determine breaker height and depth for different deepwater 
%   wave height/length ratios and sea floor slopes. Coastal Engineering
% SEE ALSO
%   refraction.m, shoaling.m, celerity.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%----------------------------------------------------------------------
%

    %if there are no wave directions assume normal to shore and prompt user
    if isempty(dst.Dir) || all(isnan(dst.Dir))                
        %get input parameters from user
        prompt = {'Mean wave direction (degTN):'};
        title = 'Define wave direction';
        numlines = 1;
        defaultvalues{1} = num2str(inp.intheta+90);
        useInp=inputdlg(prompt,title,numlines,defaultvalues);
        if isempty(useInp)       %user cancelled
            Hsi=[];Diri=[];depS=[];
            return; 
        end 
        dst.Dir = useInp*ones(size(dst.Hs));
    end 
    
    %initialise output variables
%     Hsi=zeros(size(inp.Hs)); Diri=Hsi; depS=Hsi;
    
    dep0  = dst.swl-inp.z0;  %water depth from swl to bed level 
    
    %check plots for beach profile used to estimate slopes                                 
    %[ybs,zbs,mbs] = beachprofile(inp); %get the beach profile and slopes 
    %plotprofile(ybs,zbs,mbs);

    if isnan(inp.zi)
        [depS,bS] = surfdepth(dst,inp,dep0);
    else
        depS = dst.swl-inp.zi; %water depth from swl to bed level        
        bS = profileslope(-inp.zi,0,inp.z1km,inp.ubs);
    end
    bS(depS<0) = NaN;   %depS can be negative if fzero does not find solution
    depS(depS<0) = NaN; %remove any negative water depths and related slopes   
    %call refraction model (uses celerity) to get inshore values of
    %a wave height and direction at the edge of the surf zone or at zi 
    %this is a beach/nearshore model so isshore is true
    [Hsi,Diri] = refraction(dst.Hs,dst.Tp,dst.Dir,[dep0,depS],...
                                  [inp.offtheta,inp.intheta],inp.Kf,true);
    if ~isnan(inp.zi)             
        bsi = profileslope(-inp.zi,0,inp.z1km,inp.ubs);
        Hsb = getHsb(Hsi,dst.Tp,bsi,depS,1,inp.g,inp.hboption,inp.hsbflag);
        Hsi = min(Hsi,Hsb,'includenan');%Hs wave height in depth depi
    end                                                            
end
%%
function Hsb = getHsb(Hsi,Tp,bs,depi,beta,g,hboption,hsbflag)
    %get the breaking wave height for a given depth, using one of three
    %hboption: 0=0.78Hsi; 1=SPM on a slope; 2=SPM in front of toe; 3= as 2 for overtopping
    %hsbflag:  1=SPM Hb; 2=La Roux Hb; 3=SPM+Holmes/TCP Hsb    
    switch hsbflag        
        case 0  %no wave breaking included
            Hsb = Hsi;
        case 1  %SPM formulation
            Hsb = hb_break(depi,bs,Tp,g,hboption);
        case 2  %la Roux formulation, 2007, eq.25
            alp = atan(1/bs)*180/pi;
            Hsb = depi*(-0.0036*alp^2+0.0843*alp+0.835);
        case 3  %SPM to define maximum and Holmes/TCP to get Hsb
            [hsb1,hsb2] = hs_break(Hsi,Tp,bs,depi,beta,g,hboption);
            %assign Hsb based on bed slope
            bs1 = 100;
            bs2 = 20;
            %transition values are estimates
            if bs>=bs1
                Hsb = hsb1; %Holmes - energy redistributed (lower bound)
            elseif bs<=bs2
                Hsb = hsb2; %Tucker,Carr,Pit-energy at Hb (upper bound)
            else
                Hsb = hsb1+(bs1-bs)*(hsb2-hsb1)/(bs1-bs2); %transition
            end
    end
end
%%
function [depS, mS] = surfdepth(dst,inp,dep0)
    %iterate to find the water depth at the edge of the surf zone where
    %breaking commences using fzero. return depth and bed slope
    if isa(dst,'dstable')
        dst = dst.DataTable; %convert to a table for use in parfor
    end

    depS = zeros(size(dst.Hs)); mS = depS;
    %index for nans in input data
    idx = isnan(dst.Hs) | isnan(dst.Tp) | isnan(dst.Dir) | isnan(dep0);
    %additional input values
    beta = 1;        %default breaking coefficient
    
    %initial guess of inshore depth
    depi1 = 1.3*dst.Hs; 
    %define variables to control root search (no longer used)
    Hbar = mean(dst.Hs,1,'omitnan'); %mean wave height
    Hstd = std(dst.Hs,1,'omitnan');  %standard deviation of wave height
    depcl = 2*Hbar+11*Hstd;          %profile closure depth: assume littoral 
                                     %drift takes place inside this depth    
    depi2 = ones(size(dst.Hs))*depcl;
    %if initial guess of inshore depth is greater than offshore depth
    %limit to Hsb in offshore depth 
    bsi = profileslope(depcl,0,inp.z1km,inp.ubs);
    idd = depi1>depcl;
    [~,dst.Hs(idd)] = hs_break(dst.Hs(idd),dst.Tp(idd),bsi,depcl,beta,inp.g,inp.hboption);
    depi1(idd) = depi2(idd);
    depi2(idd) = 1.3*dst.Hs(idd);

    %inshore wave and breaker height anonymous functions
    mbsfun = @(depi,j) profileslope(depi,dst.swl(j),inp.z1km,inp.ubs);
    hsifun = @(depi,j) refraction(dst.Hs(j),dst.Tp(j),dst.Dir(j),...
                    [dep0(j),depi],[inp.offtheta,inp.intheta],inp.Kf,true);
    hmxfun = @(depi,j) hsifun(depi,j)*1.27;  %highest 10%
    hsbfun = @(depi,j) getHsb(hmxfun(depi,j),dst.Tp(j),mbsfun(depi,j),...
                                      depi,beta,inp.g,inp.hboption,inp.hsbflag);                                                                                                                     
    myfunc = @(depi,j) hmxfun(depi,j)-hsbfun(depi,j);
    options = optimset('Display','off');   
    options.FunValCheck = 'on';
    nrec = length(dst.Hs);
    hpw = PoolWaitbar(nrec, 'Processing wave data');
    parfor j=1:nrec
        jfunc = @(depi) myfunc(depi,j);
        if depi1(j)>0.2 && depi2(j)>depi1(j) && ~(idx(j))
            depS(j) = findepth(jfunc,depi1(j),depi2(j),options);            
            if isnan(depS(j)) || depS(j)<1 %depS can be negative if fzero does not find solution
                mS(j) = NaN;
            else
                mS(j) = mbsfun(depS(j),j);
            end
        elseif idx(j)
            depS(j) = NaN;
            mS(j) = NaN;
        elseif depi2(j)<depi1(j)
            mS(j)=1;
        else
            depS(j) = 0.2; %default for small waves, Hs<0.154m
            mS(j) = mbsfun(depS(j),j);
        end
        increment(hpw);
    end
    delete(hpw)
end
%%
function depS = findepth(afunc,depi1,depi2,options)
    %use depi1 to find root. If this fails, iterate out to closure depth        
    [depS,~,exitflag] = fzero(afunc,depi1,options);
    if exitflag<1 || ~isreal(depS)
        nint = 5;
        dint = abs(depi2)/nint;
        for ii=1:nint
            depii = depi1+dint*ii;
            [depS,~,exitflag] = fzero(afunc,depii,options);
            if exitflag>0
                break
            end
        end
        if exitflag<1
                depS = exitflag; %errors are -1 to -6
        end
    end
end
