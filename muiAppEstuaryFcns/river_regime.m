function [hrv,wrv,etaw] = river_regime(q,s,d,tau,ros,row)
%
%-------function help------------------------------------------------------
% NAME
%   river_regime.m 
% PURPOSE
%   To compute the width and hydraulic depth of a river section given the
%   discharge and energy slope
% USAGE
%   [hrv,wrv,etaw] = river_regime(q,s,d,tau,ros,row)
% INPUTS
%   q = river discharge (m3/s)
%   s = energy slope (-)
%   d = sediment grain size (m)
%   tau = critical shear stress (Pa)
%   ros = density of sediment (kg/m3)
%   row = density of water (kg/m3)
%   [the above parameters can also be passed as a struct with fields
%    Qr, S, d50, tuacr, rhos, rhow, visc, g]
% RESULTS
%   hrv = hydraulic depth (m)
%   wrv = section width (m)
%   etaw= width regime coefficient
% NOTES
%   see Cao & Knight, 1996, Regime theory of alluvial channels based on the
%   concept of stream power and probability, Proc ICE Water Maritime and
%   Energy, 118, 3, 160-176.
% SEE ALSO
%   Used in CKFA model
%
% Author: Ian Townend
% CoastalSEA (c)June 2016
%--------------------------------------------------------------------------
%
    if nargin==1 && isstruct(q)
        [q,s,d,tau,ros,row,neu,g] = getInputs(q);
    else
        neu = 1.34*10^-6; %kinematic viscosity
        g   = 9.81;       %acceleration due to gravity        
    end

    % capture no flow case
    if q==0 
        hrv=0; wrv=0; etaw=0; 
        return 
    end
    %
    hrv = 1; %initialise depth
    urv = 1; %initialise flow velocity
    %
    sed = ros/row;
    dgr = d*(g*(sed-1)/neu^2)^(1/3);   %dimensionless grain size (eq.54)
    if dgr>60                       
        nck = 0;
        Ack = 0.17;
    elseif dgr>1
        nck = 1-0.56*log(dgr);         %eq.56
        Ack = 0.14+0.23/sqrt(dgr);     %eq.57
    else %60 microns is the lower bound for the equations in the paper
        nck = 1; %no basis for this just avoids divide by zero in urv
        Ack = 0.14+0.23/sqrt(dgr); 
    end
    dh  = 1;
    count = 0;
    %
    %natural logarithms used rather than log10 with adjustments accordingly
    %(see CK river regime.xmcd for transformation)
    while dh>0.001
        ust = sqrt(g*hrv*s);                        %shear velocity
        eta = 3*(g*row*q*s/tau/(hrv^1.5)/urv)^0.4;  %regime coefficient (eq.50)
        Ffg = ust*(g*d*(sed-1))^-0.5;               %total shear (eq.52)
        if dgr>=1
            denom= exp(0.242*(log(dgr))^1.7);        %part of eq.51
        else
            denom = 1;
        end
        Fgr = 0.04*(6*denom*Ffg+19*denom*Ack+19*Ffg-19*Ack)/denom; %from eq.51
        urv = 2*log(1/Fgr)-log(g*d*(sed-1))+2*nck*log(ust);        %from 53?
        urv = 2.457*exp(0.5*urv/(nck-1))*log(10*hrv/d);
        fri = 8*(ust/urv)^2;                        %D-W friction factor eq.55
        etah= (fri/8/g/eta^3)^(1/6);                %eq.19
        hrc = etah*q^(1/3)*s^(-1/6);                %eq.17
        count = count+1;
        if count > 100
            sprintf('Warning: Exceeded 100 iterations in river_regime.m')
            return
        end
        dh  = abs(hrc-hrv);
        hrv = hrc;
        if hrv<=0, wrv=0; etaw = []; return; end
    end
    etaw = (fri*eta^3/8/g)^0.25;                     %eq.18
    wrv  = etaw*q^0.5*s^-0.25;                       %eq.17
end
%%
function [q,s,d,tau,ros,row,neu,g] = getInputs(inp)
    %unpack input parameters from struct q
    q = inp.Qr;
    s = inp.Sr;
    d = inp.d50;
    tau = inp.taucr;
    ros = inp.rhos;
    row = inp.rhow;    
    neu = inp.visc;
    g = inp.g;
end
