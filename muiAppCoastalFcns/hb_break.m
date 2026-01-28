function Hb = hb_break(d,bs,T,g,hbflag)
%
%-------header-------------------------------------------------------------
% NAME
%   hb_break.m
% PURPOSE
%   Wave height after breaking for given water depth d
% USAGE
%       Hb = hb_break(d,bs,Tp,flag)
% INPUTS
%   d    - water depth at point of interest (m)
%   bs   - bed slope (1:bs)
%   T    - wave period (s)
%   g    - acceleration due to gravity(m/s2)
%   hbflag - choice of algorithm:, 0,1,2,3
%            0: simple ratio of 0.78 (maximum for a solitary wave)
%            1: SPM breaking on a slope
%            2: SPM breaking on a slope some distance in front of toe
%            3: SPM modified for overtopping some distance in front of toe
% OUTPUT
%   Hb   - depth limited breaking wave height
% NOTES
%   For options 1-3 usesSPM equation for depth of breaking wave height 
%   on a slope. Options 2 and 3 are at a sea wall and option 3 is as modified
%   by IHT for use with Owens overtopping model. 
%   To get energy equivalence should use Tp and Hrms (see Soulsby 1997, p69)
%   Exception is for flag=3 where Tz is used in overtopping calculations
% SEE ALSO
%   Functions hs_break.m, otop_Q.m, otop_C.m, hb_break.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%--------------------------------------------------------------------------
%
% routine specific coefficients
    kapa = 4.46;     %coefficient for 'a' in SPM eqn 2-93
    kapb = 1.56;     %coefficient for 'b' in SPM eqn 2-94
    if hbflag==3
        kapb = 1.12; %coefficient for 'b' in SPM eqn 2-94 modifed for overtopping
    end
    kp = 1.0;     %plunge length multiplier for breaking distance from toe
    %
    fa = kapa*g*(1-exp(-19./bs));
    fb = kapb*(1+exp(-19.5./bs)).^-1;
    fc = kp*(4-9.25./bs)./bs;   %coefficient for plunge length/bs (SPM eqn 7-4) 
    %
    if hbflag==0
        Hb = 0.78*d;
    elseif hbflag==1
        Hb = fb.*d./(1+fa.*d./(g*T.^2));
    elseif hbflag==2 || hbflag==3
        Hb = hb_toe(d,T);
    else
        warndlg('Error in hs_break: incorrect assignment of flag');
    end

%
%--------------------------------------------------------------------------
% nested function to calculate breaker height a plunge length in front of
% toe where the depth, d, is the depth at the toe
%--------------------------------------------------------------------------
%
    function hb = hb_toe(d,T)
        % function to estimate breaker height based on SPM eqn 2-92
        % at a depth a plunge length in front of toe with toe depth, dd 
        fun_h = @(hb,dd,TT,Fa,Fb,Fc) hb-(dd+Fc.*hb).*(Fb-Fa*hb./(g*TT.^2));
        %
        % function to find root of f(dw)
        options = optimset('TolX',1e-3);
        hb = NaN(size(d));
        for i = 1:numel(d)
        
            % Skip invalid inputs early
            if isnan(d(i)) || isnan(T(i)) || d(i) <= 0
                hb(i) = 0;
                continue
            end
        
            % Select scalar parameters (vector or scalar inputs allowed)
            Fa = fa(min(i,end));
            Fb = fb(min(i,end));
            Fc = fc(min(i,end));
            TT = T(min(i,end));
        
            % Define function of x only
            fun = @(x) fun_h(x, d(i), TT, Fa, Fb, Fc);
        
            % Bracket: [0, 5*d] is physically meaningful
            a = 0;
            b = 5*d(i);
        
            % Check function values before calling fzero
            fa_val = fun(a);
            fb_val = fun(b);
        
            if ~isfinite(fa_val) || ~isfinite(fb_val)
                hb(i) = NaN;
                continue
            end
        
            % Require a sign change for robustness
            if sign(fa_val) == sign(fb_val)
                hb(i) = NaN;
                continue
            end
        
            % Solve
            [xval,~,exitflag] = fzero(fun, [a b], options);
        
            if exitflag > 0 && isfinite(xval)
                hb(i) = xval;
            else
                hb(i) = NaN;
            end
        end
    end
end



