function ppc=pp_cline(P)
% PP_CLINE make splines of centreline parameterized to centreline length
%
%   PPC = PP_CLINE(X_C,Y_C) computes two cubic piecewise polynomials
%   representing the x and y coordinates of the centreline parameterized in
%   in the centreline arc length. The centreline is given throgh vertex
%   coordinates X_C and Y_C. This function is a helper function for XY2SN
%   and SN2XY.

assert(isnumeric(P) && ismatrix(P) && size(P,1)==2)


%% Parametric description with spline
t=cumsum([0 sqrt(sum(diff(P,1,2).^2,1))]);
ppc=spline(t,P);

%% Arc length calculation
s=cumsum([0 arrayfun(@(a,b) integral(@(t) sqrt(sum(ppval(fnder(ppc),t).^2,1)),a,b),t(1:end-1),t(2:end))]);

%% Create spline with s (real arc length) as parameter
ppc=fnxtr(spline(s,ppval(ppc,t))); 

%% Resample
s=linspace(0,s(end),numel(s));
ppc=fnxtr(spline(s,ppval(ppc,s)));

end
