function ppc=pp_cline(P)
% PP_CLINE make splines of centreline parameterized to centreline length
%
%   PPC = PP_CLINE(X_C,Y_C) computes two cubic piecewise polynomials
%   representing the x and y coordinates of the centreline parameterized in
%   in the centreline arc length. The centreline is given throgh vertex
%   coordinates X_C and Y_C. This function is a helper function for XY2SN
%   and SN2XY.
%
%   See also: sn2xy and pp_cline
%
%   Bart Vermeulen (2022). Cartesian to Curvilinear coordinate forward and backward transformation
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation
%
% MODIFICATIONS
%   Oct 2024- added check for curve fitting toolbox and alternative to use ppdiff from splinefit.  

    assert(isnumeric(P) && ismatrix(P) && size(P,1)==2)
    isok = license('test','Curve_Fitting_Toolbox');  %toolbox is licensed to use
    if isok
        addons = matlab.addons.installedAddons;
        isok = any(ismatch(addons.Name,'Curve Fitting Toolbox')); %toolbox is installed
    end

    % Parametric description with spline
    t=cumsum([0 sqrt(sum(diff(P,1,2).^2,1))]);
    ppc=spline(t,P);

    if isok         %requires Curve Fitting Toolbox
        % Arc length calculation
        s=cumsum([0 arrayfun(@(a,b) integral(@(t) sqrt(sum(ppval(fnder(ppc),t).^2,1)),a,b),t(1:end-1),t(2:end))]);

        % Create spline with s (real arc length) as parameter
        ppc=fnxtr(spline(s,ppval(ppc,t))); %requires Curve Fitting Toolbox and uses default order of 2

        % Resample
        s=linspace(0,s(end),numel(s));
        ppc=fnxtr(spline(s,ppval(ppc,s)));%requires Curve Fitting Toolbox

    else            %requires splinefit functions from Matlab Forum        
        % Arc length calculation
        s=cumsum([0 arrayfun(@(a,b) integral(@(t) sqrt(sum(ppval(ppdiff(ppc),t).^2,1)),a,b),t(1:end-1),t(2:end))]);

%         %NEED TO FIND ALTERNATIVE TO FNXTR ********************************  
%         %This runs but does not give a new surface*************************
%         % Create spline with s (real arc length) as parameter      
%         ppv = spline(s,ppval(ppc,t));
%         % Resample
%         sq=linspace(0,s(end),numel(s));
%         
%         vq=interp1(sq,ppval(ppv,sq),'spline','extrap'); 
%         ppc = spline(sq,vq); 

        % Create spline with s (real arc length) as parameter
        ppc=fnxtr(spline(s,ppval(ppc,t))); %requires Curve Fitting Toolbox and uses default order of 2
        % Resample
        s=linspace(0,s(end),numel(s));
        ppc=fnxtr(spline(s,ppval(ppc,s)));%requires Curve Fitting Toolbox

    end
end
