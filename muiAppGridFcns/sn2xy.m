function [x_out, y_out]=sn2xy(x_c_in,y_c_in,si,ni)
% SN2XY Transforms channel fitted coordinates to Cartesian coordinates
%
%   [Xi,Yi]=SN2XY(Xc,Yc,Si,Ni) Converts the channel fitted coordinate Si 
%   and Ni to the Cartesian coordinates Xi and Yi using the centreline 
%   defined through the Cartesian coordinates of its vertices Xc and Yc. 
%
%   See also: sn2xy and pp_cline
%
%   Bart Vermeulen (2022). Cartesian to Curvilinear coordinate forward and backward transformation
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation
%
% MODIFICATIONS
%   Source code for xy2sn returns a grid that is flipped in N relative to the input
%   channel coordinates. Reversed sign of Nfit array in line 110 to
%   correct. Made the same change here in line 32 to be consistent.
%   Oct 2024 - added check for curve fitting toolbox and alternative to use ppdiff from splinefit
    
    assert(isequal(class(x_c_in),class(y_c_in)) && isequal(size(x_c_in),size(y_c_in)) && isnumeric(x_c_in) && all(isfinite(x_c_in)) && all(isfinite(y_c_in)))
    assert(isequal(class(si),class(ni)) && isequal(size(si),size(ni)) && isnumeric(si))

    isok = license('test','Curve_Fitting_Toolbox');  %toolbox is licensed to use
    if isok
        addons = matlab.addons.installedAddons;
        isok = any(ismatch(addons.Name,'Curve Fitting Toolbox')); %toolbox is installed
    end

    size_out=size(si);
    if isrow(x_c_in), x_c_in = x_c_in'; end  %force a column vector 
    if isrow(y_c_in), y_c_in = y_c_in'; end  %force a column vector
    P=[x_c_in, y_c_in]';                     %2 x row vector
    
    if isrow(si), si = si'; end              %force a column vector 
    if isrow(ni), ni = ni'; end              %force a column vector 
    Pi=[si, ni]';
    
    
    %% Preprocess centreline
    ppc=pp_cline(P);
    if isok
        ppcd=fnder(ppc); %requires Curve Fitting Toolbox
    else
        ppcd=ppdiff(ppc);  %requires splinefit functions from Matlab Forum
    end
    %% compute x and y coordinates
    Tfit=ppval(ppcd,Pi(1,:));
    Nfit=[Tfit(2,:); -Tfit(1,:)];     %signs reversed, iht, Aug 2022
    Nfit=bsxfun(@rdivide,Nfit,sqrt(sum(Nfit.^2,1)));
    P_out=ppval(ppc,Pi(1,:))+bsxfun(@times,Pi(2,:),Nfit);
    
    x_out=reshape(P_out(1,:),size_out);
    y_out=reshape(P_out(2,:),size_out);
end
