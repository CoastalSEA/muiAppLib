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

assert(isequal(class(x_c_in),class(y_c_in)) && isequal(size(x_c_in),size(y_c_in)) && isnumeric(x_c_in) && all(isfinite(x_c_in)) && all(isfinite(y_c_in)))
assert(isequal(class(si),class(ni)) && isequal(size(si),size(ni)) && isnumeric(si))

size_out=size(si);
P=[x_c_in(:) y_c_in(:)]';
Pi=[si(:) ni(:)]';

%% Preprocess centreline
ppc=pp_cline(P);
ppcd=fnder(ppc);

%% compute x and y coordinates
Tfit=ppval(ppcd,Pi(1,:));
Nfit=[-Tfit(2,:); Tfit(1,:)];
Nfit=bsxfun(@rdivide,Nfit,sqrt(sum(Nfit.^2,1)));
P_out=ppval(ppc,Pi(1,:))+bsxfun(@times,Pi(2,:),Nfit);
x_out=reshape(P_out(1,:),size_out);
y_out=reshape(P_out(2,:),size_out);
