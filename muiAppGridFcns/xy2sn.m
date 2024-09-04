function [s_out,n_out]=xy2sn(x_c_in,y_c_in,xi,yi,varargin)
% XY2SN Transforms Cartesian coordinates to channel fitted coordinates
%
%   [Si,Ni]=XY2SN(Xc,Yc,Xi,Yi) Converts the Cartesian coordinates Xi and Yi 
%   to the channel fitted coordinate Si and Ni using the centreline defined 
%   through the Cartesian coordinates of its vertices Xc and Yc. 
%
%   [Si,Ni]=XY2SN(...,'ParamName',Value) Allows to set the following
%   options:
%
%   'MaximumIterations'
%   Scalar number indicating the maximum number of iterations performed by
%   the algorithm, default is 500.
%
%   'Tolerance'
%   Scalar number indicating the maximum tolerance in the estimated
%   coordinates. Default is 1/1000 of the mean sampling distance of the
%   centerline
%
%   See also: sn2xy amd pp_cline
%
%   Bart Vermeulen (2022). Cartesian to Curvilinear coordinate forward and backward transformation
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation
%
% MODIFICATIONS
%   Source code assumes a square grid: isequal(size(xi),size(yi))
%   This is not a requirement of delauneyTriangulation function so modified
%   to use x and y (IHT, Jul 2022)
%   Source code returns a grid that is flipped in N relative to the input
%   channel coordinates. Reversed sign of Nfit array in line 110 to correct

%% Hanldle input
assert(isequal(class(x_c_in),class(y_c_in)) && isequal(size(x_c_in),size(y_c_in)) && isnumeric(x_c_in))
% assert(isequal(class(xi),class(yi)) && isequal(size(xi),size(yi)) && isnumeric(xi))
assert(isequal(class(xi),class(yi)) && isnumeric(xi))
iP=inputParser;
iP.FunctionName='xy2sn';
iP.addParameter('Tolerance',1e-3,@(x) isscalar(x) && isnumeric(x) && x>0)
iP.addParameter('MaximumIterations',500,@(x) isscalar(x) && isnumeric(x) && x>=0 && mod(x,1)==0);
iP.parse(varargin{:});

size_out=size(xi);
if isrow(x_c_in), x_c_in = x_c_in'; end  %force a column vector 
if isrow(y_c_in), y_c_in = y_c_in'; end  %force a column vector
P=[x_c_in, y_c_in]';                     %2 x row vector

if isrow(xi), xi = xi'; end              %force a column vector 
if isrow(yi), yi = yi'; end              %force a column vector 
Pi=[xi, yi]';
%% Parameters
max_iter=iP.Results.MaximumIterations;
dt=sqrt(sum(diff(P,1,1).^2,2));
if strcmp('Tolerance',iP.UsingDefaults)
    ds_tolerance=1e-3*nanmean(dt);
else
    ds_tolerance=iP.Results.Tolerance;
end
clear iP

%% Preprocess centreline
ppc=pp_cline(P);
ppcd=fnder(ppc);
s=ppc.breaks;
P=ppval(ppc,s);

%% Find nearest centreline points 
DT=delaunayTriangulation(P');
try
    fc=DT.nearestNeighbor(Pi');
catch err
    if strcmp(err.identifier,'MATLAB:delaunayTriangulation:TriIsEmptyErrId')
        warning('xy2sn_legl:AddingJitter','Adding random jitter to handle collinearity')
        P=P+rand(size(P))*ds_tolerance^2;
        DT.Points=P';
        fc=DT.nearestNeighbor(Pi');
    else
        throw(err)
    end
end
T=ppval(ppcd,s);
P_fit=P(:,fc);
Tfit=T(:,fc);
citer=1;
sfit=s(fc);
last_ds=inf;
count_inc=0;
relax_iter=1;
hw = waitbar(0,'Processing conformal mapping');
while citer<=max_iter && last_ds>ds_tolerance
    PC=Pi-P_fit;
    deltas=dot(PC,Tfit,1)./sqrt(sum(Tfit.^2,1)); %Calculate projection of points onto vector along centreline (dot product)
    maxds=nanmax(deltas);
%     disp(num2str(maxds));
    if maxds>last_ds, count_inc=count_inc+1; else count_inc=max(0,count_inc-.1); end
    if count_inc>1
        relax_iter=max(0.1,0.5*relax_iter);
%         disp(['relax: ',num2str(relax_iter)]);
        count_inc=0;
    end
    last_ds=maxds;
    sfit=sfit+relax_iter*deltas;
    P_fit=ppval(ppc,sfit);
    Tfit=ppval(ppcd,sfit);
    citer=citer+1;
    waitbar(citer/max_iter,hw)
end
delete(hw)
citer=citer-1;
if citer>=max_iter
    warning(['Maximum iterations reached, maximum delta_s: ',num2str(maxds)]);
end

Nfit=[Tfit(2,:); -Tfit(1,:)];              %signs reversed, iht, Aug 2022
nfit=dot(PC,Nfit,1)./sqrt(sum(Nfit.^2,1)); %Calculate projection of points onto vector along centreline (dot product)

s_out=reshape(sfit,size_out);
n_out=reshape(nfit,size_out);

end
