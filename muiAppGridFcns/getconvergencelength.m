function len = getconvergencelength(xdata,ydata,x0)
%
%-------function help------------------------------------------------------
% NAME
%   getconvergencelength.m
% PURPOSE
%   least squares fit using fminsearch (curve fitting via optimization) to
%   find the convergence length of a channel from x-y data 
% USAGE
%   len = getconvergencelength(xdata,ydata)
% INPUTS
%   xdata - independent variable, e.g. x-coordinates for ydata
%   ydata - dependent variable to be fitted with best fit exponential
%   x0 - initial guess (see fminsearch documentation for details)
% OUTPUT
%   len - convergence length of ydata
% SEE ALSO
%   used in cf_property_plots and gd_grossproperties
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    if nargin<3
        x0 = [ydata(1),1/(max(xdata)/4)]; %initial guess needs to be of right order of magnitude
    end
    
    if all(size(xdata)~=size(ydata))
        ydata = reshape(ydata,size(xdata));
    end
    fun = @(x)sseval(x,xdata,ydata);
    
    options = optimset('MaxIter',5000,'MaxFunEvals',5000,'TolFun',1e-3,'TolX',1e-3);
    bestx = fminsearch(fun,x0,options);
    lambda = bestx(2);
    len = -1/lambda;
    function sse = sseval(x,xdata,ydata)
        A = x(1);
        alambda = x(2);
        sse = sum((ydata - A*exp(-alambda*xdata)).^2);
    end
end