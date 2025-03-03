%% PLinterface
% Abstract class with properties and methods to manipulate point and line 
% datasets. The class inherits the handle class.
%% Description
% *PLinterface* is used as a superclass to provide data handling 
% functionality for classes that interactively manipulate points and lines.
% Classes that inherit *PLinterface* define menus as a subset of the options 
% defined in the _setFigure_ function, or overload this
% function if a more customised set of actions is required.
% Any class using the interface must define the Abstract properties for
% outPoints, outLines and isXYZ.

%% PLinterface properties
% The properties in the abstract class include the following:
%%
% * *Figure* - handle to class figure
% * *Axes* - handle to class axes
% * *Points* - current set of points
% * *pLines* - current set of lines
% * *NaNpnts* - x,y,(z) NaN points used to terminate lines

%%
% Abstract properties to to be defined by subclasses:
%%
% * *outPoints* - saved points for output
% * *outLines* - saved lines for output
% * *isXYZ* - logical flag, true if points are lines are to include z values

%% PLinterface methods
% The methods in the abstract class include the following:
%%
% Methods to add, delete and update grids:
%%
% *setFigure*
% - initialise a figure and panel for PL classes
%%
%   obj = setFigure(obj,figtitle,tag,position);
%%
% where _obj_ is an instance of any class that inherits *PLinterface*,
% _figtitle_ is text for the figure title, _tag_ is text for the figure Tag 
% and _position_ is the figure Position (optional with default of   
% [0.37,0.58,0.26,0.34].

%%
% Static methods:
%%
% *checkDirection*
% - check direction of new points relative to existing line and flip if
% necessary so that all points are in the same directio
%%
%   newpnts = obj.checkDirection(endpnt,newpnts,isstart);
%%
% where obj is an instance of any class that inherits from PLinterface, 
% _endpnt_ is the end point of a line, _newpnts_ are the points to be
% added to the line and _isstart_ is true if the end point is the first
% point in the line being extended.

%%
% *isPointNearLine*
% - calculate the distance from the point to points on line and
% also return sorted indices and distances to all Points
%%
%   [isNear,idP,distances] = obj.isPointNearLine(points,point,tol);
%%
% where obj is an instance of any class that inherits from PLinterface,
% _points_ defines the line, _point_  the point to be tested and _tol_ the
% tolerance to be allowed in checking proximity. Returns _isNear_ if the
% point is within tolerance along with _idP_ and _distances_ ranking the
% indices and distances of all _points_ from _point_.

%%
% *lineLength*
% - find length of line between 2 points
%%
%    distance = obj.lineLength(pnt1,pnt2);
%%
% where obj is an instance of any class that inherits from PLinterface,
% _pnt1_ and _pnt2_ define 2 points and the function returns the _distance_
% between then.

%%
% *setInterval*
% - prompt user to set point spacing interval for the contour sampling
% along a line.
%%
%   cint = obj.setInterval();

%%
% *setLevel*
% - prompt user to set the level for the contour to be extracted and
% returns the level.
%%
%   zlevel = obj.setLevel();

%%
% *setSmoothinInputs*
% - prompt user for the smoothing input parameters to define the method, 
% window size, degree and mimimum number of points to smooth.
%%
%   inp = obj.setSmoothingInputs();

%% See Also
% Some classes that use *PLinterface* and related functions are detailed in 
% the help for <matlab:doc('grid_point_line') Points and Lines>. These work
% in conjuction with the <matlab:doc('gdinterface') GDinterface> and 
% functions detailed in the help for <matlab:doc('grid_class_fcns') Grid classes and functions>.
