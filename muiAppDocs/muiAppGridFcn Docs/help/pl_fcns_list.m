%% Point and Line Classes and Functions
% Functions used to manipulate points and lines can be found in the
% _muiAppGridFcns_ folder and include the following:

%% Classes
% * *PL_Sections*: provides methods to manipulate points, lines and sections 
% using datasets in classes that inherit <matlab:doc('gdinterface') GDinterface>.
% * *PL_Boundary*: Class to extract contours and generate model boundaries
% * *PL_CentreLine*: Class to extract valley/channel centre-line
% *PL_PlotSections*: display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
% * *PL_SectionLines*: Class to create a set of cross-sections lines based on 
% the points in a centre-line and any enclosing boundary.

%% Functions
% * *gd_boundary.m* - UI to interactively generate a contour line (e.g. to define a boundary).
% * *gd_centreline* - create a centreline of a channel to trace the
% shortest path between start and end points whilst finding the deepest
% points (i.e. a thalweg).
% * *gd_cplines2plines.m*
% - convert cell array of plines to an array of points that define plines.
% * *gd_curvelineprops* - for each point from idL to the end use the centre-line coordinates and 
% direction to find the lengths and directions along the centre-line.
% * *gd_findline.m* - find which _pline_ a given _point_ lies on when there are multiple lines
% separated by NaN values in a _plines_ struct array.
% * *gd_getcontour.m* - extract a contour at a defined level.
% * *gd_getpline.m* - interactively select a line on a plot and return the line point coordinates.
% * *gd_getpoint.m* - interactively select a single point in a UI figure and return the point
% coordinates.
% * *gd_lines2points.m* - convert _lines_  as x,y vectors in various formats (see *gd_points2lines*) to a 
% _points_ array of structs for each _point_ with x, y (and z) fields.
% * *gd_linetopology* - interactively define line connectivity.
% * *gd_orderplines.m* - amend the order of the _plines_ in an array of _plines_.
% * *gd_plines2cplines.m*
% - convert an array of _plines_ to a cell array of _plines_.
% * *gd_plotgrid* - create pcolor plot of gridded surface (from Grid Tools).
% * *gd_plotsections* - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure (from Grid Tools).
% * *gd_points2lines.m* - convert a _points_ or plines array of structs with x, y (and z) fields 
% to a [Nx2] or [Nx3] array, or a single stuct, or matrix with column vectors
%  in the x, y (and z) fields.
% * *gd_sectionlines.m* - extract the section lines that are normal to the channel centre line
% and extend to the bounding shoreline.
% * *gd_selectpoints.m* - UI to interactively create a specified number of x,y point on a grid.
% * *gd_setpoint.m* - interactively create a single point in a UI figure and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. elevation).
% * *gd_setpoints.m*
% - interactively create a set of points in a UI figure and return the point
% coordinates. Includes an option to enter an additional value at the
% selected points (e.g. elevation).
% * *gd_smoothlines.m* - smooth one or more line segments using either moving average, or the
% Savitzky-Golay method.
