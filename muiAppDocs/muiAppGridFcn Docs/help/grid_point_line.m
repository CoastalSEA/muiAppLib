%% Points, Lines, Sections and Segments
% A short summary of the conventions used and the functions available, to
% manipulate points, lines, sections and segments in Grid Tools.

%% Conventions
% A grid held in any class that inherits <matlab:doc('gdinterface') GDinterface> 
% and obtained from the dstable or using _getGrid_ will return an xyz
% struct where x and y are column vectors and z is a matrix. The 
% <matlab:doc('grid_class_fcns') Grid Tools> function _gd_plotgrid_ can 
% be used to plot the grid. The *GD_Sections* class and associated functions
% make use of the following conventions:
%%
% * *Points*: a single pair of x,y coordinates or a set of points. Use
% _gd_pnt2vec_ and gd_vec2pnt to convert between points and other formats.
%%
% <html>
% <ul><ul>
% <li><i>point</i> - a struct with fields x, y (and z) containing the coordinates of a point.</li>
% <li><i>points</i> - a row struct array of xy(z) <i>point</i> structs .</li>
% </ul></ul>
% </html>

%%
% * *Lines*: a single line with 2 or more points, or a set of lines.
%%
% <html>
% <ul><ul>
% <li><i>line</i> - x, y (and z) column vectors terminated with a Nan. 
% Can be a matrix, struct or table.</li>
% <li><i>lines</i> - x, y (and z) column vectors with line segments defined by NaN terminators.</li>
% </ul></ul>
% </html>

%%
% * *pLines*: a set of 2 or more _points_ that define a line.
%%
% <html>
% <ul><ul>
% <li><i>pline</i> - a <i>points</i> array terminated with a NaN <i>point</i>.</li>
% <li><i>plines</i> - cell array of <i>pline</i> segments concatenated as a row struct array, 
% with individual <i>pline</i> segments defined by NaN terminators.</li>
% </ul></ul>
% </html>

%%
% * *Segments*: multiple sets of lines held in a single cell array.
%%
% <html>
% <ul><ul>
% <li><i>cplines</i> - a cell array of <i>plines</i> (each of any length and any number of lines).</li>
% </ul></ul>
% </html>

%%
% * *Sections*: a special case of a _pline_ with only 2 _points_ (start and
% end) and a NaN _point_.
%%
% <html>
% <ul><ul>
% <li><i>pslice</i> - a <i>pline</i> defined by two end <i>points</i>.</li>
% <li><i>pslices</i> - multiple <i>pslice</i> segments concatenated as a row struct array, 
% with individual <i>pslice</i> segments defined by NaN terminators.</li>
% </ul></ul>
% </html>

%%
% Hence a _line_ and _lines_ are always a set of column vectors in a struct,
% table or matrix. Whereas the other formats are all different arrays of a _point_ set.

%%
% <html>
% <table border=1><tr><td><u>Note</u>:<br> 
% <li>To extract a <i>point</i> from <i>points</i> use:   point = points(i);</li>
% <li>To extract the ith <i>pline</i> of points from a <i>plines</i> use:</li>
% &emsp; &emsp; idN = [0,find(isnan([points(:).x]))];<br>
% &emsp; &emsp;   for i=1:length(idN)-1<br>
% &emsp; &emsp; &emsp; pline{1,i} = points(idN(i)+1:idN(i+1));<br>
% &emsp; &emsp;   end<br>
% <li>????o extract a <i>point</i> from <i>plines</i> use: point = plines(i)???? ;</li>
% </td></tr></table>
% </html>

%% Classes
% The Apps that inherit <matlab:doc('gdinterface') GDinterface> can make use of 
% *GD_Sections* to manipulate points and lines, and extract boundaries and
% cente-lines from an imported digital terrain model. This extends the
% capability of the functions provided as part of the standard
% <matlab:doc('grid_menu_help') Grid tools Menu>. 

%%
% *GD_Sections*: provides methods to manipulate points, lines and sections 
% using datasets in classes that inherit <matlab:doc('gdinterface') GDinterface>.

%%
% Properties for sections:
%%
% * _Boundary_ - enclosing boundary used to crop length of sections.
% * _ChannelLine_ - centre-line of channel can be multiple lines to define
% a network.
% * _ChannelProps_ - properties used when extracting the channel network.
% * _SectionLines_ - sections defined at right angles to points on
% centre-line extending to Boundaray line (if cropping applied).
%%
% All properties are held as _lines_ and the GD_Sections instance is saved 
% to the Sections property in the calling class (e.g. EDBimport).

%%
% Methods to manipulate sections:
%% 
% * setBoundary - load a grid and interactively define a boundary line with
% options to autogenerate and edit or digistise the line.
% * setChannelNetwork - load a grid and extract the centre lines of the channels in
% a network, or digitise them manually, and combine into a vector of points.
% * setSectionLines - use the centre-line to generate a set of sections at right
% angles and interactively edit them, or digitise them manually.
% * setSections - using the SectionLines data extract the widths as
% a function of elevation and distance along the centre-line.
% * viewsections - view boundary channel network and cross-sections line work.

%%
% Static methods to manipulate sections:
%%
% * editLines - edit linework for selected line type (can be Boundary,
% Channel Network or Section Lines).
% * loadLines -  load linework for selected line type from a shapefile.
% * getSection - retrieve existing class object or create a new instance.
% * sectionMenuOptions - default set of menu options for use in Model UIs
% such as EstuaryDB. Can be used to call any of the methods listed above


%% Functions
% Functions that derive and manipulte sections from a Grid 
% such as contour boundaries, channel centre-lines can be found in the _muiAppGridFcns_ 
% folder. Use the Matlab(TM) _help_ function in the command window to get 
% further details of each function.

%%
% *gd_boundary.m*
% - UI to interactively generate a contour line (e.g. to define a boundary).
%%
%   blines = gd_boundary(grid,paneltxt,outype,isdel);
%%
% where _grid_ is a struct of x, y, z, _paneltext_ is a character string
% used for title, _outype_ is the format of output (see *gd_pnt2vec*), _inlines_
% is struct or table of x,y vectors to be edited, or the format 
% of output if no lines are being input, and _isdel_ is a logical flag true to delete figure 
% on completion (optional - default is false). The output, _blines_, 
% contains _lines_ of the extracted contour.

%%
% *gd_centreline* - create a centreline of a channel to trace the
% shortest path between start and end points whilst finding the deepest
% points (i.e. a thalweg).
%%
%   cline = gd_centreline(grid,mobj,props,clines);
%
%%
% where _grid_ is a struct of x, y, z, _mobj_ is a muitoolbox
% model instance, _props_ is a struct for maximum water level and depth exponent 
% (maxwl and dexp), _clines_ is optional and passes any previously defined
% centre-lines. The output, _cline_, is a _line_ of the extracted centre-line.

%%
% *gd_digitisepoints.m*
% - UI to interactively digitise x,y,z points on a grid and edit
% z elevations, if required
%%
%   lines = gd_digitisepoints(grid,paneltxt,outype,isxyz,isdel)
%%
% where _grid_ is a struct of x, y, z, _paneltext_ is a character string
% used for title, _outype_ is the format of output (see *gd_pnt2vec*), _isxyz_
% is a logical flag true to input z values (optional, default is false)
% and _isdel_ - logical flag true to delete figure 
% on completion (optional - default is false). Output is a _lines_, or 
% _plines_ array, depending on _outype_.

%%
% *gd_editlines.m*
% - UI to interactively digitise x,y,z points on a grid and edit
% z elevations, if required
%%
%   lines = gd_editlines(grid,paneltxt,inlines,isdel);
%%
% where _grid_ is a struct of x, y, z, _paneltext_ is a character string
% used for title, _outype_ is the format of output (see *gd_pnt2vec*), _inlines_
% is struct or table of x,y vectors to be edited, or the format 
% of output if no lines are being input, and _isdel_ is a logical flag true to delete figure 
% on completion (optional - default is false). Output is a _lines_, or 
% _plines_ array, depending on _outype_.

%%
% *gd_findline.m*
% - find which _pline_ a given _point_ lies on when there are multiple lines
% separated by NaN values in a _plines_ struct array.
%%
%   lineIndex = gd_findline(plines,queryPoint);
%%
% where _plines_ - is a struct array of x,y points defining lines and each line is 
% separated by a NaN x,y _point_, and _queryPoint_ is an x,y struct of the
% _point_ to be tested. The output _lineIndex_ is an index to the line in
% _plines_ on which the query _point_ lies.

%%
% *gd_getcontour.m*
% - extract a contour at a defined level.
%%
%  clines = gd_getcontour(grid,zlevel)
%%
% where _grid_ - struct of x, y, z, _zlevel_ is the level of contour to be extracted
% and _isplt_ is a logical flag true to create plot (optional, default is
% false). The output, _clines_, is a set of _lines_ of the extracted contour.

%%
% *gd_getpoint.m*
% - interactively select a single point in a UI figure and return the point
% coordinates.
%%
%   point = gd_getpoint(ax,promptxt)
%%
% where _ax_ is a figure axes used to interactively select point, and _promptxt_
% is a prompt to be used for point being selected. The output is a _point_ struct.

%%
% *gd_plotgrid*
% - create pcolor plot of gridded surface (from Grid Tools).
%%
%   ax = gd_plotgrid(hfig,grid);  %hfig is figure handle and ax is axes handle

%%
% *gd_plotsections*
% - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure (from Grid Tools).
%%
%   gd_plotsections(grid);  %where grid is a struct of x, y, z


%%
% *gd_points2lines.m*
% - convert a _lines_ array of structs with x, y (and z) fields to a [Nx2] or [Nx3] 
% array, or a single stuct with column vectors in the x, y (and z) fields
%%
%   lines = gd_point2lines(points,outype);
%%
% where _outype_ determines the format of the output, _lines_ as follows:
%%
% <html>
% <ul><ul>
% <li>outype=0: array of structs with x, y (and z) fields defining _points_.</li>
% <li>outype=1: [Nx2] or [Nx3] array.</li>
% <li>outype=2: struct with x, y (and z) column vectors</li>
% <li>outype=3: table with x, y (and z) column vectors.</li>
% </ul></ul>
% </html>

%%
% *gd_sectionlines.m*
% -extract the section lines that are normal to the channel centre line
% and extend to the bounding shoreline.
%%
%   [slines,clines] = gd_sectionlines(obj,cobj,paneltxt,isdel);
%%
% where _obj_ is an instance of GD_Sections with a Boundary and ChannelLine 
% defined and _cobj_ is an instance of EDBimport class with grid or geoimage
% _paneltext_ is a character string used for title, _outype_ is the format
% of output (see *gd_pnt2vec*), and _isdel_ is a logical flag true to delete figure 
% on completion (optional - default is false). The output, _sleines_ and 
% _clines_, sets of _lines_ of the extracted sections and centreline (can
% be smoothed as part of workflow).

%%
% *gd_selectpoints.m*
% - UI to interactively create a specified number of x,y point on a grid.
%%
%   [points,hfig] = gd_selectpoints(grid,paneltxt,promptxt,inlines;npts,outype,isdel);
%%
% where _grid_ is a struct of x, y, z, _paneltext_ is a character string
% used for title, _promptxt_ is a cell array of prompts to be used for each 
% point being defined (if a single cell the text is apended with a count
%  for each point, othewise the cell array should have a length of npts),
% _inlines_ - struct, table or matrix of x,y column vectors to show on base
% plot, _npts_ is the number of points to be selected, _outype_ - format of 
% output (see *gd_pnt2vec*) and _isdel_ - logical flag true to delete figure 
% on completion (optional - default is false). Output is a _points_ array and
% _hfig_ a handle to the UI figure;

% *gd_setcplines.m*
% - %   converts a set of points to a set of _plines_ based on NaN separators,
% plots the graphical lines, or digistise a set of points and return as a
% _pline_.
%%
%   cplines = gd_setcplines(ax,promptxt,linepoints);
%%
% where _ax_ is a figure axes used to interactivly select points, _promptxt_
% is a prompt to be used for point being defined, and _linepoints_ is a
% set of _points_ to be set as lines x,y struct with each line
% define by a NaN separator (optional - user prompted to define points if not set)
% cplines - a cell array of _plines_, or anew digitised _pline_.

%%
% *gd_setpoint.m*
% - interactively create a single point in a UI figure and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. elevation).
%%
%   point = gd_setpoint(ax,promptxt,isxyz)
%%
% where _ax_ is a figure axes used to interactively select points, _promptxt_
% is a prompt to be used for point being defined, and _isxyz_
% is a logical flag true to input z values (optional, default is false).
% The output is a _point_ struct.

%%
% *gd_setpoints.m*
% - interactively create a set of points in a UI figure and return the point
% coordinates. Includes an option to enter an additional value at the
% selected points (e.g. elevation).
%%
%   points = gd_setpoints(ax,promptxt,isxyz);
%%
% where _ax_ is a figure axes used to interactively select points, _promptxt_
% is a prompt to be used for point being defined, and _isxyz_
% is a logical flag true to input z values (optional, default is false).
% The output is a _points_ array.


%%
% *gd_smoothlines.m*
% - smooth one or more line segments using either moving average, or the
% Savitzky-Golay method.
%%
%   smoothlines = gd_smoothlines(lines,method,win,deg,npnts)
%%
% where _lines_ is a struct of x,y vectors to be smoothed, with each line
% segment separated by a NaN, _method_ can be 'movmean' for moving average,
% or 'sgolay' for Savitzky-Golay smoothing, _window_ is the size of window 
% to use (1 or 2 elements), _degree_ is the  Savitzky-Golay degree
% (<window) and _npnts_ is the minimum number of points required to apply smoothing

%%
% *gd_lines2points.m*
% convert _lines_  as x,y vectors in various formats (see *gd_pnt2vec*) to a 
% _points_ array of structs for each _point_ with x, y (and z) fields.
%%
%   [points,outype] = gd_lines2points(lines);


















