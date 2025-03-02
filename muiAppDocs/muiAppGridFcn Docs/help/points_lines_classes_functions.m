%% Points and Lines
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
% _gd_points2lines_ and gd_lines2points to convert between points and other formats.
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
% Hence a _line_ and _lines_ are always a set of column vectors in a struct,
% table or matrix. Whereas the other formats are all different arrays of a _point_ set.

%%
% <html>
% <table border=1><tr><td><u><b>Notes</b></u>:<br><br>
% <li>To extract a <i>point</i> from <i>points</i> use: &emsp; point = points(i);</li><br>
% <li>To extract the ith <i>pline</i> of points from a <i>plines</i>
% use:</li><br>
% &emsp; &emsp; idN = [0,find(isnan([points(:).x]))];<br>
% &emsp; &emsp;   for i=1:length(idN)-1<br>
% &emsp; &emsp; &emsp; pline{1,i} = points(idN(i)+1:idN(i+1));<br>
% &emsp; &emsp;   end<br><br>
% <li>To switch between plines and cplines use <i>gd_plines2cplines</i>and <i>gd_cplines2plines</i>;</li><br>
% <li>To switch between point, or points, and line, or lines, use
% <i>gd_points2lines</i> and <i>gd_lines2points</i>.</li><br>
% <li>When using lines, some functions use row vectors and others use
% column vectors. <br>&emsp; To switch an xy struct from column to row or vice versa use:<br>
%  &emsp; &emsp; lines = structfun(@transpose,lines,'UniformOutput',false); <br><br>
% </td></tr></table>
% </html>

%% Point and Line classes
% The Apps that inherit <matlab:doc('gdinterface') GDinterface> can make use of 
% *PL_Sections* to manipulate points and lines, and extract boundaries and
% cente-lines from an imported digital terrain model. This extends the
% capability of the functions provided as part of the standard
% <matlab:doc('grid_menu_help') Grid tools Menu>. 

%%
% *PL_Sections*: provides methods to manipulate points, lines and sections 
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

%%
% To manipulate points and lines, a number of classes make use of the 
% PLinterface abstract class, see <matlab:doc('plinterface') PLinterface>
% for further details. The most generic is *PL_Editor* and the input and
% outputs defined for this class are also used in the other _PL_ classes
% that provide bespoke functionality.
%%
% *PL_Editor*:
% figure to interactively digitise points on a grid and add
% elevations if required.
%% 
%   [lines,points] = PL_Editor.Figure(grid,promptxt,inlines,,isxyz,isdel);
%% 
% where _grid_ is a struct of x, y, z (e.g. as used in getGrid in the
% GDinterface); _promptxt_ character string used for initial prompt in title; _inlines_ 
% a struct of points and lines or a struct of x,y vectors to be edited, or
% an output format flag, _isxyz_ is a logical flag true to input z values -
% optional, default is false; and _isdel_ is a logical flag true to delete
% figure on completion - optional, default is false. The ouput format
% depends on the _outype_ flag. If _outype_=0: array of structs with x, y and z 
% fields defining selected points, _outype_=1: Nx2 or Nx3 array,
% _outype_=2: struct with x, y (and z) vector fields, and _points_ = [] if 
% user closes figure, or no points defined.

%%
% *PL_Boundary*: Class to extract contours and generate model boundaries
%%
%   lines = PL_Boundary.Figure(grid,promptxt,inlines,isdel);
%%
% where the inputs are as for *PL_Editor*, with the addition of _linlines_ for any
% existing boundary lines.

%%
% *PL_CentreLine*: Class to extract valley/channel centre-line
%%
%   [lines,props] = PL_CentreLine.Figure(grid,promptxt,inlines,props,isdel);
%%
% where the inputs are as for *PL_Editor*, with the addition of _linlines_ for any
% existing centre-line and _props_ struct for maximum water level, depth exponent
% and the sampling interval along the centre-lines (maxwl,dexp,cint). The
% variable _props_ is included in the output to capture any interactive updates.

%%
% *PL_PlotSections*:
% display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
%%
%   PL_PlotSections.Figure(grid,promptxt,isdel); 
%%
% where the inputs are as for *PL_Editor*.

%%
% *PL_SectionLines*: Class to create a set of cross-sections lines based on 
% the points in a centre-line and any enclosing boundary.
%%
%   [slines,clines] = PL_SectionLines.Figure(grid,promptxt,setlines,isdel);
%%
% where the inputs are as for *PL_Editor*, with the addition of _setlines_,
% a struct with fields Boundary, ChannelLine, ChannelProps and SectionLines, 
% or an instance of the *PL_Sections* class. The output includes _slines_
% for the sections defined and _clines_ for the updated centre-line.

%% Point and Line Functions
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
% used for title, _outype_ is the format of output (see *gd_ppoints2lines*), _inlines_
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
% *gd_cplines2plines.m*
% - convert cell array of plines to an array of points that define plines.
%%
%   plines = gd_cplines2plines(cplines);
%%
% where _cplines_ is a cell array of plines and _plines_ is a struct array 
% with x and y fields defining a set of points.

%%
% *gd_curvelineprops*
% - for each point from idL to the end use the centre-line coordinates and 
% direction to find the lengths and directions along the centre-line.
%%
%   [clinedir,ncplines,cumlen] = gd_curvelineprops(cplines,idL)
%%
% where _cplines_ is a cell array of plines and _idL_ is the index of the
% start point on any of the lines in _cplines_. The output includes 
% _clinedir_, the mean direction of line at each point in x,y space, 
% _ncplines_ an updated version of _cplines_ starting from point _idL_ and
% _cumlen_, the cumulative length from the defined start point.

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
% *gd_getpline.m*
% - interactively select a line on a plot and return the line point coordinates.
%%
%   [pline,H] = gd_getpline(ax,promptxt,tagname,ispoints);
%%
% where _ax_ is the figure axes used to interactively select point, 
% _promptxt_ is the prompt to be used for point being defined, _tagname_ is
% a character vector of text to be used as tag for plotted points, _ispoints_ 
% is true to return as an array of point structs, otherwise returns an xy 
% struct of points (optional - default is true). Returns _pline_ as a 
% struct array of points defining selected line or a struct with x and y 
% fields, depending on value of _ispoints_ and _H_ is a handle to selected line

%%
% *gd_getpoint.m*
% - interactively select a single point in a UI figure and return the point
% coordinates.
%%
%   point = gd_getpoint(ax,promptxt);
%%
% where _ax_ is a figure axes used to interactively select point, and _promptxt_
% is a prompt to be used for point being selected. The output is a _point_ struct.

%%
% *gd_lines2points.m*
% - convert _lines_  as x,y vectors in various formats (see *gd_points2lines*) to a 
% _points_ array of structs for each _point_ with x, y (and z) fields.
%%
%   [points,outype] = gd_lines2points(lines);

%% 
% *gd_linetopology*
% - interactively define line connectivity.
%%
%   [cumlen,G,hf,hg] = gd_linetopology(grid,plines);
%%
% where  _grid_ is a struct of x, y, z, _plines_ is a struct of x,y vectors 
% defining one or more lines. The outputs include _cumlen_, the cumulative 
% lengths along lines from first point in line set, _G_ a directed graph of 
% the channl network, _hf_ a handle to figure used to define links and _hg_ 
% a handle to figure of the resultant network.

%%
% *gd_orderplines.m*
% - amend the order of the _plines_ in an array of _plines_
%%
%   plines = gd_orderlines(ax,plines);
%%
% where _ax_ is the figure axes being used to edit points and _plines_ is 
% a struct of x,y vectors defining one or more lines. The output is a copy
% of _plines_ in the revised order.

%%
% *gd_plines2cplines.m*
% - convert an array of _plines_ to a cell array of _plines_.
%%
%   cplines = gd_plines2cplines(plines);
%%
% where _plines_ is a struct array with x and y fields defining one or more
% plines, and _cplines_ is a cell array of plines.

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
%   [xlines,hf] = gd_plotsections(grid,cplines,inp)
%%
% where _grid_ is a struct of x,y,z values that define grid, _cplines_ is a
% cell array of section lines and _inp_ is a struct with fields for zmax, 
% sint and method, optional if empty user is prompted to enter the values 
% in an input dialog. the output includes _xlines_ for the interpolated 
% distances and elevations along section lines and _hf_, the handle to plot of sections

%%
% *gd_points2lines.m*
% - convert a _points_ or plines array of structs with x, y (and z) fields 
% to a [Nx2] or [Nx3] array, or a single stuct, or matrix with column vectors
%  in the x, y (and z) fields.
%%
%   lines = gd_points2lines(points,outype);
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
% - extract the section lines that are normal to the channel centre line
% and extend to the bounding shoreline.
%%
%   [s_lines,c_lines] = gd_sectionlines(obj,cobj,paneltxt,isdel);
%%
% where _obj_ is an instance of GD_Sections with a Boundary and ChannelLine 
% defined and _cobj_ is an instance of EDBimport class with grid or geoimage
% _paneltext_ is a character string used for title, _outype_ is the format
% of output (see *gd_points2lines*), and _isdel_ is a logical flag true to delete figure 
% on completion (optional - default is false). The output, _s_lines_ and 
% _c_lines_, sets of _lines_ of the extracted sections and centreline (can
% be smoothed as part of workflow).

%%
% *gd_selectpoints.m*
% - UI to interactively create a specified number of x,y point on a grid.
%%
%   points = gd_selectpoints(grid,paneltxt,promptxt,inlines;npts,outype,isdel);
%%
% where _grid_ is a struct of x, y, z, _paneltext_ is a character string
% used for title, _promptxt_ is a cell array of prompts to be used for each 
% point being defined (if a single cell the text is apended with a count
%  for each point, othewise the cell array should have a length of npts),
% _inlines_ - struct, table or matrix of x,y column vectors to show on base
% plot, _npts_ is the number of points to be selected, _outype_ - format of 
% output (see *gd_points2lines*) and _isdel_ - logical flag true to delete figure 
% on completion (optional - default is false). Output is a _points_ array.

%%
% *gd_setpoint.m*
% - interactively create a single point in a UI figure and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. elevation).
%%
%   point = gd_setpoint(ax,promptxt,isxyz);
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


















