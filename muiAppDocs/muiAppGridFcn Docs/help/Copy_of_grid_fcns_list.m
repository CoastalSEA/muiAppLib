%% Grid Functions
% Functions used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and include the following:
%%
% * *gd_ax_dir*
% - check direction of grid axes and reverse if descending, OR
% find grid orientation using ishead and direction of x-axis, OR
% check a grid axis direction by prompting user.
% * *gd_basin_hypsometry*
% - compute area and volume hypsometry from gridded elevation data.
% * *gd_basin_indices*
% - get the indices of the grid x-axis that fall within the basin or channel,
% when the mouth is offset from the grid origin. (NB: assumes basin/channel
% is aligned iwth the x-axis). Also returns the index of mouth position on 
% the x-axis.
% * *gd_basin_properties*
% - use the basin hypsometry from gd_basin_hypsometry to compute several 
% along-channel/x-axis morphological properties.
% * *gd_basin_volumes.m*
% - compute area and volume totals over x-z from hypsometry of gridded
% elevation data.
% * *gd_centreline.m*
% - create a centreline of a channel using function _a_star_ to trace the
% shortest path between start and end points whilst finding the deepest
% points (i.e. a thalweg).
% * *gd_colormap*
% - check if Mapping toolbox is installed to use land/sea colormap, or call
% _cmap_selection_ if not available (see <matlab:doc('psfunctions') Plotting and statistical functions> 
% in the <matlab:doc('muitoolbox') muitoolbox>).
% * *gd_convergencelength* 
% - least squares fit using fminsearch to % find the convergence length of 
% a channel from a distance-width xy data set.
% * *gd_digitisepoints*
% - creates figure to interactively digitise points on a grid and add
% elevations if required.
% * *gd_dimensions*
% - get the grid dimsnions for a grid struct (as used in GDinterface).
% * *gd_getpoint.m*
% - interactively select a point on a plot and return the point
% coordinates.
% * *gd_grid_line*
% - create a grid from scattered data points input as xyz tuples.
% * *gd_gross_properties*
% - compute the gross properties of a gridded bathymetry.
% * *gd_lineongrid_plot*
% - plots a defined line onto a countour or surface plot of a grid (e.g a
%   channel centre-line).
% * *gd_plan_form* 
% - compute planform variation along the x-axis at specified planar levels.
% * *gd_plotgrid*
% - create pcolor plot of gridded surface.
% * *gd_plotsections*
% - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
% * *gd_pnt2vec.m*
% - convert an array of structs with x,y (and z) fields to a [Nx2] or [Nx3] 
% array of points, or a single stuct with vectors for the x, y (and z)
% fields.
% * *gd_property_plots*
% - plots displayed on Proprety tab or stand-alone figure in Apps that use 
% GDinterface, such as ChannelForm and ModelSkill.
% * *gd_readshapefile.m*
% - read the x and y coordinates from a shape file. Lines are concatenated
% and separated by NaNs in single x and y vectors. Suitable for reading
% boundaries or sections into a single array.
% * *gd_section_properties*
% - compute the width, cross-sectional area and prism along channel.
% * *gd_selectpoints*
% - accept figure to interactively create a specified number of x,y points
% on a grid.
% * *gd_setpoint*
% - interactively select a single point on a plot and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. for elevation).
% * *gd_setpoints.m*
% - interactively create a set of points on a plot and return the point
% coordinates. Includes an option to enter an additional value at the
% selected points (e.g. elevation).
% * *gd_startendpoints*
% - accept figure to interactively select start and end points on a grid.
% * *gd_subdomain*
% - accept figure to interactively select a subdomain of a grid.
% * *gd_property_plots* 
% - plots displayed on Proprety tab in ChannelForm model and on a figure 
% in ModelSkill.
% * *gd_xy2sn*
% - map grid from cartesian to curvilinear coordinates with option to return 
% the elevations on the source cartesian grid, or as a curvilinear grid.
% * *gd_sn2xy*
% - map grid from curvilinear to cartesian coordinates.
% * *gd_user_function*
% - function for user to define bespoke use of grids and grid tools.
% * *getsubgrid*
% - extract a subdomain from a grid and return the extracted grid and the 
% source grid indices of the bounding rectangle.
%%
% Functions from Matlab(TM) Exchange Forum
%%
% * *a_star*
% - implements the A* search algorithm to find the shortest path given
% constraints (inaccessible cells) and a cost function (e.g. water depths).
% Author: Alex Ranaldi, 2022, https://github.com/alexranaldi/A_STAR
% * *InterX* 
% - intersection of two curves. MATLAB Central File Exchange, 
% Author: NS, 2010, https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections.
% * *xy2sn* 
% - Bart Vermeulen,2022, Cartesian to Curvilinear 
%   coordinate forward and backward transformation. 
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation 
% * *sn2xy* 
% - as above.

%%
% Further details can be found in <matlab:doc('grid_class_fcns') Grid classes and functions>
% 