function gd_user_function(obj,mobj)
%
%-------function help------------------------------------------------------
% NAME
%   gd_user_function.m
% PURPOSE
%   function for user to define bespoke use of grids and grid tools
% USAGE
%   gd_user_function(obj,mobj)
% INPUTS
%   obj - instance of a class that inherits GDinterface abstract class
%   mobj - instance of model class (ie a class that inherits muiModelUI)
% OUTPUT
%   user defined
% NOTES
%   provided as an option on the default grid tools menu to allow the user
%   to rapidly add bespoke tools
% SEE ALSO
%   called from GDinterface.getMenuOptions
%   
% Author: Ian Townend
% CoastalSEA (c) Oct 2022
%--------------------------------------------------------------------------
%

    %add code or call function: 
    % userfunction(obj,mobj);      %<< rename to required function
    warndlg('Add code, or function call, to ''gd_user_function''')