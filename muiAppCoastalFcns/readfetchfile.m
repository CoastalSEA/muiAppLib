function [Fdir,Flen]  = readfetchfile()
%
%-------function help------------------------------------------------------
% NAME
%   readfetchfile.m
% PURPOSE
%   Read the file that contains the fetch lengths as a function of direction
% USAGE
%   [Fdir,Flen]  = readfetchfile()
% INPUTS
%   none
% OUTPUT
%   Fdir - fetch directions (degTN)
%   FLen - fetch lenghts (m)   
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%

    [fname, path]=uigetfile('*.txt','Wind-Wave Fetch file');
    fid = fopen([path,fname], 'r');
    if fid<0
        errordlg('Could not open file for reading','File write error','modal')
        Fdir = []; Flen = [];
        return;
    end

    ncols = 2;
    headSpec = '%s';
    %read header
    header = textscan(fid,headSpec,ncols); %concatenates header to single string
    %read numeric data - Dir Fetch
    dataSpec = '%f %f';
    data = textscan(fid,dataSpec);
    Fdir = data{:,1};
    Flen = data{:,2};
    fclose(fid);
end