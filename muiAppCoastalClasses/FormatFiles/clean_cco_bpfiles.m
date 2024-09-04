function clean_cco_bpfiles
%
%-------function help------------------------------------------------------
% NAME
%   clean_cco_bpfiles.m
% PURPOSE
%   function to clean up CCO beach profile data that has 8 columns
%   interspersesd with 7 column data
% USAGE
%   clean_cco_bpfiles
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
userprompt = 'Select data file(s)>';
[fname, path]=uigetfile('*.txt',userprompt,'MultiSelect','on');
if isequal(fname,0)
    return; 
elseif ischar(fname)
    fname = {fname};
end

for kk = 1:numel(fname)
    % Read the file
    fid = fopen([path,fname{kk}],'r');
    data = textscan(fid,'%s','Delimiter','\n');
    fclose(fid); 

    %do things with the data
    txtdata = data{1};
    delimiter = sprintf('\t'); %tab 
    %find number of columns in header
    ncols = numel(split(txtdata{1},delimiter));
    for i=2:length(txtdata)
        %check whether each row has correct number of columns
        rowtxt = split(txtdata{i});
        ndata= numel(rowtxt);
        if ndata~=ncols
            %remove extra column - penultimate column of input row
            rowtxt(ncols) = [];
            txtdata(i) = join(rowtxt,delimiter);
        end
    end
    
    % Save as a text file
    fnam = sprintf('clean-%s',fname{kk});
    fid = fopen([path,fnam],'w');
    fprintf(fid,'%s\n', txtdata{:});
    fclose(fid);
end
warndlg('All files cleaned and written to output')