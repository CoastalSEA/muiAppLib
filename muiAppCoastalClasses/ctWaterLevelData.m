classdef ctWaterLevelData < muiDataSet                   
%
%-------class help------------------------------------------------------===
% NAME
%   ctWaterLevelData.m
% PURPOSE
%   Class to illustrate importing a data set, adding the results to dstable
%   and a record in a dscatlogue (as a property of muiCatalogue)
% USAGE
%   obj = ctWaterLevelData()
% SEE ALSO
%   inherits muiDataSet and uses dstable and dscatalogue
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%    
    properties  (Transient)
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        % importing data requires muiDataSet propertiesm DataFormats and
        % FileSpec to be defined in class constructor.
        %Additional properties:  
        datum = 0 %datum used for CD to OD conversion
        datumnotset = true; %flag to re-use entered datum for all files
    end

    methods 
        function obj = ctWaterLevelData()            
            %class constructor
            %initialise list of available input file formats. Format is:
            %{'label 1','function name 1';'label 2','function name 2'; etc}
            obj.DataFormats = {'Date-Record format','wl_daterec_format';...
                               'Channel Coastal Observatory format','wl_cco_format';...
                               'BODC NTSLF format','wl_bodc_ntslf_format';...
                               'BODC ODV format','wl_bodc_odv_format';...
                               'GESLA format','wl_gesla_format';...                               
                               'ShoreCast data format','wl_scast_format'};
            %define file specification, format is: {multiselect,file etnsion types}
            obj.FileSpec = {'on','*.txt;*.csv'};
        end
%%
        function tabPlot(obj,src)
            %generate plot for display on Q-Plot tab
            tabDefaultPlot(obj,src);
        end 
%%
        function wldata = addDatum(obj,data,isreqd)
            %get OD to CD offset and add data column for OD
            %data - table of variable data (excluding datatime column)
            %isreqd - logical flag to indicated if adjustment required
            %wldata - insert adjusted variable as first column in table
            if ~isreqd
                answer = questdlg('Adjust datum?','Import WL data',...
                                             'Yes','No','No');
                if strcmp(answer,'No'), wldata = data; return; end
            end
                
            if obj.datumnotset
                prompt = {'OD relative to CD (should be positive):',...
                              'Use this value for all files to be loaded'};
                title = 'Datum correction';
                numlines = 1;
                default = {num2str(obj.datum),'Yes'};
                ok = false;
                while ~ok
                    answer = inputdlg(prompt,title,numlines,default);
                    if isempty(answer) %user cancels
                        %force user to give a positive answer because data load
                        %formats assume a column for OD and CD
                    elseif str2double(answer{1})>=0
                        obj.datum = str2double(answer{1});
                        if strcmp(answer{2},'Yes')
                            obj.datumnotset = false;
                        end
                        ok = true;
                    end
                end
            end
            %water levels to chart datum are in first column of table
            %to be consistent with CCO data need to make first column ODN 
            %and the rest moved to right
            zodn = data{:,1}-obj.datum;
            zodn(abs(data{:,1})==99) = 99;
            wldata = addvars(data,zodn,'Before',1);
%             wldata = data(1,1:2);
%             wldata = horzcat(wldata,zodn,data(1,3:end));
        end
    end
end