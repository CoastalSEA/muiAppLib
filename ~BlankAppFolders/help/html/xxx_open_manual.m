function xxx_open_manual()
%find the location of the asmita app and open the manual
appname = 'ModelName';
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},appname));
fpath = [appinfo(idx(1)).location,[filesep,'doc',filesep,'ModelName manual.pdf']];
open(fpath)
