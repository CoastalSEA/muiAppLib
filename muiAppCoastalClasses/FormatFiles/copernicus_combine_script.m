% copernicus_combine_script.m
% Computes combined Hs, representative Tp and mean direction from Copernicus NWSHELF bulk fields
% Requires: MATLAB with netcdf support

[fname,path,~] = getfiles('PromptText','Copernicus netcdf file','FileType','*.nc');
if fname==0, formatfile = []; return; end  %user cancelled
ncfile = [path,fname];
% Read variables (use try/catch to allow missing fields)
vars = {'VHM0','VHM0_SW1','VHM0_SW2','VHM0_WW',...
        'VMDR','VMDR_SW1','VMDR_SW2','VMDR_WW','VPED',...
        'VSDX','VSDY',...
        'VTM01_SW1','VTM01_SW2','VTM01_WW',...
        'VTM02','VTM10','VTPK'};
for v=1:numel(vars)
    try
        data.(vars{v}) = ncread(ncfile, vars{v});
    catch
        data.(vars{v}) = [];
    end
end

dur = ncread(ncfile,'time');
date = datetime('1970-01-01 00:00:00');
myDatetime =  date+seconds(dur);
myDatetime.Format = 'dd-MM-yyyy HH:mm:ss';

% Ensure arrays align: assume time x lat x lon or similar; operate elementwise
% Example: compute combined Hs using RSS if component Hs exist
Hs_sw1 = data.VHM0_SW1; Hs_sw2 = data.VHM0_SW2; Hs_ww = data.VHM0_WW;
if ~isempty(Hs_sw1) || ~isempty(Hs_sw2) || ~isempty(Hs_ww)
    % replace missing with zeros
    if isempty(Hs_sw1), Hs_sw1 = 0*data.VHM0; end
    if isempty(Hs_sw2), Hs_sw2 = 0*data.VHM0; end
    if isempty(Hs_ww),  Hs_ww  = 0*data.VHM0; end
    Hs_comb_RSS = sqrt(Hs_sw1.^2 + Hs_sw2.^2 + Hs_ww.^2);
else
    Hs_comb_RSS = data.VHM0; % fallback to provided total if present
end

% Energy-weighted peak period (approx)
Tp_sw1 = data.VTM01_SW1; Tp_sw2 = data.VTM01_SW2; Tp_ww = data.VTM01_WW;
% set defaults and masks
E_sw1 = (Hs_sw1).^2; E_sw2 = (Hs_sw2).^2; E_ww = (Hs_ww).^2;
E_sum = E_sw1 + E_sw2 + E_ww;
% avoid divide by zero
E_sum(E_sum==0) = NaN;
Tp_energy = (E_sw1.*Tp_sw1 + E_sw2.*Tp_sw2 + E_ww.*Tp_ww) ./ E_sum;

% Representative mean direction: circular energy-weighted mean
dir_sw1 = deg2rad(data.VMDR_SW1); dir_sw2 = deg2rad(data.VMDR_SW2); dir_ww = deg2rad(data.VMDR_WW);
Cx = E_sw1.*cos(dir_sw1) + E_sw2.*cos(dir_sw2) + E_ww.*cos(dir_ww);
Cy = E_sw1.*sin(dir_sw1) + E_sw2.*sin(dir_sw2) + E_ww.*sin(dir_ww);
theta_mean = atan2(Cy, Cx); % radians
theta_mean_deg = mod(rad2deg(theta_mean),360);

% Fallback: if component directions missing, use VPED or Stokes drift
if all(isnan(theta_mean_deg(:))) || isempty(data.VMDR)
    if ~isempty(data.VPED)
        theta_mean_deg = data.VPED;
    elseif ~isempty(data.VSDX) && ~isempty(data.VSDY)
        theta_mean_deg = mod(rad2deg(atan2(data.VSDY, data.VSDX)),360);
    end
end

%plot results
t = myDatetime;
hf = figure('Tag','PlotFig');
ht = tiledlayout(hf,3,1,"TileSpacing","compact","Padding","compact");
nexttile
plot(t,Hs_comb_RSS(:))
ylabel('Wave height')
nexttile
plot(t,Tp_energy(:))
ylabel('Wave period')
nexttile
plot(t,theta_mean_deg(:))
ylabel('Wave direction')
xlabel('Time')
sgtitle(fname)

%plot differences with sea surface values in Copernicus file
H = Hs_comb_RSS(:)-data.VHM0(:);
T = Tp_energy(:)-data.VTPK(:);
D = theta_mean_deg(:)-data.VMDR(:);
t = myDatetime;
hf = figure('Tag','PlotFig');
ht = tiledlayout(hf,3,1,"TileSpacing","compact","Padding","compact");
nexttile
plot(t,H(:))
ylabel('Height difference')
nexttile
plot(t,T(:))
ylabel('Period difference')
nexttile
plot(t,D(:))
ylabel('Direction difference')
xlabel('Time')
sgtitle('Differences: Recombined-Source values')