%   Script to test xy2sn and sn2xy - Juernjakob Dugge (2015). jdugge/xy2sn, 
%   https://github.com/jdugge/xy2sn, which requires 
%   arclength - John D'Errico, 2012, https://www.mathworks.com/matlabcentral/fileexchange/34871-arclength 
%   distance2curve - John D'Errico, 2013, http://www.mathworks.de/matlabcentral/fileexchange/34869-distance2curve 
%   interparc - John D'Errico, 2012, https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc 

% Generate test data
dataX = sin( (0:0.01:1)' *2*pi) + (rand(101,1) - 0.5) * 0.4;
dataY = (0:0.02:2)' - 0.4*sin((0:0.01:1)' *4*pi)  + ...
    (rand(101,1)-0.5) * 0.4;

% figure; scatter(dataX, dataY);

% Provide centerline path
centerlineX = sin( (0:0.05:1)' *2*pi);
centerlineX = [-0.2; centerlineX; 0.2];

centerlineY = (0:0.1:2)' - 0.4 * sin( (0:0.05:1)'*4*pi);
centerlineY = [0.1; centerlineY; 1.9];

% hold on
% plot(centerlineX, centerlineY);
% hold off

% Transform to flow-oriented S-N coordinate system
[S, N, ~] = xy2sn(dataX, dataY, centerlineX, centerlineY);

% Plot points in S-N coordinate system
% figure; scatter(S, N, [], S, 'f')
% axis equal


% Plot points in X-Y coordinate system, using S for colour
% figure; scatter(dataX, dataY, [], S, 'f')
% axis equal

%% Generate a grid in the S-N coordinate system and transform it to X-Y
L = S(end); mnn = floor(min(N)*10)/10; mxn = ceil(max(N)*10)/10;
Srange = (0:0.02:1); lenS = length(Srange);
Nrange = (mnn:0.1:mxn); lenN = length(Nrange);
gridlinesS=repmat( Srange*L, lenN, 1);
gridlinesN=repmat( Nrange',1,lenS);
% figure; 
% plot(gridlinesS,gridlinesN, 'k')
% axis equal
% hold on
% plot(gridlinesS',gridlinesN', 'k')
% hold off

[gridlinesX, gridlinesY] = sn2xy(gridlinesS, gridlinesN, ...
    centerlineX, centerlineY);

% Reshape resulting vector back to matrix structure
gridlinesX = reshape(gridlinesX, size(gridlinesS));
gridlinesY = reshape(gridlinesY, size(gridlinesS));

% Plot data and grid in XY and SN coordinates
subplot(3,2,[1 3])
plot(gridlinesX,gridlinesY, 'k')
axis equal
hold all
plot(gridlinesX',gridlinesY', 'k')
colormap('winter')
scatter(dataX, dataY, [], S, 'f')
hold off
subplot(3,2,[5:6])
plot(gridlinesS,gridlinesN, 'k')
hold all; axis equal; plot(gridlinesS',gridlinesN', 'k')
scatter(S, N, [], S, 'f')
clear centerlineX centerlineY dataX dataY gridlinesN gridlinesS ...
      gridlinesX gridlinesY L lenN lenS mnn mxn S N Nrange Srange
      