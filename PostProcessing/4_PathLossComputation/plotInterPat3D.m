function [ hInterPat3D ] = plotInterPat3D( patAz, patEl, ...
    INTER_METHOD, FLAG_INTER_IN_DB, numPtsPerDim)
%PLOTINTERPAT3D Plot in 3D the antenna pattern interpolation results for a
%specified method.
%
% Inputs:
%   - patAz, patEl
%     The antenna patterns, for the Azimuth and Elevation sweeps,
%     respectively; Each of which is a struct containing fields:
%       - azs
%         The azimuth angles in degree from set [0, 360).
%       - els
%         The elevation angles in degree from set [0, 360).
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%     All of these fields contains a column vector with each row
%     corresponding to a sweep sample.
%   - INTER_METHOD
%     A string for the interpolation method to use. This should be
%     supported by antPatInter.m.
%   - FLAG_INTER_IN_DB
%     Set this to be true if the interpolation should be carried out using
%     the amplitude in dB (instead of using the raw linear amplitudes).
%   - numPtsPerDim
%     A scalar specifying how many points should be generated for each
%     angle dimension.
% Output:
%   - hInterPat3D
%     The output figure handler.
%
% Yaguang Zhang, Purdue, 10/03/2017

if nargin < 4
    % By default, we will interpolate directly using the linear amplitudes.
    FLAG_INTER_IN_DB = false;
end
if nargin < 5
    % For illustration, by default we will use 100 points for each axis.
    numPtsPerDim = 100;
end

[azs, els] = deal(linspace(0,360,numPtsPerDim));
[azs, els] = meshgrid(azs, els);

% We need to carry out an interpolation accordingly to the antenna S21
% measurement results.
if FLAG_INTER_IN_DB
    patAzDb = patAz;
    patAzDb.amps = 10.*log10(patAzDb.amps);
    patElDb = patEl;
    patElDb.amps = 10.*log10(patElDb.amps);
    ampsDb = antPatInter(patAzDb, patElDb, azs, els, INTER_METHOD);
else
    amps = antPatInter(patAz, patEl, azs, els, INTER_METHOD);
    % Plot in dB.
    ampsDb = 10.*log10(amps);
end

% We will also plot the sweep data, just for reference.
AZS = [patAz.azs; zeros(length(patEl.azs),1)];
ELS = [zeros(length(patAz.els),1); patEl.els];
AMPS = [patAz.amps; patEl.amps];
AMPSDB = 10.*log10(AMPS);

% Shift all the amplitudes in dB to nonegative values.
minAmpDb = min([ampsDb(:);AMPSDB(:)]);
ampsDb = ampsDb - minAmpDb;
AMPSDB = AMPSDB - minAmpDb;

% Convert from the polar coordinate system to the Cartesian system for
% plotting. We have
%    x   = amp * cosd(el) * sind(360-az)
%     y  = amp * cosd(el) * cosd(360-az)
%      z = amp * sind(el)
x = ampsDb .* cosd(els) .* sind(360-azs);
y = ampsDb .* cosd(els) .* cosd(360-azs);
z = ampsDb .* sind(els);

X = AMPSDB .* cosd(ELS) .* sind(360-AZS);
Y = AMPSDB .* cosd(ELS) .* cosd(360-AZS);
Z = AMPSDB .* sind(ELS);

hInterPat3D = figure('units','normalized', ...
    'outerposition',[0.1 0.05 0.8 0.9]); 
colormap jet; hold on;
plot3k([x(:),y(:),z(:)], 'ColorData', ampsDb(:));
hRef = plot3(X,Y,Z, '.-.', 'Color', ones(1,3).*0.7);
hold off;
xlabel('x (to antenna''s right-hand side)');
ylabel('y (to front)');
zlabel('z (to top');
title({'Interpolated Antenna 3D Radiation Pattern'; ...
    '(Amplitude in dB Relative to the Minimum Value)'});
legend(hRef,'Ref', 'Location', 'southeast');
axis equal; view(135,30);
grid on;

end
% EOF