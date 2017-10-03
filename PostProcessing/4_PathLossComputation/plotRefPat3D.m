function [ hPat3DRef ] = plotRefPat3D( patAz, patEl)
%PLOTREFPAT3D Plot in 3D the reference antenna pattern.
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
% Output:
%   - hPat3DRef
%     The output figure handler.
%
% Yaguang Zhang, Purdue, 10/03/2017

AZS = [patAz.azs; zeros(length(patEl.azs),1)];
ELS = [zeros(length(patAz.els),1); patEl.els];
AMPS = [patAz.amps; patEl.amps];
AMPSDB = 10.*log10(AMPS);

% Shift all the amplitudes in dB to nonegative values.
minAmpDb = min(AMPSDB(:));
AMPSDB = AMPSDB - minAmpDb;

% Convert from the polar coordinate system to the Cartesian system for
% plotting. We have
%    x   = amp * cosd(el) * sind(360-az)
%     y  = amp * cosd(el) * cosd(360-az)
%      z = amp * sind(el)
X = AMPSDB .* cosd(ELS) .* sind(360-AZS);
Y = AMPSDB .* cosd(ELS) .* cosd(360-AZS);
Z = AMPSDB .* sind(ELS);

hPat3DRef = figure('units','normalized', ...
    'outerposition',[0.1 0.05 0.8 0.9]);
colormap jet;
plot3k([X,Y,Z], 'ColorData', AMPSDB);
xlabel('x (to antenna''s right-hand side)');
ylabel('y (to front)');
zlabel('z (to top');
title({'Reference Antenna 3D Radiation Pattern'; ...
    '(Amplitude in dB Relative to the Minimum Value)'});
axis equal; view(135,30);
grid on;

end
% EOF